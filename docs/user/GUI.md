# SASSI-EDU graphical user interface (`sassi-gui`)

The SASSI-EDU GUI is the main window of requirements section 5 (menu bar, toolbars, tabbed document
area, Command Entry, status bar, Options dialogs, plots, help) as a **local web application** shown
in your browser. It replaces the Tk GUI of the original decision D-GEN-05, because the system Tk on
macOS is 8.5.9 and deprecated (ARCHITECTURE section 9).

Everything runs on your own machine. The server listens on `127.0.0.1` only, and Plotly.js comes from
the `plotly` Python package you have installed, so the GUI also works offline.

**One interpreter, one language (UI-01, rule L17).** Each menu item, dialog OK and toolbar button
sends **command text** to the same `sassi.prep.Interpreter` that runs the console (`sassi`) and
`.pre` files. The Command History therefore records a replayable session. A GUI session gives the
same model and the same results as the matching `.pre` file. The command each item submits is given
below with the item.

**New to SASSI?** The GUI opens on the **Learn** tab: a guided course for structural engineers who
design nuclear facilities but have not used SASSI, the tutorial examples, and a command explainer
([section 13](#13-learn-the-guided-course)).

---

## 1. Starting and stopping

```bash
# from the project root, with the project virtual environment
.venv/bin/sassi-gui                                   # opens http://127.0.0.1:8765/ in your browser
.venv/bin/sassi-gui --model-dir ~/ssi/projects --port 8800
.venv/bin/python -m sassi.ui.server --no-browser      # the same program as a module
```

| Option | Meaning |
|---|---|
| `--port N` | Port on 127.0.0.1 (default 8765; `0` picks a free port). If the port is busy, a free port is used and its URL is printed. |
| `--no-browser` | Do not open the browser; open the printed URL yourself. |
| `--model-dir DIR` | Initial working directory. It stays a readable folder for the file dialogs during the whole session. |
| `--settings-dir DIR` | Folder for `SASSIini.xml`, `SASSIdb.xml` and `SASSIani.xml`. The default is the per-user folder: `~/Library/Application Support/SASSI-EDU` on macOS, `%APPDATA%\SASSI-EDU` on Windows, `~/.config/sassi-edu` on Linux. The environment variable `SASSI_EDU_SETTINGS_DIR` also sets it. |
| `--verbose` | Log every HTTP request on the console. |

To stop the GUI, use **Model > Exit** or press Ctrl-C in the terminal. Both save `SASSIini.xml`
and stop a module run that is still going (its worker process and any module it started).
Models and plots are not saved automatically: use **Model > Save** (`SAVE`) first.

You can reload the browser page at any time. The models, plots and Command History live in the
server, and the page restores the open plot tabs.

The same GUI also runs with no installation at all, entirely in a web browser: see
[section 16](#16-web-version-github-pages).

---

## 2. The main window (requirements 5.2)

```
+--------------------------------------------------------------------------------------+
| SASSI-EDU User Interface (ACS SASSI V3 methodology)          Model 0: demo - /path    |
| Model  File  Plot  Modules  Options  View  Learn  Help                               |
| [main toolbar: new open save output converters capture plots animations options .. ]|
| [3D Plot toolbar: view centre wireframe shrink G M P labels DOF mass pause debug]   |
| Command History | Learn | Lesson | SITE out  || Model 0 - Model Plot | Results | Lesson    |
|                                           ||                          | panel     |
|   work: history, lessons, editors,        ||  plots and the Results   | (Learn,   |
|   module output                           ||  browser (split view)    |  docked)  |
| live argument hint: NGEN,itim,[step],n1,... (Learn)                                  |
| Command Entry [ N,1,0,0,0                                                     ] [?] |
| status text                     Model 0 | 30 nodes, 14 elements | cwd ...  [====]   |
+--------------------------------------------------------------------------------------+
```

* **Command Entry.** Type a command and press Enter. Up and Down recall earlier commands.
  Commands are executed in order, one at a time. The **?** button next to it switches the live
  argument hint on and off (section 13.4).
* **Command History** (the first tab) shows the six message classes in their colours (5.6):
  command echo (black), output confirmation (blue), comments (green), information (purple, which
  cannot be hidden), warnings (amber) and errors (red). It follows new messages while it is
  scrolled to its end, also while another tab is in front. **Copy** copies its text.
* **Status bar.** It shows the progress of `INP` (current line / total lines), of module runs
  (for example `ANALYS: frequency 14/22 (7.007 Hz)`) and of `PROCFRAME`. A **Cancel** button appears
  while a module is running.
* **Tabs and the split view.** Plots open as tabs named like `Model 0 - Model Plot` or `Spectrum Plot -
  <title>`. Bringing a plot tab to the front sends `ACTIVATEPLOT,<id>`, and closing it sends `CLOSEPLOT`,
  because the plot setting commands act on the active plot. File editors (`File Editor - <path>`),
  module outputs (`SITE output`), Results, Help and Verification also open as tabs. With **View > Split
  View** (on by default) the plots and the Results browser open in a pane of their own to the right of
  the work -- Command History, Learn, the lesson, File Editors, module output -- so a lesson step or a
  listing stays on the screen beside its plots. Each pane has its own tabs; the right pane closes with its
  last tab; drag the divider to share the width (double-click: equal halves; kept by the browser). The
  keyboard rotation of a 3D plot (section 4) acts after a click in the plot pane, so that the same keys
  still scroll the Command History or a lesson. A window narrower than about 760 px puts the plots below
  the work. With Split View off every tab shares one tab bar, as in ACS SASSI. The selected tab is
  scrolled into view when the tab bar is full, and a plot follows the size of its pane.
* **Dialogs** open over the window. They fit windows down to about 760 px wide: two-column parts
  fall back to one column, long rows wrap, wide tables scroll inside their box and a tall dialog
  scrolls. Escape or the close box cancels; drag a dialog by its title bar. Questions (remove a
  group, unsaved file, Exit) are dialogs of the page, not browser pop-ups.

---

## 3. Menus (requirements 5.3)

A command that needs arguments but is typed without them opens its dialog. For example, `SPECPLOT`
opens Line Selection, `CUTPLOT` opens the cut dialog and `SHOWDOF` opens Boundary Conditions. Inside
a `.pre` file the same commands fail, as in batch mode (D-UI-08).

### 3.1 Model

| Item | What it does | Command text |
|---|---|---|
| New | makes the lowest unused model number the active model | `ACTM,<n>` |
| Open (Ctrl+O) | **Load Model** dialog: the `SASSIdb.xml` tree of groups and models. **Add Group** asks for a name; **Add Model** asks for name, folder and title; **Remove Group / Remove Model** ask first whether to remove the entry and then whether to delete the model files from disk (only files named after the model, `<name>.*`, `<name>_*`). **Open** (or a double-click) loads the highlighted model into the active model. | `MDL,<name>,<path>` + `RESUME` (only `MDL` and `TIT` for a model never saved) |
| Save | saves the active model to `<model>.sdb` | `SAVE` |
| Input | file dialog (`.pre`, `.txt`, `.mac`; *show all files* for others), then runs the file | `INP,<file>` |
| Open Example... | the examples gallery as a dialog (section 13.5): copy a tutorial example into a fresh workspace, then open it in the File Editor or run it | `CD,<workspace>`, `DMODEL,0` (+ `INP,<example>.pre` for Run all) |
| Converters > SASSI .hou | **SASSI .hou to .pre Converter**: Input File Name (`<<` browses), Output .pre File Name (optional: also writes the converted model as commands), Save Converted Data to Model Number (blank = the active model; the dialog proposes the lowest unused number), a note that the converters have had limited testing. **Convert** submits. | `CONVERT,SSI,<model>,<file.hou>,[<prefile>]` (reads the `.sit` / `.poi` decks alongside) |
| Converters > ANSYS .cdb | **ANSYS .cdb to .pre Converter**: as above plus **Enter Value for Gravity** (required: 32.2 ft/s², 386.4 in/s² or 9.81 m/s²) and an optional damping ratio (SASSI-EDU extension; blank = the materials' damping). See [ANSYS.md](ANSYS.md). | `CONVERT,ANSYS,<model>,<file.cdb>,<gravity>,[<prefile>],[<damp>]` |
| Converters > GT-STRUDL Database | greyed (not included) | |
| Output | **Output (WRITE)**: file name and folder of the `.pre` file. Options > Write adds the `MDL` / `AFWRITE` lines. | `WRITE,<file>,<path>` |
| Export to ANSYS | writes the active model as ANSYS APDL input (`<model>.inp` in the model folder) and downloads it to this computer (in the web version too); other names or folders: type `ANSYS,<file>,<dir>` (run `ANSYSREFORMAT` first for beam end releases, manual WARNING) | `ANSYS` |
| Export to STRUDL | greyed (not included) | |
| Exit | lists the models changed since their last SAVE and asks; then saves `SASSIini.xml`, cancels a running module and stops the server | |

Paths and texts that contain a comma are put in double quotes (`MDL,m,"/Users/me/Smith, J"`), so
that the folder stays one argument when the history is replayed.

### 3.2 File

| Item | What it does | Command text |
|---|---|---|
| Open | text editor tab for any text file of the readable folders; a new name creates the file (spec 04 section 4.5). The tab has **File > Save** (also Ctrl+S / Cmd+S), **Input > Connect to Command Entry** and **Run (INP)**. While an editor is connected, every accepted command (typed, or sent by a menu, dialog or toolbar) is appended to it; one editor at a time. **Run (INP)** saves the buffer and runs it with INP; that INP line is not appended to the file itself. Closing an editor with unsaved changes asks first. | `INP,<file>` for Run |
| Export Image | the active plot as an image file, rendered by the headless renderer (D-UI-10): PNG when the name contains `.png`, otherwise BMP | `CAPTUREPLOT,<file>` |
| Export Table | the active spectrum, time-history or soil-property plot as CSV: the x column and one column per line on the union grid (D-UI-17) | (no command; file written by the server) |
| Results Browser | the model folder with result types (TF, RS, histories, decks, listings); **Plot** opens Line Selection for the file, **Open** opens it in an editor (extension) | |
| Upload to Workspace..., Download..., Clear Saved Files... | web version only (section 16): copy text files in and out of the browser's workspace; delete the workspace the browser keeps between visits | |

### 3.3 Plot

| Item | What it does | Command text |
|---|---|---|
| Model > Elements | element plot of the active model (section 6) | `MODELPLOT` |
| Model > Nodes | node plot of the active model | `NODEPLOT` |
| Cuts | **Select Cut to Display**: Cut Number (the defined cuts and their element counts are listed) and Model Number (default the active model) | `CUTPLOT,<cut>,<model>` |
| Spectrum TFU-TFI | **Line Selection** dialog (section 6.1) | `READSPEC` + `SPECPLOT` (+ `PLOTTITLE`, `XTITLE`, `YTITLE`, `AXES`, `PLOTRANGE`) |
| Time History | **Line Selection** for histories | `READTH` + `THPLOT` (+ settings) |
| Soil Layers | layer column and property table of the active model | `LAYERPLOT` |
| Soil Properties | **Select Dynamic Soil Property**: the `DYNP` properties of the model and the built-in curves Clay, Sand and Rock (marked "built-in": SOIL uses them by their label without `INP` while the model does not define that label; editing one stores model `DYNP` points of that label); **New** asks for a name; the 11-point table (strain %, G/Gmax, strain %, damping %) can be edited. Ok submits the edited points and plots. | `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>` ... + `SOILPROPPLOT,<label>` |
| Non Uniform Soil Field | greyed (not active in V3) | |
| Process Animation Frame List | **Parse Frame Data** (section 6.3): List File Name (`<<` a list file, **Folder** a frame folder written by a module restart option), Frame Storage Dir, Data Description, Plot Type | `PROCFRAME,<list>,<dir>,<description>,<type>` |
| Bubble / Vector / Contour / Deformed Shape | **Load Frame Data**: the processed animations of `SASSIani.xml` (Remove Animation asks, then offers to delete the frame store), Start, End, Stride, and the colour-map range and column (Bubble, Contour) or the scale factor (Vector, Deformed) | `BUBBLEPLOT` / `CONTOURPLOT,<dir>,<start>,<end>,<stride>,<min>,<max>,<col>`; `VECTORPLOT` / `DEFORMPLOT,<dir>,<start>,<end>,<stride>,<scale>` |

### 3.4 Modules

| Item | What it does | Command text |
|---|---|---|
| Location | **Module Directories**: per module `built-in` (the Python module of SASSI-EDU, default) or the path of an external executable (refused when the file does not exist). LIQUEF and PINT are greyed. Ok saves `SASSIini.xml`. | |
| Extension | **File Extension Options**: per module the input and output file extensions. Ok saves the edits of all modules; they are passed to external executables (the built-in modules use the manual's deck names). | |
| EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE, ANALYS, COMBIN, MOTION, STRESS, RELDISP | runs the module on the active model in a worker process (section 7), with its listing streamed into an output tab and the Command History, progress in the status bar and **Cancel** | `RUN<MODULE>` |
| LIQUEF, PINT | greyed (not in V3) | |
| NONLINEAR | Option NON (tier P2): runs like the other modules (worker process, listing in the `NONLINEAR output` tab, progress, Cancel). `RUNNONLINEAR` writes `<model>.eql` from the NONLINEAR tab of Options > Analysis (section 5.2) and needs `<model>.hou` and the RELDISP `.THD` files of the panel corners and spring ends. It is greyed only while `RUNNONLINEAR` is a placeholder. | `RUNNONLINEAR` |
| ANSYS Eq. Static Load | **ANSYS Static Load Converter** dialog (Option A, section 7.1). **Ok** submits the changed settings; **Run** submits them and then runs LOADGEN like a module (`LOADGEN STATIC output` tab) | `LOADGEN`, `LGFILE`, `LGNODE`, `LGTIME`, `LGMAP`, `LGOPT` (changed ones); Run: + `RUNLOADGEN,STATIC` |
| ANSYS Dynamic Load | **ANSYS Dynamic Load Converter** dialog (section 7.1), Ok / Run as above | `LOADGENDYN`, `LGFILE`, `LGNODE`, `LGMAP`, `LGOPT`; Run: + `RUNLOADGEN,DYNAMIC` |
| ANSYS Super Element Utilities | greyed (tier P2) | |

### 3.5 Options

| Item | What it does | Command text |
|---|---|---|
| Model | **Model Options** (section 5) | `MOPT,<incomp>,<matrix>,<mass>,<force>` |
| Write | **Extended Write Options**: MDL command in *.Pre, extend integer fields (recorded, no effect on the decks), AFWR command in *.Pre, simulation commands to `*-Sim.Pre` or `*.Pre`. Session settings: they change what WRITE writes. | (no command in ACS SASSI) |
| Check | **Check Options**: show warnings / errors, Suppress Error Window, Break Check at N messages. Session only, not saved (UI-07). | (no command) |
| Analysis | **Analysis Options**, the 12 tabs (section 5) | the record, list and request commands of the changed fields |
| Windows Settings | the settings window of the active plot: Graph Plot Options (spectrum, time history), Soil Layer Windows Setting, Soil Properties Window Settings, Window Options (3D plots and animations); section 6 | `PLOTTITLE`, `AXES`, `PLOTRANGE`, `WINDOWSETTINGS,...` |
| Colors | the six Command History colours (saved in `SASSIini.xml`). Plot palettes are set with `COLOR,<Palette>,<Num>,<R>,<G>,<B>`. | |
| Font | greyed (tier P2) | |
| Shader Options | node / bubble size, vector and displacement scale, element outline thickness and shrink (fractions: 0.06 = 6 %); saved in `SASSIini.xml` | `SHADEROPTIONS,<points>,<linew>,<shrink>,<scale>` |
| Reset Plot | resets the view and the rotation centre of the active 3D plot | `RSTVIEW` + `RSTCENTER` |

### 3.6 View

| Item | What it does |
|---|---|
| Check Errors | the **CHECK: Errors and Warning for - &lt;model&gt;** window showing `<model>.err` (requirements 5.5): module headers in run order, errors red, warnings amber, the totals line. **Run CHECK** submits `CHECK` and refreshes it. The window also opens by itself after a CHECK or AFWRITE that reports messages, unless Options > Check > Suppress Error Window is set. |
| Command Window | brings the Command History tab back and clears it |
| Command Display > Command Echo / Output Confirmation / Comments / Warnings & Errors | show or hide these message classes in the Command History (information is always shown); saved in `SASSIini.xml` |
| Toolbars > Main Toolbar / Plot Toolbar | show or hide the toolbars (not saved, UI-07) |
| Split View | plots and the Results browser in a pane of their own beside the work (on by default; section 2); kept by the browser |
| Run Summary | the summary window of the active model's run (section 7.2) |
| Results Browser | as File > Results Browser |

### 3.7 Learn

The guided course (section 13): **Start Page** (the Learn tab), **Continue: <lesson>**, the lessons by
part (Fundamentals, Design applications, Advanced) with the steps completed (`3/7`), **How SASSI Works**,
**Explainer Videos** (section 13.7), **Examples Gallery**, **Open Example...**, **Explain a Command...**, **Lesson Layout** (full width or side
panel), **Course Workspace Folder...**, **Free Disk Space...** and **Reset Course Progress**.

### 3.8 Help

| Item | What it does |
|---|---|
| Help (F1) | the **Help** tab (section 8): the SASSI-EDU documentation rendered from the Markdown files in `docs/` and `examples/README.md`, with links, contents, search and the live command index |
| About | version, build, methodology, and the Python, NumPy, SciPy, Plotly and Matplotlib versions, platform and settings file |
| Verification | the **Verification** tab: runs `VERIFY` in a worker process and shows the result table (section 9) |

---

## 4. Toolbars (requirements 5.9)

Every button submits the same command text as the matching menu item, or opens the same dialog.
Toolbar visibility is not saved (UI-07).

**Main toolbar:** New (`ACTM`), Open the Model Database (Load Model), Save (`SAVE`), Output model to
.pre (WRITE dialog), Output model to ANSYS (`ANSYS` dialog), .hou converter, ANSYS converter, GT
STRUDL converter (greyed), Capture (Export Image, enabled when a plot is open), Element plot
(`MODELPLOT`), Node plot (`NODEPLOT`), Cut plot, Spectrum (Line Selection), Time History (Line
Selection), Soil Layer (`LAYERPLOT`), Soil Properties, Process Frames (Parse Frame Data), Bubble,
Vector, Contour, Deformed (Load Frame Data), Analysis Options, AFWRITE (`AFWRITE`), and the
extensions Check Errors window and Results browser, and **Learn** (the Learn tab, section 13).

**3D Plot toolbar** (enabled for the plots each command applies to; a lit button shows a setting that
is on): Change View (dialog: rotations about X, Y, Z in degrees, pan X / Y, zoom 1 = fit:
`CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>`), Reset View (`RSTVIEW`), Change Centre (dialog:
`CNGCENTER,<X>,<Y>,<Z>`), Reset Centre (`RSTCENTER`), Wireframe (`WIREFRAME`), Shrink (`SHRINK`),
colours by Group / Material / Property (`ELECOLOR,1|2|3`), Node / Element / Group labels (`NODENUM`,
`ELENUM`, `GROUPNUM`), Show DOF (Boundary Conditions dialog: `SHOWDOF,<X ... ZZ>` or `SHOWDOF,NONE`),
Show Mass (`SHOWMASS`), Show Soil (`SHOWSOIL`: the free-field soil layers drawn around the foundation, a display
aid; [User Guide §11.4](USER_GUIDE.md#114-plots)), Pause / Start (animations, section 6.3) and Debug (`DEBUG`: the view values on
the plot). The keyboard rotates the active 3D plot by 5 degrees per press: Insert / Delete (X),
Home / End (Y), PageUp / PageDown (Z), each as a `CNGVIEW`.

---

## 5. Options dialogs (requirements 5.4)

**Options > Model** submits `MOPT,<incomp>,<matrix>,<mass>,<force>`.
**Options > Write** (Extended Write Options) and **Options > Check** have no command in ACS SASSI. They
set the session settings directly. Check Options are not saved (UI-07).

**Options > Analysis** is one window with the tabs `EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE |
ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`. A changed tab is marked with `*`. Every
field of the requirements 5.4 table is on its tab and is bound to the one command that stores it;
OK submits only what changed, and reopening the dialog shows the stored values (WRITE writes the same
commands):

* Record fields (`SITE <nl>`, `HOUSE <coh>` ...). OK submits the complete command line, for example
  `SITE,0,1,0,12,3,1,0,1,,1,0,0.01,1024,1` (a blank Frequency 2 stays blank: NFFT/2).
* Indexed entries. A selector picks the key: Spectrum Number for `RSIN/RSOUT/ACCIN/ACCOUT/TPSD`,
  Layer Number for `SPRO/SACC/SRS/SSTR/SSAF`, the wave page for `WAVE`, Input Motion Number for `ME`
  and `AMP`.
* Lists (`DAMP`, `TOPL`, `AMP`). OK submits the clear form first, then the values, e.g. `TOPL,0` and
  `TOPL,1,2,3`. Separators are blank, tab, `,`, `;` or Enter.
* Request tables (`NOUT`, `EOUT`, `RDND`). They have Add/Delete rows; OK submits `NOUT,0` and then
  one command per row. Node and element lists accept `1, 3-6 10`. A wide table scrolls inside its box.
* String fields (`THFILE`, `THTIT`, `EQTIT`, `RELFILE`). **<<** browses; **Edit** opens the file (a built-in
  `@` file opens read-only). File fields that can take a built-in input (THFILE, the SOIL-only file, RSIN,
  ACCIN, TPSD) also have a **Library** button: it lists the built-in files that fit the field (records, the
  Ricker load pulse, spectra, target PSDs; command `LIBRARY`) and puts the `@` name into the field, which
  then sits on its own line below the buttons. A blank file with a built-in default shows that default as
  its placeholder, for example `built-in: RG 1.60, 0.30 g (@rg160h_030g.acc)`, `built-in: 5 Hz Ricker pulse
  (@ricker_5hz.th)` for Foundation Vibration, `default: <model>_eq1.rso` for RSOUT; the placeholder follows
  the values being edited (type of analysis, File Contains Pairs, time step, spectrum number) and says when
  no default applies. CHECK reports every default used as Warning EDU-29 (User Guide section 4.6).
* Records of the later options (Option NON `EQL`, `P`, `S`, `BBCI`/`BBCX`/`BBCY`; SOIL-NON `NLSOIL`,
  `NLSLAYER`; sections 5.1 and 5.2). Their commands check the values themselves, so OK first runs the
  commands on a copy of the model; an error refuses the whole commit with the command's own message,
  for example `EQL: <disp> = 2: the equivalent-linear displacement factor must lie in (0, 1]`. Tables of
  such entries (panels, springs, nonlinear soil layers) submit a deletion for a removed row (`PDEL,3`,
  `DELSPR,3`, `DELNLS,3,3,1`) and the complete command for a new or changed row; unchanged rows submit
  nothing.

**Shared variables are one value** (spec 07 section 3). The time step, NFFT, frequency step, frequency
set and control layer belong to SITE; the SSI gravity and ground elevation to HOUSE; the type of
analysis to ANALYS; the control-motion scaling and records to MOTION. A shared variable shown on
several tabs is the same field: edit it on one tab and every tab shows the new value.

Defaults are the new-model defaults of the 5.4 "Default" column (D-UI-03). A field left blank keeps
its documented meaning; for example a blank EQUAKE Number of Frequencies means "the number of
records of RSIN 1", and a blank SITE **Frequency 2** means NFFT/2 (what AFWRITE writes), so it follows
the Nr. of Fourier Components. The SITE frequency numbers show their value in Hz next to the field,
computed as `n * df` with `df = 1/(dt * NFFT)` or the frequency step.

**The dialog shows what the analysis uses.** Some pages are not stored until you change them:

* **SITE wave pages.** Without any `WAVE` command, AFWRITE writes the vertical SV field (SH for
  "SH- and L-Waves"), so that page shows the field on and the others off. When you switch on
  another page, OK stores the SV (or SH) field too, because once any `WAVE` entry exists AFWRITE
  writes only the stored ones. Once a `WAVE` entry exists, every page without an entry shows "No
  Wave Field".
* **SOIL layer output** (Accelerations, Response Spectrum, Stresses/Strains, Spectral Amplification).
  A layer without a request shows "not requested", because AFWRITE writes no row for it and SOIL
  computes nothing there. The first change on such a layer starts from the 5.4 values (Compute
  Maximum; Save RS; compute stresses and strains), and OK stores the request (`SACC,3,1,0` ...).
* **Stochastically Simulated Incoherency Input** (HOUSE tab). The input is stochastic only when a
  Horizontal or Vertical SEED is non-zero and the Random Phase Angle is above 0 (D-INC-02). Choosing
  it sets the Random Phase Angle to 180 (5.4 default). OK refuses the choice with a message while
  both SEEDs are 0. Choosing Deterministic (Median) sets the SEEDs and the phase to 0.

**Validation (UI-06).** OK checks each record with the same CHECK rules as the commands. An invalid
value refuses the whole commit and shows the Chapter 10 message, for example
`SITE: Error 47 : Illegal Number of Layers for Halfspace Simulation [<nl> = 3 (0 or 4..20)]`. Fix the
value and press OK again. The check covers the whole record, so a record that was incomplete before
(SOIL without Number of Values, POINT without a radius) must be completed when you edit it.

**Enable/disable rules (D-UI-04)**

| Situation | Effect |
|---|---|
| Coherent motion | the incoherency fields are disabled |
| 2D dimension | Incoherent is disabled; OK stores Coherent |
| Unlagged coherency model 2-7, or Multiple Excitation on | Use Wave Passage is switched on and locked, with a note. It returns to its previous value when the reason goes away. |
| Coherent | Free-Field Load / Free-Field Motion is disabled (enabled only when incoherent) |
| Interpolation option 6 | the smoothing parameter is disabled |
| Type of analysis | the vibration response type is shown only for Foundation Vibration |
| Ansys Model Input off | Ansys Model Type is greyed (P2) |
| SOIL-NON Damping Type not 3 (Rayleigh) | the Mass and Stiff Matrix Mult. are disabled |
| LOADGEN Data to Add | Acceleration disables the displacement boxes, Displacement / Disp. for Soil Module the acceleration and mass boxes; a rotational box needs its check box; the mass file box follows the Mass Type |
| LOADGENDYN damping ratio zeta > 0 | Alpha and Beta are disabled (computed from zeta at f1 and f2) |

The AFWRITE tab sets `AOPT`; LIQUEF and PINT are greyed.

Helpers: **From mesh (RADIUS)** on the POINT tab runs `RADIUS` and copies the average central-zone
radius into the field. The **Edit** buttons open `<model>.pin`, `SRSSTF.txt`, `CONTTRS.txt` or
`Frames.txt` in a File Editor tab. The **...** button next to the SOIL Dynamic Soil Property opens
Select Dynamic Soil Property. The field's list proposes the model's DYNP labels and the built-in curves.
While the model has no `SPRO` entry, every SOIL layer page shows the default profile AFWRITE writes (the
TOPL layers with Sand or Rock by Vs, the SITE half-space last; the "Without SPRO" line summarises it);
changing one layer stores the whole profile with it, as the implicit wave field is stored with the first
WAVE page.

### 5.1 SOIL tab: Nonlinear Soil Behavior (SOIL-NON)

The group of spec 05a section 6.6, bound to `NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,
<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>`: **Nonlinear Time Domain** (`<Opt>`), Subincrements
per Timestep, Displacement and Force Convergence Error, Equilibrium Iterations, **Bedrock Interface**
Rigid / Viscoelastic (the manual calls Viscoelastic "not applicable to this version"; SASSI-EDU
implements it), **Damping Type** 1 frequency independent, 2 visco-elastic, 3 Rayleigh (0, the default,
is used as 1), and the Rayleigh **Mass / Stiff Matrix Mult.** (enabled for type 3). Every field is
written when one changes, e.g. `NLSOIL,1,0,0,0,0,0,3,0.1,0`. A model without `NLSOIL` shows the
defaults (all 0: SOIL runs SOIL-EQL) and OK without a change submits nothing.

**Nonlinear Soil Layers (NLSLAYER)** is a table below the columns: one row per SOIL sublayer (SPRO
numbering, 1 = top) with Layer, Curve Fit (fit to the sublayer's DYNP G/Gmax curve), Beta, S exponent,
Reference Strain (%) and Viscosity. **Add the SOIL sublayers without a row** adds a Curve Fit row for every
SPRO sublayer that has none (the half-space too: CHECK Error 125 counts it). OK submits `NLSLAYER,<n>,...`
for new and changed rows and `DELNLS,<n>,<n>,1` for removed ones. AFWRITE writes both commands to
`<model>.nls`.

### 5.2 NONLINEAR tab (Option NON)

The tab of spec 05d section 3.8:

* **Global Modeling Options** = `EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>`: Disp.
  Factor (EDF), Damping Cutoff %, Damping Scale Factor (0 = 1), Use Non-linear Panels / Springs (bits 1
  and 2 of `<NonLinOpts>`), Include Elastic Damping. Material Parameter and Use Non-linear Beams are
  greyed, as in the manual. While `<NonLinOpts>` was never given, the check boxes show the element types
  of the P and S records (what RUNNONLINEAR uses) and OK keeps the argument blank (`EQL,0.75,,7,0,1`);
  a click on a type stores it.
* **Backbone Curve Data**: the Backbone Curve number (spin box; **Defined curves** lists the curves of
  the model), Type (1 CMS, 2 CMB, 3 TAK, 4 GMR: titles only), Yield Num., **New Curve** (the next free
  number) and **Delete Curve**. **Backbone Curve Points** is the X / Y grid of the curve (origin omitted,
  point 1 = cracking point) with Add / Delete rows. OK redefines a changed curve completely
  (`DELBBC,1`, `BBCI,1,<yield>,<type>`, `BBCX,1,<n>,<yield>,...`, `BBCY,...`), so the stored curve is the
  one shown and no yield-number warning appears; a deleted curve is `DELBBC,<n>`. BBC (file), BBCGEN and
  BBCP still define curves from the Command Entry.
* **Panel Data** (`P,<num>,<group>,<bbc>,<disp>,<force>`) and **Spring Data**
  (`S,<num>,<group>,<elem>,<bbc>,<disp>,<force>`) are tables below the columns: Disp. Opt. 1 shear strain /
  2 bending (experimental), Force Opt. 1 CMS / 2 CMB (not included) / 3 TAK (experimental); spring Dof.
  1 X / 2 Y / 3 Z and Force Opt. 4 GMR. Removed rows are `PDEL,<n>` / `DELSPR,<n>`. The commands warn
  about a group that is not a SHELL (SPRING) group; PLIST checks the panels.
* **Beam Data** is greyed (not available in this version); it shows the number of stored B records.

The panel model itself is built with the commands (PNLGEN, WALLFLR, PANELIZE, EDGE ...; SHEAR and
BBCGEN for the curves; Help > Option NON). **Modules > NONLINEAR** runs `RUNNONLINEAR`, which writes
`<model>.eql` from this tab.

---

## 6. Plots

| Plot | Drawn as | Interaction |
|---|---|---|
| Element (MODELPLOT) | Plotly `mesh3d` faces for SOLID / SHELL / PLANE, lines for BEAMS, SPRING and GENERAL, element colours by group / material / property (ElemPalette), red interaction nodes, mass markers (X red, Y green, Z blue) | toolbar toggles; mouse: drag to rotate (turntable about the vertical Z axis), wheel to zoom, right-drag to pan; a click on a node or face prints its number in the Command History |
| Node (NODEPLOT) | element-connected nodes: black, interaction red, fixed DOF green square, mass, selected blue square | `NODESEL`, Window Options hide/show |
| Cut (CUTPLOT) | wireframe of the model with the cut elements filled red | any model number |
| Spectrum / Time History | one line per line object (SpecLines / THLines colours), legend, log axes, minor grids, markers, stippling; the axes span the data extent unless PLOTRANGE sets them | Graph Plot Options (Windows Settings) |
| Soil Layer | layer column proportional to thickness plus the property table, with a Halfspace row | Soil Layer Windows Setting (Start / End layer, columns) |
| Soil Properties | damping (left axis, red) and G/Gmax (right axis, green) against shear strain % on a log axis | Soil Properties Window Settings |
| Bubble / Contour / Vector / Deformed | animated frames from `PROCFRAME` stores: jet colour bar, bubble sizes, X/Y/Z vectors, `x' = x + Scale u` | Pause/Start, `−`/`+` frame steps and the slider while paused, Window Options |

The 3D view follows the plot's camera: rotation, zoom and pan of `CNGVIEW` and the rotation centre of
`CNGCENTER`. The view is orthographic, and zoom 1 fits the whole axis box (model, axes and tick labels,
with a 10 % margin) around the rotation centre. With the mouse, a drag rotates the scene as a turntable
about the vertical (global Z) axis, so the model stays upright; a view rolled about the viewing axis or
looking straight along Z rotates freely instead. The wheel zooms. Mouse rotation and zoom are display
actions and are not recorded: they are kept when the window or the panels are resized, when a toggle
re-draws the plot (labels, colours, shrink ...) and between animation frames, and any view command
(Change View, Reset View, Change Centre / `CNGVIEW`, `RSTVIEW`, `CNGCENTER`, `RSTCENTER`) replaces them by
the commanded view.

### 6.1 Line Selection (Plot > Spectrum TFU-TFI / Time History)

The dialog of spec 06 section 5.2:

* **Lines in memory**: one row per line object, `<n>: <name>` with its kind and number of points.
  Check the lines to plot.
* **Title / Axis Labels** and **X / Y Axis Options**: Min, Max (blank = data extent), Logarithmic,
  Show Ticks (minor grid lines).
* **Input Line File**: File Name (`<<` browses; the dialog shows the file's type, data columns and
  points), Starting Number (the next free line number; a number you type is kept), Lines in file
  (spectra: columns after the frequency column; transfer functions `f amp phase` load 1 = the
  amplitude) or Line Number and Pair (histories). **Add Line(s)** submits `READSPEC,<file>,<n>,<start>`
  (or `READTH,<file>,<pair>,<num>`); the new lines are **checked**.
* The dialog **starts from the plot it would extend**: the active plot when it has the same kind,
  otherwise the most recent plot of that kind. Its lines are checked and its title, axis labels, log
  axes, minor ticks and ranges are filled in. (Results > Plot of a file starts empty.)
* **Ok** loads a file still waiting in File Name, then opens a **new** plot with the checked lines (at
  most 50), as in ACS SASSI, and submits the settings that make it look as shown (a new plot starts
  from the session defaults, so only the differences are submitted): `SPECPLOT,1,2,3`,
  `PLOTTITLE,...`, `XTITLE,...`, `AXES,1,1,<minorX>,<minorY>,<logX>,<logY>`,
  `PLOTRANGE,<xmin>,<xmax>,<ymin>,<ymax>`. **Replace plot n** closes the plot the dialog started
  from first (`ACTIVATEPLOT` + `CLOSEPLOT`). Ranges are checked before anything is submitted (min
  below max; a logarithmic axis needs min > 0).

### 6.2 Windows Settings of a plot

* **Graph Plot Options** (spectrum and time-history plots, spec 06 section 5.3): the lines on the plot
  (also the inputs of the calculations), titles and axes as in Line Selection, Data Points of the
  highlighted line (`MARKERS`) and Line Stippling (`STIPPLE`). **Average**, **SRSS**, **Broaden**
  (Peak Difference %, Broaden %), **Linear Combin** (coefficients of the first three checked lines),
  **Addition** and **Subtraction** create a new line (`AVERAGE,<dest>,...` etc.; the Post Processing
  Results line number, blank = the next free number); no file is written (use `WRITESPEC` /
  `WRITETH`). OK re-plots when the checked lines changed (`CLOSEPLOT` + `SPECPLOT,...`) and submits
  the changed settings.
* **Soil Layer Windows Setting**: Start Layer, End Layer (-1 = deepest) and the columns shown
  (`WINDOWSETTINGS,START|END|SHOW,...`).
* **Soil Properties Window Settings**: logarithmic axes and ticks (`AXES`), the modulus and damping
  lines (`WINDOWSETTINGS,SHOW,MODULUS|DAMPING,<0|1>`).
* **Window Options** (3D plots and animations, spec 06 section 1.4): Model Display Volume
  (`WINDOWSETTINGS,VOLUME,...`), Show Element Group (`HIDEGROUP` / `SHOWGROUP`), Colormap Value Range
  (Bubble, Contour), Output Direction X / Y / Z / All (Vector), Hide/Show Elements or Nodes (applied at
  once by the button: `WINDOWSETTINGS,HIDEELEM,<group>,<list>` ...; the requests are stored in the
  model), Scale Factor and Frame Pause, Show Undeformed Shape (Deformed), Title, and **Show all**
  (`WINDOWSETTINGS,SHOWALL`).

### 6.3 Animations (requirements 5.8)

1. Ask a module for frames: Options > Analysis, MOTION tab, Post Processing Options **Restart for TF**
   (complex transfer functions of every node, folder `TFU/`), **Restart for ACC** or **Restart for RS**;
   RELDISP tab **Restart For Frame Generation** (relative displacements, folder `THD/`); STRESS tab
   Restart for Nodal Stress / Soil Pressure Contours. Then AFWRITE and run the module. A shorter
   MOTION *Total Duration to be Plotted* gives fewer history frames.
   **Steady-state motion at one frequency** (SASSI-EDU extension): `HARMFRAME,<FILE8>,<f>,<folder>`
   writes the harmonic motion of every node at the computed frequency closest to *f*,
   u(t) = Re(H e^{iωt}) per unit control motion, as 24 frames over one period (`HARM_<ωt°>_<k>`; also
   from the `TFU/` frames; `[NFrames]`, and `[Ref]` 0 relative to the free field or a node number).
   It is the clearest picture of an SSI mode: the structure deforming while the foundation slides
   and rocks. No module rerun is needed: FILE8 holds every node.
2. **Plot > Process Animation Frame List**: List File Name = the frame folder (**Folder** button) or a
   list file; Frame Storage Dir (created); a description; the plot type. `PROCFRAME` writes the frame
   store (`frame_*.npy`, `index.json`) and the `SASSIani.xml` entry.
3. **Plot > Bubble / Vector / Contour / Deformed Shape**: choose the animation in Load Frame Data, the
   frame range and stride, the colour range (Bubble, Contour) or scale (Vector, Deformed).

**The shortcut for a deformed shape** (SASSI-EDU): after an analysis, **Plot > Deformed Shape** (or the
toolbar's Deformed button) also offers *Or animate the analysis results*: the FILE8-type files of the active
model (FILE8, or FILE8001 ... of the load cases of a vibration analysis), the computed frequencies with the
one of the largest deformation preselected (where the motion differs most from node to node: the SSI
resonance of a structure), and total motion or relative to the free field. **Animate** submits the three
steps above as command text, `HARMFRAME,FILE8,<f>,HARM_<f>`, `PROCFRAME` into `HARM_<f>_ani` and `DEFORMPLOT`
with the automatic scale (the largest displacement drawn as 15 % of the model size), so the Command History
replays it. Without results the dialog says to run the analysis first.

The animation tab has **Pause / Start**, `−` / `+` (one frame, while paused) and a frame slider.
**Pause** keeps the frame on the screen and records it: `PAUSE` + `WINDOWSETTINGS,FRAME,<k>`, so a
replayed session (and `CAPTUREPLOT`) shows the same frame; moving the slider while paused submits
`WINDOWSETTINGS,FRAME,<k>` too. The toolbar's Pause/Start button and the Pause key do the same.

---

## 7. Module runs and Cancel (UI-04)

`Modules > <MODULE>` runs `RUN<MODULE>` in a **worker process**, so the GUI stays responsive. A
`RUN<MODULE>[,model]` typed in Command Entry (full name or abbreviation) runs the same way, with its
own output tab and Cancel (a new run of a module replaces the finished output tab of its previous run;
the listing file keeps the text). Commands sent in the same request after a typed `RUN<MODULE>` wait
until the run has ended, then run in order. The worker interpreter holds a copy of your models, and
a module run only writes files in the model folder. The tier-P2 runs work the same way: Modules >
NONLINEAR (`RUNNONLINEAR`, which first writes `<model>.eql`) and the Run of the two LOADGEN dialogs
(`RUNLOADGEN,STATIC` / `RUNLOADGEN,DYNAMIC`, which first writes `<model>.lgn`; section 7.1). A typed
`RUNLOADGEN,[STATIC|DYNAMIC],[model]` (or `RUNLOADGEN,<model>,[STATIC|DYNAMIC]`) also runs as a job, in
the tab `LOADGEN STATIC output` or `LOADGEN DYNAMIC output`. The run:

* checks its prerequisites as the command does: MDL, the deck from AFWRITE, no CHECK error and the
  upstream files (requirements 2.6); a missing prerequisite ends the run with a message such as
  `RUNEQUAKE: ex01.equ not found -- Run AFWRITE first (with EQUAKE enabled in AOPT)`;
* streams the listing into the `<MODULE> output` tab and the Command History;
* shows the module's progress (for example `ANALYS: frequency 14/22 (7.007 Hz, schur)`) in the status
  bar.

**Cancel** (status bar or output tab) terminates the worker process and anything it started. Only one
module runs at a time. After a successful run, `RUN<MODULE>` is added to the replay history at the
place where the run started, before any command you typed while it was running. Inside a `.pre`
file, macro or FOREACH loop, `RUNxxx` runs synchronously in the main interpreter.

**Modules > Location** can point a module at an external executable instead of the built-in Python
module, for example to cross-check against the commercial program. That executable is run without a
shell, in the model folder, with the three-line batch protocol on its standard input:
`model`, `model.<input ext>`, `model_<module>.out` (the extensions come from Modules > Extension).
The same preconditions apply as for the built-in module: MDL set, the input deck written by AFWRITE,
and no CHECK error for that module in the last AFWRITE.

**The activity panel.** While an input file runs (`INP`: Model > Input, Run (INP) in the File Editor, an
example's *Run all*, a lesson step) or a module runs, a panel at the bottom right shows what is happening:
the file and its progress (line *k* of *n*), the command being executed in plain words (from the command
explainer, e.g. *Generating nodes*, *Writing the module input files*), the module running with its own
progress (*ANALYS: frequency 14/22 (7.0 Hz)*; an animated bar for a module that reports none), one chip
per module run (green with its time, *x n* when it runs several times, as in nonlinear iterations; red
when it failed), the elapsed time and the warning and error counts. A clean run fades out after a few
seconds; with errors the panel stays, and **Command History** jumps to the first error. **x** hides it (the
run goes on). The panel only reads the session's events (it sends no command). Long module listings are
added to the Command History in batches, so the page stays responsive during long runs.

**The computation being done.** Below the module's progress the panel shows the step the module is in and
its mathematics, with the sizes of the model being solved: a strip of the module's steps (ANALYS:
*Impedance*, *Dynamic stiffness*, *Condense + factorise*, *Solve*, repeated at every frequency; MOTION:
*FFT*, *Interpolate*, *Convolve*, *Spectra*; SITE, POINT, HOUSE, FORCE, STRESS, RELDISP, SOIL and EQUAKE
likewise), the equation in the notation of the Theory Manual -- for example
S = C_ff − C_fn C_nn⁻¹ C_nf and LU(S + X_ff) -- and the numbers of this run (*frequency 14/97 · 3.589 Hz ·
sparse LU of C_nn (960 DOFs) · dense LU of S + X (495 × 495)*). **Theory §…** opens that section of the
Theory Manual in Help; **∑** in the panel's title bar hides or shows the equations (kept by the browser).
The modules report the steps as they enter them (a display aid: nothing of the run depends on it).

### 7.1 ANSYS Eq. Static Load and ANSYS Dynamic Load (Option A, LOADGEN)

The two Modules-menu dialogs of spec 04 section 15.6 (manual 6.4.15), with the SASSI-EDU settings of
[OPTION_A.md](OPTION_A.md). The manual's Ok writes the LOADGEN input and starts LOADGEN; here **Ok**
stores the settings as commands (so a session replays from its `.pre` file, L17) and **Run** stores them
and then runs `RUNLOADGEN` as a module job (section 7). A box left blank keeps its default, shown in grey
(`<model>.hou`, the model directory, `<model>_LGS.inp` ...). Reopening a dialog shows the stored values;
the file boxes they share (Path, HOUSE Module Input, Ansys Path) are the same `LGFILE` entries.

**ANSYS Static Load Converter** (ANSYS Eq. Static Load):

| Group | Fields | Command |
|---|---|---|
| Data to Add From ACS SASSI to the ANSYS model | Displacement / Acceleration / Disp. and Accel. / Disp. for Soil Module; Use Multiple File List Inputs | `LOADGEN <data>, <multi>` |
| SSI Model and Results Input | Path, HOUSE Module Input, Displacement Results (+ Rotational Disp. and its box), Trans. Acceleration Results (+ Rotational Accel. and its box), History Source RESULTS / FILE8 | `LGFILE,SSIPATH / HOUSE / DISP / DISPROT / ACC / ACCROT,<name>`; `LOADGEN <rotdisp>, <rotacc>, <source>` |
| Ansys Model and Data Input | Path | `LGFILE,ANSYSPATH,<dir>` |
| Mass Data for Internal Load | Mass Type Lumped / Master Node, Generate Mass Data, Lumped node file, Master Node Mass file, Master Nodes | `LOADGEN <masstype>, <genmass>`; `LGFILE,LUMPED / MASTER`; `LGNODE,M,...` |
| ANSYS Output File | ADPL File | `LGFILE,APDL,<file>` |
| Critical Times | Criterion (V, VX ... MZ, ACC, DISP, TIME, STEP), Number of Peaks, Minimum Separation, Moment Point, Node and DOF, Times / Steps | `LGTIME,...` |
| Nodes, Numbering and Output | D Nodes, ANSYS Numbering IDENTITY / PAIRS / COORD with its file and tolerance, APDL Significant Digits, Data Check Only, Relative Displacements Start at Rest | `LGNODE,D,...`; `LGMAP,...` (IDENTITY: `LGMAP`); `LGOPT,<digits>,<opmode>,<rest>` |

**ANSYS Dynamic Load Converter** (ANSYS Dynamic Load): SASSI Model and Results Input (Path, HOUSE Module
Input, Ground Acceleration File with its format and scale factor), Ansys Model and Data Input (Path),
Raleigh Damping Coeff. (Alpha, Beta; or a damping ratio zeta at f1 and f2, then Alpha and Beta are
computed and written blank), ANSYS Output File (ADPL File = `LGFILE,APDLDYN`), Method (REL / ACC with the
Reference Node, History Source, rotations) and the node lists D and A, numbering and output options.
Commands: `LOADGENDYN,...`, `LGFILE`, `LGNODE`, `LGMAP`, `LGOPT`; Run adds `RUNLOADGEN,DYNAMIC`.

OK runs the commands on a copy of the model first: an invalid value (a criterion ACC without node, a
descending node range, method ACC without reference node ...) refuses the commit with the command's
message. Run needs the SSI results LOADGEN reads (HOUSE, MOTION, RELDISP or FILE8; OPTION_A.md section 3);
a missing input ends the job with LOADGEN's message in its output tab. A changed node list is written as
WRITE writes it: `LGNODE,M,0` (clear), then the nodes.


### 7.2 The run summary window

When a run that ran modules ends -- an input file (`INP`, an example's **Run all**) or a module run from the
Modules menu or Command Entry that produces results (MOTION, STRESS, RELDISP, SOIL, EQUAKE, COMBIN,
NONLINEAR) -- the **Run Summary** window shows the model's key inputs, its key outputs and their graphs:

* **Modules**: one chip per module listing in the model folder, with its status, time and warnings.
* **Key inputs**, from the decks the modules read (`<model>.sit`, `.hou`, `.anl`, `.mot` ...): the model
  (nodes, elements by type, interaction nodes and DOFs, embedment, method, the structure and excavated-soil
  masses of the HOUSE listing), the soil profile (layers, Vs and damping ranges, half-space, waves, control
  point), the frequencies (number and range, NFFT, Δt, Δf), the analysis, the input motion (file, scaling, peak,
  duration) or the load history, the spectrum damping, the output requests; SOIL and EQUAKE inputs for a
  free-field run.
* **Key outputs**, from the result files: the transfer functions (`.TFI`: the largest amplitude and its
  frequency, the SSI resonance), the peak accelerations (`.ACC`, also as a multiple of the input's), the
  in-structure response spectra at 5 % damping when computed (`.RS`: peak, frequency, ZPA), the relative
  displacements (`.THD`) and the largest element forces (`.THS`); the SOIL / EQUAKE spectra and histories.
* **Graphs**: the transfer-function amplitudes, the largest acceleration history with the input over it, the
  response spectra with the input motion's spectrum (dashed, computed for the window), the soil profile.
  **Open as plot** opens a graph's result files as a plot tab (`READSPEC` / `READTH` and `SPECPLOT` / `THPLOT`
  in the Command History).

With several models in memory (an example with a fixed-base reference, or FV and subtraction variants) a
selector switches between them. The window does not open by itself while a lesson is open (the lesson has its
own result buttons) or after a single intermediate module (SITE, POINT, HOUSE, FORCE, ANALYS); **Show this
summary after every run** switches it off (kept by the browser), and **View > Run Summary** opens it at any
time.
---

## 8. Help (F1)

The Help tab shows the documentation of the project **as it is on disk**, rendered to HTML by the
server (a regenerated Verification Manual appears at once). LaTeX formulas (`$...$`, `$$...$$`, fenced
`math` blocks) are typeset with KaTeX:

| Group | Documents |
|---|---|
| Start here | `docs/index.md`, the documentation map (the Home page) |
| User guides | `docs/user/*.md`: the User Guide, this GUI guide, the ANSYS interface and any later guide |
| Theory | `docs/theory/*.md`: the Theory Manual |
| Verification | `docs/verification/*.md`: the Verification Manual (tables and figures) and the impedance study |
| Reference | `docs/reference/*.md`: the Command Reference |
| Examples | `examples/README.md`: the tutorial examples |
| Commands | **Command index (live)**: every command of the running interpreter with abbreviations, tier, availability and the first line of its description, with a filter |

* The left column lists the documents and the **Contents** (headings) of the open document, with a
  filter. A click on a heading, on a link to a section (`#...`) or on a link to another document
  (`../theory/THEORY_MANUAL.md#3-flexible-volume-substructuring`) opens the page and scrolls to the
  heading. **Back**, **Forward** and **Home** move through the pages you opened.
* **Search all documents** (Enter or Search) lists the sections whose text contains the words, with
  a snippet; a click opens the section.
* Links to the specification (`docs/spec/...`) open as well; external links open in a new browser tab;
  a link to a folder or a source file is shown as text.
* The renderer (`sassi/ui/markdown.py`) supports the Markdown the documents use: headings,
  paragraphs, lists, tables, block quotes, code, emphasis, links and images. It never passes raw HTML
  through: everything else is shown as text.

---

## 9. Verification (Help > Verification)

Choose a tier (P0, P1, P2, all) or a problem id (`VP-30`) and press **Run**. This runs
`VERIFY,<selection>` in a worker process. The table lists each check with its computed value,
reference, error, tolerance and pass/fail, plus the notes of each problem (requirements 6.1).
**List problems** shows the registered problems with their tier and modules (double-click one to
select it).

---

## 10. Typical workflow

1. **Model > Input** (`INP,model.pre`), or type commands in Command Entry. Set the name and folder
   with `MDL,<name>,<path>`, or use Model > Open with a database entry.
2. **Plot > Model > Elements** to check the geometry. Use the 3D toolbar for wireframe, shrink,
   colours by group, material or property, labels, fixed DOFs and masses. Interaction nodes are red.
3. **Options > Analysis**: set SITE (layers, `TOPL`, waves, frequency set), POINT (radius), HOUSE,
   ANALYS and MOTION (output nodes `NOUT`, damping ratios, control motion file). On the AFWRITE tab,
   choose the modules.
4. **AFWRITE** (toolbar). If CHECK reports messages, the **Check Errors** window opens by itself
   (unless Options > Check > Suppress Error Window is set; View > Check Errors opens it any time).
   Modules with errors get no deck (D-AFW-01).
5. **Modules > SITE, POINT, HOUSE, ANALYS, MOTION** (and RELDISP, STRESS). Each run streams its
   listing into an output tab and the Command History.
6. **File > Results Browser**. Click **Plot** on a `.TFU`, `.RS` or `.ACC` file to open the Line
   Selection dialog for it, or **Open** to read a listing.
7. **Plot > Spectrum TFU-TFI**: add file columns as lines, check them, set log axes, then OK. Open the
   dialog again to add more lines to the same view. **Options > Windows Settings** on the plot opens
   Graph Plot Options, with Average, SRSS, Broaden, Linear Combin, Addition and Subtraction. Their
   results are new lines; save them with `WRITESPEC`.
8. **File > Export Image** / **Export Table** for reports. **Model > Save** before Exit.

To keep the session as a `.pre` file, open a file with **File > Open** and press
**Input > Connect to Command Entry**. Every accepted command (typed, or sent by a menu or dialog) is
then appended to that file. Only one editor can be connected at a time.

Tutorial example 1 is a good first session: copy `examples/` somewhere, start
`sassi-gui --model-dir <copy>`, type `INP,ex01_surface_stick.pre` in Command Entry (about 20 s: it
runs every module), then plot `ex01/00041TR_X.TFU` and `ex01/00041TR_X01.RS` with Plot > Spectrum
TFU-TFI.

---

## 11. Files, settings and safety

* **Readable folders.** The file dialogs, the Results browser and the File Editor work only inside
  the model folders of the models in memory, the working directory, the start folder
  (`--model-dir`) and the course workspace folder (Learn). A name that leads outside them (`../..`, absolute paths, symbolic links) is
  refused. To work in another folder, use **CD** in the file dialog, or type `CD,<path>` or
  `MDL,<name>,<path>`. Commands typed by you (`INP`, `READSPEC` ...) follow the normal path rules of
  the interpreter (L15).
* **Paths with commas.** Menus and dialogs put a path or text that holds a comma in double quotes
  (`MDL,m,"/Users/me/Smith, J"`), so the folder stays one argument. Do the same when you type such
  a command.
* **Network.** The server binds to `127.0.0.1` only, refuses requests whose `Host` header is not a
  loopback name, and requires a per-session token that only its own page knows. Another web site
  therefore cannot send commands to it. The page's Content-Security-Policy allows scripts from the
  server only, and the plot mode bar has no cloud upload button.
* **Help pages.** Only Markdown documents of the help set are rendered, and `/docs/<path>` serves only
  the PNG, JPEG and GIF figures under `docs/`.
* **No shell.** Module runs start a Python worker process, or the executable chosen in Modules >
  Location, with an argument list.
* **Settings files** (requirements 5.10):

| File | Holds | Written by |
|---|---|---|
| `SASSIini.xml` | module locations, extensions, Command Display filters, message colours, shader options, the course workspace folder (Learn) | OK of Location / Extension / Colors / Shader Options (any `SHADEROPTIONS`), the Command Display toggles, Exit and Ctrl-C |
| `SASSIdb.xml` | Load Model groups: model name, folder, title | the Load Model dialog. The second Remove question deletes only the files named after the model (`<name>.*`, `<name>_*`). |
| `SASSIani.xml` | processed animations | `PROCFRAME`; Remove Animation in Load Frame Data, which can also delete the frame store files |

---

## 12. JSON API (for scripts and tests)

All `/api/` requests carry the header `X-SASSI-Token: <token>`; the token is in the page's
`<meta name="sassi-token">` tag. See `sassi/ui/api.py` for the full list.

| Request | Purpose |
|---|---|
| `POST /api/command {"line": "..."}` or `{"lines": [...]}` | execute command text; returns per-line `ok` and the messages. A typed `RUN<MODULE>` returns the started `job` and the `deferred` lines. |
| `GET /api/state`, `GET /api/events?since=N&timeout=T` | session state; long poll of messages, plot events, progress, jobs, and `check` (a CHECK with messages) |
| `GET /api/model[?number=N]` | nodes (global coordinates), fixities, interaction flags, elements by group, masses |
| `GET /api/options/<ANALYSIS or tab or MODEL/WRITE/CHECK or LOADGEN/LOADGENDYN>` | dialog layout and current values |
| `POST /api/options/<name> {"records": {...}, "indexed": {...}, "lists": {...}, "requests": {...}, "strings": {...}, "xrecords": {...}, "xtables": {...}, "xindexed": {...}, "dry_run": false}` | validated dialog commit; returns the command text. `xrecords` / `xtables` / `xindexed` hold the command records outside the option records (`EQL`, `NLSOIL`, `LOADGEN`, `LGFILE` ...; the tables `P`, `S`, `NLSLAYER`; the curves `BBC`, a `null` curve is deleted), see `sassi/ui/cmdrecords.py`. LOADGEN / LOADGENDYN with `"run": true` then start `RUNLOADGEN` and return its `job` (or `run_error`). |
| `POST /api/run/<MODULE> {"args": [...], "model": n}`, `GET /api/jobs/<id>?since=k`, `POST /api/jobs/<id>/cancel` | module runs (`args`: `["STATIC"]` or `["DYNAMIC"]` for LOADGEN only) |
| `GET /api/files?dir=`, `GET /api/file?name=`, `GET /api/fileinfo?name=`, `POST /api/file {"name", "text"}` | files (path safety) |
| `GET /api/plots`, `GET /api/plot/<id>[?frame=k]`, `GET /api/lines` | plot state (active plot, its lines, titles and settings, the 2D defaults) and plot data (`sassi.plotting.state.plot_data`) |
| `GET /api/help`, `GET /api/help/doc?name=<path>`, `GET /api/help/search?q=` | the help documents, one document rendered to HTML (title, headings), the sections containing a text |
| `GET /docs/<path>` | a figure of the documentation (no token: images only) |
| `GET /api/check`, `POST /api/verify {"select"}`, `GET /api/about` | windows of the View and Help menus |
| `GET /api/lessons`, `GET /api/lessons/<id>` | Learn: the course outline by part (and lesson files that do not parse); one lesson rendered for the lesson panel (HTML of every piece, command blocks, sections, actions) |
| `POST /api/lessons/<id>/open {"confirm"}` | a fresh workspace `<course root>/<id>/`; returns the commands the browser submits (`CD`, fresh model, setup), or `needs_confirm` with the unsaved models |
| `POST /api/lessons/<id>/action {"verb", "args"}`, `POST /api/lessons/<id>/progress {"step"}` | a step button resolved to command text or a GUI request (file, dialog, help page, explainer); a step has run in the current workspace |
| `GET /api/examples`, `POST /api/examples/<name>/prepare {"run", "confirm"}` | the examples gallery; copy an example into a fresh workspace and return the commands (`CD`, fresh model[, `INP`]) and the copied `.pre` |
| `GET /api/explain?line=<command line>&cursor=<k>` | the command explainer: name, syntax, meaning, tier, status, what it configures, the arguments with the given values; `cursor_arg` for the live hint |
| `POST /api/settings/course {"course": {"root"}}` | Learn > Course Workspace Folder (saved in `SASSIini.xml`) |
| `GET /static/katex/<file>`, `GET /static/katex/fonts/<file>.woff2` | KaTeX (LaTeX typesetting of the lessons and Help pages), served from `sassi/ui/static/katex` |
| `GET /api/course/workspaces`, `POST /api/course/workspaces/delete {"names"}` | Learn > Free Disk Space: the course workspaces with their size; delete those not in use (only folders with the course marker) |

---

## 13. Learn: the guided course

The Learn features turn the GUI into a guided educational resource for **structural engineers who design
nuclear (and other safety-related) facilities but are new to SASSI and to soil-structure interaction**.
They start from what you know (fixed-base ANSYS models, design spectra, ISRS, member forces) and show
what SASSI does, step by step, with the commands, their meaning and their technical basis. Everything
they do is command text sent to the same interpreter (rule L17): the Command History of a lesson is a
replayable `.pre` file.

### 13.1 The Learn tab (start page)

The **Learn** tab opens when the GUI starts (the check box *Show this page when the GUI starts* turns
this off) and from **Learn > Start Page** or the toolbar button with the graduation cap. It has:

* a **welcome** for the learner, with **Continue** (the lesson and step you were at), **Start the
  course** / **Next lesson**, **Open an example**, and the course workspace folder (**Change...**,
  **Free disk space...**);
* **How SASSI works in 2 minutes**: the module chain EQUAKE / SOIL → SITE → POINT → HOUSE → ANALYS →
  MOTION / RELDISP / STRESS. A click on a module shows what it computes, its input and output files, the
  ANSYS analogy where there is one, and buttons for the lesson that teaches it, the Theory Manual section
  and the explanation of its `RUN<MODULE>` command;
* **The course**: the lessons grouped by part, each with its number, title, summary, time, objectives,
  example, prerequisites and a progress bar (steps completed), **▶ Explainer** (the lesson's video,
  section 13.7) and **Start**, **Continue** or **Review**;
* **Examples**: one card per `examples/*.pre` (section 13.5);
* **While you work**: how to explain commands and try your own variations.

Progress (the steps you completed and the step you were at) is kept **per browser** (`localStorage`); a
private window or cleared site data simply starts without progress. **Reset progress** on a lesson card,
or **Learn > Reset Course Progress**, forgets it; the workspaces on disk are not touched.

### 13.2 The lesson panel

A lesson opens **full width, in a tab of its own** (*Lesson 3 · step 4/8*). While a step runs, the lesson
stays in front: the module output tabs and the plots its commands make open behind it, and every line
is marked as it runs. A step button (plot, listing, file) brings its result to the front, and the
**◀ Back to lesson** button at the bottom right of the tab area (or the Lesson tab) returns. **×** in the
lesson header (or closing the Lesson tab) closes the lesson; **Learn > Continue** reopens it.

**Side panel** in the lesson header (or **Learn > Lesson Layout**) switches to the other layout: a panel
docked to the right of the tab area, so that the plots, the listings and the Command History stay
visible next to it (at 1280 px the panel is 430 px wide, at 900 px 370 px; drag its left edge to resize
it, double-click the edge for the default width, **»** collapses it to a strip). **Full width** switches
back. The layout is remembered in this browser.

**Equations** in the lessons (and in the Help pages) are written in LaTeX and typeset with KaTeX, which
the GUI serves itself (`sassi/ui/static/katex`, no internet connection needed). A formula KaTeX cannot
typeset shows its LaTeX source with a red outline (the error is its tooltip);
`tests/unit/test_lessons.py` checks every formula of the course with the same KaTeX in Node.js.

**Concept figures** sit in the lesson text where they teach: a framed drawing with a title, sliders and
buttons for the parameter it is about, a line of computed values, its formulas and a caption. They are
computed in the browser from the formulas they show (`sassi/ui/static/figures.js`, no server request and
no command in the Command History), with defaults from the lesson's model; what is illustrative is
labelled in the figure's note. Animated figures have **Play / Pause**; they pause while scrolled out of
view and start paused when the system asks for reduced motion. On some charts a click or a drag sets the
frequency (or, in lesson 10, adds a computed frequency). The settings of a figure are kept while the
lesson stays open.

| Lesson, step | Figure | Computed from |
|---|---|---|
| 1, step 1; 5, introduction; 11, introduction | free field + structure − excavated soil = SSI system, built up term by term | manual Eq. 2.1 (schematic) |
| 1, step 3 | control motion → H(f) → floor motion → oscillators → ISRS, building up as the record plays | IFFT[H A]; Nigam-Jennings oscillators |
| 1, step 4; 4, step 9 | a one-mode stick on a fixed base and on sway-rocking springs and dashpots, with wave fronts for the radiated energy and \|H(f)\| of both | 3-DOF steady state with the lesson 4 impedance of the mat |
| 1, step 5; 5, step 10 | SSI = kinematic interaction + inertial interaction | the exact superposition of the two problems (schematic drawing) |
| 2, step 6 | shear wave in the layered column: displacement profile, outcrop and within amplification | exact layer solution (SHAKE recursion, complex modulus) |
| 3, step 7 | spring and dashpot under harmonic motion: phase lag and force-displacement ellipse | K(ω) = k + iωc, E_D = πcωU², ξ = ωc/2k |
| 5, step 9; 11, step 5 | the FV, FI-FSIN, FI-EVBN and FFV interaction nodes on the excavation mesh | the INTGEN sets (150, 105, 114, 132 nodes) |
| 7, step 3 | vertically propagating SV, SH and P waves (X, Y, Z input) | incident + reflected pulse in a half-space |
| 7, step 8 | envelope of soil cases and ±b peak broadening | BROADEN window maximum |
| 9, steps 1 and 8 | backbone, Masing loop, secant stiffness and ξ = E_D/(4πE_S) against the amplitude | hyperbolic soil / BBCGEN panel backbone, Masing rule |
| 10, step 3 | computed and interpolated transfer function, CRITFREQ flags, adding the flagged frequency | MOTION interpolation option 1, CRITFREQ |

* **Header**: part, lesson number and time, the title, a progress bar and one pill per page: **i** the
  introduction, **1 ... n** the steps. A green pill with ✓ is a completed step; a pill underlined in
  green has run in the current workspace. **⟲** resets the lesson.
* **Introduction**: summary, *You will be able to* (the objectives), time, example, *Builds on* (the
  prerequisite lessons, as links), the introduction text, the workspace folder and **Start ▶**.
* **A step**: the narrative; the commands of the step as a code block, where every line has **?** (or a
  click on the line: the command explainer, section 13.4) and **✎** (copy the line into Command Entry, to
  edit and run your own variation); blocks marked *Syntax / variation* are shown but not run. Then the
  step's **buttons** (*After running this step* / *See the results*) and its sections as labelled panels:
  **What this does**, **Why it matters**, **Technical basis**, **In ANSYS terms**, **Try this** and
  **Check yourself** with **Show answer**. Links to `docs/...` open in the Help tab at the heading.
* **Run step** submits the step's commands one by one, exactly as if typed in Command Entry: each line
  is marked ▶ while it runs, ✓ when it succeeded and ✗ at an error (hover for the message; the run stops
  there). A `RUN<MODULE>` line runs in a worker process with its listing streamed into the
  `<MODULE> output` tab and the Command History (section 7); the next line waits for it. **Stop** stops
  after the current command and cancels a running module.
* **Run steps 1–k** runs, in order, the steps up to the current one that have not run in this workspace
  yet (for example after reopening a lesson at step 5).
* **◀ Previous / Next ▶**; a step without commands counts as completed once read.

**Opening a lesson** creates a fresh workspace `<course root>/<lesson id>/` (default course root:
`<GUI start folder>/sassi-course/`, **Learn > Course Workspace Folder...**), copies the lesson's example
`.pre`, `examples/data/` and the files the example reads (e.g. `sassi/data/dynp_library.pre`), and runs
`CD,<workspace>`, a fresh model (`ACTM,0`, `DMODEL` of the models in memory) and the lesson's setup
commands. When a model in memory has changes that were not saved and it is not a model of a course
workspace, the GUI asks first. A lesson already open in this server session (after a page reload) is
shown as it is, without a new workspace. The GUI only ever deletes a workspace folder that it created
(marker file `.sassi-course-workspace`). Opening a lesson (or **⟲** Reset) also closes what the previous
lesson opened: its plots (`ACTIVATEPLOT,<id>` + `CLOSEPLOT`, in the Command History), its finished module
output tabs and its File Editor tabs of course-workspace files (a file with unsaved edits asks first);
tabs that were open before the lesson stay.

Workspaces stay on disk after a lesson, so you can come back to the results (most lessons leave a few
MB; lessons 8, 9 and 11 leave several tens of MB). **Learn > Free Disk Space...** lists the workspaces with
their size and deletes those that are not in use (the working directory, a model in memory or the open
lesson keep theirs); a lesson recreates its workspace the next time it is opened.

### 13.3 Step buttons (actions)

| Button | What it does | Command text |
|---|---|---|
| Plot the model / nodes / soil layers | Plot > Model > Elements / Nodes, Plot > Soil Layers | `MODELPLOT` / `NODEPLOT` / `LAYERPLOT` |
| Plot soil curves | Plot > Soil Properties of a DYNP property | `SOILPROPPLOT,<label>` |
| Plot *files* (spectra, transfer functions) | reads each file (first data column) as a new line and plots them; log frequency axis when the lesson asks; files of the same name are labelled with their folder | `READSPEC,<file>,1,<n>` ..., `SPECPLOT,<n>,...` (+ `LINENAME`, `AXES,,,,,1`) |
| Plot history *files* | reads each history and plots them | `READTH,<file>,<pair>,<n>` ..., `THPLOT,...` |
| Open *file* / Open the *MODULE* listing | a File Editor tab (the listing is `<model folder>/<model>_<MODULE>.out`) | |
| Open Options > Analysis > *TAB* | the Options dialog at that tab | |
| Read: *document* | the Help tab at the section | |
| Explain *command* | the command explainer | |
| Animate: *title* | stores the frames a step wrote (a `HARMFRAME` folder or restart frames) in `<folder>_ani` and plays them on the active model, undeformed shape in grey; the scale is the lesson's (the same for animations to compare) or 15 % of the model size for the largest displacement | `PROCFRAME,<folder>,<folder>_ani,<title>,3`, `DEFORMPLOT,<folder>_ani,1,<N>,<stride>,<scale>`, `WINDOWSETTINGS,UNDEFORMED,1`, `CNGVIEW,-90,0,0` (a front view, when the lesson asks for one), `WINDOWSETTINGS,TITLE,<title>` (or `VECTORPLOT` / `BUBBLEPLOT` / `CONTOURPLOT`) |

A button that shows results of a step that has not run in this workspace yet asks whether to run the step
first. File names are relative to the lesson workspace and may not leave it.

### 13.4 The command explainer

**Click a command line** in a lesson or in the Command History (the echo lines `> ...`), press an
**Explain** button, or use **Learn > Explain a Command...**: a panel shows

* the command's full name (an abbreviation such as `GROU` is resolved and named as such), a summary in
  engineering terms, the syntax, what it configures (model nodes, SITE module options, a module run ...),
  the tier and whether this build implements it, and the Command Reference category;
* a table of its arguments: position, name, **the value given in the line** (a coded value is decoded,
  e.g. `GROUP,1,2` → `2 = BEAMS`), the meaning and the default; values beyond the command's arguments are
  flagged;
* the dialog that edits the command (**Open Options > Analysis > SITE tab** ...) and **Command
  Reference**.

The argument names of the module option commands are those of their records (the names WRITE and the
dialogs use); the meanings come from the curated table `sassi/ui/command_args.json`, written from the
argument definitions of `docs/spec/07` to `11` and checked against the command handlers (where the two
differ, the table describes what SASSI-EDU does and says so). Every command used by the tutorial
examples is covered.

**Live argument hint.** While you type in Command Entry, the bar above it shows the command's arguments
with the one under the cursor highlighted, and its meaning and default (`SITE,0,1,0,12` → `nl`: number of
generated half-space sublayers ...). **details** opens the full explanation. The **?** button next to
Command Entry switches the hint on and off (remembered by the browser); pressed again while text is
typed, it opens the full explanation. Nothing is executed by an explanation.

### 13.5 Examples gallery and Model > Open Example...

One card per tutorial example (`examples/*.pre`): a picture of the model, title and *What you learn* (from
the comment header of the `.pre`), topic, modules and run time (from `examples/README.md`), the physics or
model paragraph, and:

* **the picture** (the default isometric view: element colours by group, interaction nodes red, a few lumped
  masses as dark diamonds; example 4, which has no structure, shows its soil column): a click copies the
  example into a fresh workspace as **Load into workspace** does, then submits its model commands (the `.pre`
  without `CHECK`, `AFWRITE`, `RUN<MODULE>`, `WRITE`, `FCOPY` ..., up to the first change of model: no module
  runs) and `MODELPLOT` -- `LAYERPLOT` for example 4 -- so the 3D model view opens in a second and the Command
  History rebuilds the model. The pictures are `sassi/ui/static/examples/<name>.png` (480 x 300 px), written by
  `python -m sassi.ui.thumbnails` (run it again after changing an example);
* **Open guided lesson** (or **Lesson n** when several lessons use the example);
* **Load into workspace**: copies the example (with `data/` and its support files) into a fresh folder
  `<course root>/examples/<name>/`, submits `CD` there and a fresh model, and opens the `.pre` in the File
  Editor: read it, change it, and press **Run (INP)**;
* **Run all**: the same, then `INP,<name>.pre` (the run time is on the card).

The originals in `examples/` are never changed. **Model > Open Example...** shows the same cards in a
dialog.

### 13.6 Lesson files

The lessons are Markdown files in `sassi/ui/lessons/` (format: `docs/internal/lesson_format.md`). The GUI
re-reads them on every request, so a lesson being written appears in the course at once; a file that
cannot be parsed is listed under the course with its error instead of hiding the others.

### 13.7 Explainer videos

Every lesson has a narrated **motion-graphics explainer video** (about 3-4 minutes) for newcomers: the one
big idea of the lesson in plain words and everyday pictures, how SASSI computes it, and what it changes
for a design, with the lesson's own results (the curves are the lesson's runs, `python -m sassi.ui.video_data`). They open in
a browser tab of their own from **▶ Explainer** on a lesson card, **▶ Watch the explainer video** on a
lesson's introduction, or **Learn > Explainer Videos** (the list of all eleven); the files are
`sassi/ui/static/videos/NN.html` and also open directly from disk.

The videos are HTML, SVG and JavaScript; the **voiceover is recorded** with ElevenLabs (one small MP3 per
sentence, `web/voice_videos.mjs`), with the browser's text-to-speech as an alternative voice (⚙) and as the
fallback for a sentence without a recording. The player has captions (**c**), a
transcript with the chapters (**t**; a click on a sentence jumps there), play / pause (**space**),
previous / next sentence (**← →**) and chapter (**shift ← →**), voice on / off (**m**, captions only),
full screen (**f**), and ⚙ for the voice and the speed (1× to 1.6×; 1.15× by default). Authoring contract:
`docs/internal/explainer_videos.md`; `node web/test_videos.mjs` plays every beat of every video in
headless Chrome.

## 14. How the GUI is tested

* `tests/unit/test_ui_api.py`, `test_ui_dialogs.py`, `test_ui_jobs.py`, `test_ui_server.py`: the API,
  the dialog layouts and commits, module runs and Cancel, the HTTP server and its safety rules. The
  command-record bindings (SOIL Nonlinear Soil, NONLINEAR tab, the two LOADGEN dialogs) are pinned by
  round trips (open, change, OK = the command text, reopen, WRITE, posting the values back emits nothing),
  the refused values with the commands' messages, the NONLINEAR and LOADGEN jobs (listing, progress,
  replay history, typed `RUNLOADGEN`) and the 2D PLANE model in the 3D model view.
* `tests/unit/test_ui_flows.py`: the menu, toolbar and dialog flows as command / API sequences: every
  Options > Analysis tab round-trips (open, change, OK, reopen, WRITE), Model / Write / Check, Line
  Selection (added lines, pre-fill state, replace), converters with an output `.pre`, Load Model and
  Exit, Process Animation Frame List on a frame folder and the four animations, the menus and toolbars
  of requirements 5.3 / 5.9, regression guards of the front end and a JavaScript syntax check.
* `tests/unit/test_ui_markdown.py`: the Markdown renderer (blocks, inline, safety) and the help set
  (every link between documents reaches an existing heading).
* The front end was also exercised by hand in a browser against tutorial example 1, at 1024 px and at
  800 and 760 px window widths: every menu item, toolbar button and dialog of sections 3 and 4. The
  NONLINEAR and SOIL tabs, the two LOADGEN dialogs (Ok and Run), Modules > NONLINEAR and a 2D PLANE
  model plot were checked the same way at 1280 and 900 px.
* `tests/unit/test_ui_learn.py` and `test_ui_learn_explain.py`: the Learn API (course outline and lesson
  rendering, workspace preparation and setup, the unsaved-model question, Run step with a module run
  streamed through a worker job, every step action and its path safety, the examples gallery, loading and
  running an example, the course folder setting) and the command explainer (abbreviations, record
  arguments, coded values, unknown commands, the cursor argument, the curated table against the command
  catalogue and the record fields, coverage of every command of the examples), plus front-end guards. The
  Learn tab, the lesson panel (lessons 1 and 2: Run step, Run steps 1–k with five module runs, every kind
  of step button, Show answer, Reset, collapse, page reload), the explainer from lesson lines, history
  lines and Command Entry, Model > Open Example... (Load into workspace, Run (INP)) and the Learn menu were
  exercised in a browser at 1280 and 900 px.
* `tests/unit/test_web.py`: the web version (section 16): inline module runs and the bridge (a typed
  `RUN<MODULE>` queued by its request and run afterwards, its listing pushed and shown once, the deferred
  commands and the replay history, Cancel before the start), the event poll that never waits, the files of
  the built site (nothing of `reference/`, `docs/spec/`, `docs/internal/`, `tests/` or the course
  workspaces; every file the page names exists) and, when `npm install` was run in `web/`, the built bundle
  in Pyodide under Node (`web/test_pyodide.mjs`: lesson 1 driven like the page; without `--quick` every
  lesson headless).

---

## 15. Notes and limitations

* Plots are drawn by Plotly in the browser. `CAPTUREPLOT` (File > Export Image) still uses the
  headless matplotlib renderer, so exported images look like the batch images (UI-03).
* Mouse controls follow Plotly: drag to orbit, scroll to zoom, right-drag to pan. The camera of the
  plot (`CNGVIEW`) sets the initial view; mouse rotation is a display action and is not recorded.
* Tier P1/P2 features (incoherency, multiple excitation, Option NON, SOIL-NON, Option A, binary
  databases, Font) appear as in the manual. When they are not available in this build they are greyed
  or marked, and the interpreter prints the tier message when such a command is stored. Option NON (the
  NONLINEAR tab and Modules > NONLINEAR), SOIL-NON (SOIL tab) and Option A (the two LOADGEN dialogs) are
  available and editable; nonlinear beams, the Material Parameter and ANSYS Super Element Utilities stay
  greyed.
* A module set to an external executable in Modules > Location must find its input deck; for NONLINEAR
  and LOADGEN the decks (`.eql`, `.lgn`) are written by the built-in `RUNNONLINEAR` / `RUNLOADGEN`, so an
  external program of these two modules needs a deck from an earlier built-in run.
* Modules > Extension values are stored and passed to external executables. The built-in modules use
  the deck extensions of the manual.

---

## 16. Web version (GitHub Pages)

The GUI also exists as a static web site that anyone can open in a browser, with nothing to install.
Python runs in the page itself: [Pyodide](https://pyodide.org) (CPython compiled to WebAssembly, with
NumPy and SciPy) in a Web Worker runs the same `sassi` package, the same API and the same front end as
`sassi-gui`. No server computes anything and nothing you type leaves your computer.

* **First start.** The first visit downloads about 27 MB (Python, NumPy, SciPy, Plotly.js and SASSI-EDU)
  and keeps it on your computer, so the next visits download nothing -- the site even works offline -- and
  only start Python (about 3 s). Python, NumPy and SciPy are kept by the Python engine itself, also in a
  browser where the site's service worker cannot run. A new version of the site downloads only the files
  that changed (usually the SASSI-EDU bundle, about 2 MB). A progress panel shows what is loading, with the
  disclaimer; the Command History then says how Python started ("came from this computer: nothing was
  downloaded", or what was downloaded). A current desktop browser (Chrome, Edge, Firefox, Safari) is
  needed; a private window may refuse to keep files, and then downloads them at every visit.
* **Same course and tools.** The Learn tab, the lessons with their module runs, the examples, plots,
  Help (with the formulas and figures), the command explainer, the Options dialogs and Command Entry
  work as described above. Module runs show their listing and the status-bar progress while they run.
* **Files are kept in the browser.** The workspace is `/home/pyodide/work`; the course workspaces are under
  it. It and the settings folder are saved in the browser's storage (IndexedDB) after every change and
  when you leave the page, and restored at the next visit, so decks, results and saved models survive a
  reload; the models in memory do not (Model > Save or Output keeps one). **File > Clear Saved Files...**
  deletes what is saved and starts afresh (your course progress is kept). In a private window the files
  live in the tab only.
* **Your computer's files.** Every file dialog has **From your computer…**: it copies text files from your
  computer into the folder the dialog shows and selects the first, so Model > Input, File > Open and the
  file fields of the Options dialogs read your own decks directly. **File > Upload to Workspace...** does
  the same for several files, **File > Download...** saves a text file of the workspace (a listing, a deck,
  a `.pre`, a spectrum) on your computer, and Model > Output, Export Table, Export to ANSYS and Export Image
  download what they write. Binary files (FILE1 ... FILE8, `.sdb`) cannot be transferred.
* **Differences.** A running module cannot be cancelled (it runs in the page's Python engine; other
  requests wait for it, Stop in the lesson panel acts after the current command). Model > Exit and
  Modules > Location (external executables) are not shown. File > Export Image saves the active plot as
  a PNG with Plotly instead of `CAPTUREPLOT`. Calculations take 1.2 to 1.6 times as long as with the
  local `sassi-gui`; Help > Verification runs, but slowly.

Building, testing and publishing the site: `web/README.md` in the project (`python web/build.py`, then
push to GitHub with *Settings > Pages > Source: GitHub Actions*).
