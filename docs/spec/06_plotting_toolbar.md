# 06 — Plotting Operations (3D model plots, 2D line plots, soil plots, animations) and Toolbar Reference

Source: ACS SASSI V3 User Manual, Chapter 7 "Plotting Operations and Examples" (7.1–7.3) and Chapter 8 "Toolbar
Reference" (8.1–8.2).
Text lines 8044–8690 of `reference/acs-sassi.txt` (printed pages 191–215; Read-tool PDF pages 193–217). I viewed
every page in this range as an image to transcribe the dialog boxes, plot screenshots and toolbar tables.

Chapter 7 is short and mostly descriptive. Exact syntax, defaults and data semantics for most features live in
other parts of the manual, so I read those as well. They are cited inline and summarised only as far as this
section needs them:

| Topic | Where in the manual | Used for |
|---|---|---|
| Plot commands 9.14.1–9.14.52 (ADDITION … YTITLE2) | lines 12230–13020, printed pp. 291–307 | command syntax for every plot, toolbar mapping, line-object maths |
| Plot menu tree (Section 6 menu table) and 6.3.1–6.3.12 | lines 3930–3960, 4405–4451 | menu names |
| 6.2.2 Export Image, 6.2.3 Export Table | lines 4387–4396 | image and CSV export |
| 6.5.5–6.5.9 Windows Settings, Colors, Font, Shader Options, Reset Plot | lines 7849–7944 (also in spec `05d`) | palettes, fonts, shader defaults |
| Section 3 module descriptions, Tables 3.1 and 3.2 (images on PDF pp. 49–50) | lines 2040–2290 | result text files and frame-file naming schemes that feed the plots |
| MOTION, STRESS and RELDISP post-processing "Restart" options | lines 6946–6977, 7318–7395, 7451–7482 | which frames feed which animation |
| FRAMECOMBIN, FRAMESEL, MODFRAMES (9.7.11, 9.7.12, 9.7.24) | lines 10890–10924, 11070–11081 | frame-file header and list-file conventions |
| Binary database commands 9.18 (ACCDBANI, DISPDBANI, THSDBANI, MAXDBFRAME, BINFRAMEOUT, LOAD*DB) | lines 14015–14325 | the "(Binary)" animation types |
| L, TOPL commands (9.4.18, 9.2.44) | lines 10278–10290, 9613–9624 | soil layer plot data |
| Graphing macro example 5.6.2 | lines 3585–3620 | READSPEC → SRSS → WRITESPEC → SPECPLOT → CAPTUREPLOT workflow |
| Model-size limits (Section 1.2) | lines 534–585 | soil-curve limits (100 curves, 11 points) |
| Commands SAVE, WRITE, AFWRITE, ANSYS, CONVERT | 9.2.29, 9.2.47, 9.2.3, 9.11.1, 9.11.4 | main-toolbar mapping |

Tags: **[CORE]** needed for computational correctness · **[IO]** file or input format · **[UI]** user interface,
plotting or convenience · **[ADV]** advanced (binary databases, incoherency, nonlinear, Option A/AA/PRO).
"(inferred)" means the manual does not state the item outright. I derived it from a screenshot or from
consistency with other sections, or I propose it as standard practice. Every "(inferred)" item also appears in
Section 13, Open questions.

> Plotting is mostly [UI], but several items are [CORE] for a *verifiable* program. These are the line-object
> maths (AVERAGE, SRSS, LINECOMBIN, ADDITION, SUBTRACTION, BROADEN, with interpolation and extrapolation), the
> WRITESPEC/WRITETH resampling rules, the reading of SSI result columns, and the mapping from frame data to
> deformed geometry, colour and vectors. Results the engineer reports, such as broadened ISRS and SRSS spectra,
> come out of these plotting tools.

---

## 0. Architecture overview

### 0.1 Plot families [UI]

The UI has three plot families. All of them open as **tabs** in the main window, next to the "Command History"
tab. Screenshots show the tab captions below (window title "SSI Submodeler"):

| Family | Plot (menu name) | Command | Tab caption seen in screenshots |
|---|---|---|---|
| 3D static | Element Plot | `MODELPLOT` | `Model 0 - Model Plot` |
| 3D static | Node Plot | `NODEPLOT` | `Model 0 - Node Plot` |
| 3D static | Cut Plot | `CUTPLOT,[Cut],[Model]` | `Model 0 - Cut # 1 Plot` |
| 2D line | Spectrum TFU-TFI (Spectrum Plot) | `SPECPLOT,<Line1>,…,<Line50>` | `Spectrum Plot` |
| 2D line | Time History Plot | `THPLOT,<Line1>,…,<Line50>` | `Time History Plot` |
| 2D special | Soil Layer Plot | `LAYERPLOT` | `Soil Layer Plot` |
| 2D special | Soil Properties Plot | `SOILPROPPLOT,<PropName>` | `Soil Property Plot` |
| 3D animated | Process Animation Frame List (pre-processing step) | `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]` | (no tab; progress bar) |
| 3D animated | Bubble | `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | `Model 0 - Bubble Plot` |
| 3D animated | Vector | `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | `Model 0 - Vector Plot` |
| 3D animated | Contour | `CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | `Model 0 - Contour Plot` |
| 3D animated | Deformed Shape | `DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | `Model 0 - Displacement Plot` |

The Plot menu also lists "Non Uniform Soil Field", marked "(not active)". Do not implement it, or show it disabled.

### 0.2 Equivalence rule: menu = toolbar = command [UI][IO]

* Every plot can be opened in three ways: from the **Plot menu**, from the **main toolbar**, or by typing its
  **command**. The manual says repeatedly that choosing the menu item "is the same as" the command.
* A plot is the **active plot** when its tab is in front. Most plotting commands act on the active plot, for
  example `PLOTTITLE`, `PLOTRANGE`, `CNGVIEW`, `WIREFRAME`, `PAUSE`, `CAPTUREPLOT` and `CLOSEPLOT`.
* Commands can come from a `.pre` file (batch). The manual stresses that plotting commands allow "batch image
  manipulation using a .pre file". Example from 5.6.2 (lines 3594–3603):
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
  **Implementation requirement [CORE]:** plots must render correctly with **no interactive window**, so that
  `CAPTUREPLOT` works in batch mode. Use an off-screen or headless renderer.
* If a plot command that needs arguments gets none, it opens the matching popup dialog (Cut, Line Selection,
  Select Dynamic Soil Property, Load Frame Data, Parse Frame Data). The manual warns that for the animation
  commands, "Partial argument entry may cause inconsistent results". **Recommended:** when an animation command
  gets some but not all of its arguments, either reject it with an error or fill the missing ones from the
  defaults stored in `SASSIani.xml`. (inferred)

### 0.3 Session state that plots depend on [CORE][UI]

| State | Scope | Notes |
|---|---|---|
| **Models in memory** | by reference number (`Model 0`, …) | One model is *active*. Every 3D plot except Cut Plot draws the active model. |
| **Hide requests** (hidden groups, elements, nodes) | **stored in the model data**, shared by all plots of that model | This is new compared with PREP: you do not re-enter hide requests for each plot. Whether they are saved to disk by `SAVE` is not stated (see Section 13). |
| **Selected nodes** (`NODESEL`) | model-level; applies to "all current and future plots" | Up to 20 nodes per command. They accumulate across commands, and repeating a node number deselects it. |
| **Line objects** | session-wide, by integer reference number | Loaded with READSPEC/READTH or the Line Selection dialog, or created by line maths. **Cleared when the UI closes.** Line *name* and *markers* are global line properties that show in every graph. |
| **Dynamic soil properties** | model data | Named G/Gmax–γ and D–γ curve pairs. They can be created and edited from the Soil Properties plot popup. |
| **Animation database** `SASSIani.xml` | per user, system-dependent default location | Lists processed animations: description, directory, type, frame count, defaults. |
| **Colour palettes, fonts** | UI settings | See spec `05d` §5 (Colors, Font). The `COLOR,<Palette>,<Num>,<R>,<G>,<B>` command edits them. |
| **Shader options** | applies to the current active plot (window) | Point size, outline thickness, shrink, scale factor. |

---

## 1. 3D plot controls (7.1.1) — common to all 3D plots and animations

### 1.1 Camera model [UI]

* Six numbers define the view. `CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>` sets them; `DEBUG` mode displays them
  on the plot:
  * `rX`, `rY`, `rZ`: rotations about the X, Y and Z axes (units not stated; degrees inferred);
  * `px`, `py`: horizontal and vertical screen pan;
  * `zoom`: zoom constant.
* The **centre of rotation** is a separate 3D point, set by `CNGCENTER,<X>,<Y>,<Z>`. `RSTCENTER` resets it to the
  default: the **centre of the axis-aligned bounding box of all nodes used by elements**, in the global system.
* **Default view** (`RSTVIEW`, Options ▸ Reset Plot, or the toolbar "Reset current plot view"): the "top down view
  of the model when the plot is first opened". The default zoom and position come from the same bounding box of
  element-connected nodes. The screenshots show an oblique view on opening, so "top down" may mean "the initial
  view". (see Section 13)
* A small **axis triad** (X red, Y green, Z blue, labelled) is drawn in the lower-left corner of every 3D plot.
  (from screenshots)
* A parallel (orthographic) projection is recommended, because engineering plots should not distort. (inferred)

### 1.2 Mouse controls [UI]

| Gesture | Action |
|---|---|
| Right click + drag | Rotate the model about the centre of rotation |
| Left click + drag (vertical) | Zoom in and out |
| Middle click + drag | Pan the model across the plot area |
| Shift + middle click + drag | Draw a zoom box; on release of the middle button, zoom to that box |
| Left double click | Highlight the **node** under the cursor and print its node number in the Command History |
| Right double click | Highlight the **element** under the cursor and print its element number in the Command History |

Picking (inferred): for a node, take the nearest *visible* projected node within a pixel tolerance, about 5 px.
For an element, cast a ray through the pixel and take the nearest visible face it hits.

### 1.3 Keyboard controls [UI]

| Key | Action |
|---|---|
| Insert / Delete | Rotate about the X axis (+ / −) |
| Home / End | Rotate about the Y axis |
| PageUp / PageDown | Rotate about the Z axis |
| Pause | Start or stop the animation (same as `PAUSE`) |
| `+` / `-` | Step forward or back one frame. The manual says to use this **only while the animation is paused**. |

The rotation step per key press is not stated. Use 5° per press. (inferred)

### 1.4 3D "Window Options" dialog (Windows Settings for a 3D plot) [UI]

This dialog opens from Options ▸ Windows Settings or the `WINDOWSETTINGS` command whenever a 3D plot is active.
Its title in the screenshot is **"Window Options"** (PDF p. 194). The same dialog serves every 3D plot. Controls
that do not apply to the active plot type are **greyed out**: the screenshot of an Element plot shows Colormap
Value Range, Output Direction, Scale Factor and Show Undeformed Shape disabled.

| Group | Field | Type | Default / example | Enabled for | Behaviour |
|---|---|---|---|---|---|
| Model Display Volume | X Min, X Max, Y Min, Y Max, Z Min, Z Max | 6 real edit boxes | Bounding box of the model, e.g. X 0.00–50.00, Y 0.00–80.00, Z 0.00–50.00 | all 3D | Geometry outside the box is hidden (a clipping box). By default it is always the box around the model. |
| Show Element Group | Check-list of groups | check list box | All checked. Labels look like `Group 1 SHELL - 70`, `Group 4 SOLID - 70` | all 3D | Uncheck to hide a group, check to show it. To toggle many at once, highlight several rows and press **Space**. The meaning of the number after the dash is not stated; it is probably the element count (inferred). |
| Colormap Value Range | Min, Max | 2 real edit boxes | Data range of the animation, e.g. 0.00 / 0.00 when not applicable | **Bubble and Contour only** | Sets the values at the ends of the colour map. Values below Min take the min colour (**dark blue**); values above Max take the max colour (**dark red**). |
| Output Direction | radio buttons: **X**, **Y**, **Z**, **All** | radio | X | **Vector only** | Draws the vectors in the chosen direction only, or all three at once. |
| Hide/Show Elements (or Hide/Show Nodes for node-based plots) | radio **Hide** / **Show** | radio | Hide | all 3D | Chooses the action of the button below. |
|  | Group | edit box | empty | element plots | **Only one group** per press. |
|  | Elem. Numbers (or Node Numbers) | edit box | empty | all 3D | Either a **space-separated list**, e.g. `3 7 12`, or a **single dash range**, e.g. `100-250`. Mixing a list and a range in one entry is not allowed: one list *or* one range per press. |
|  | **Hide/Show Elem.** button | button | — | all 3D | Applies the request **immediately**. It stays in effect even if the user then closes the dialog with Cancel. Press it as many times as needed. |
| Animation Options | Scale Factor | real | 1.00 | **Vector and Deformed Shape** | Multiplies vector length and displacement. |
|  | Frame Pause (ms) | integer | **33** | all plots, including non-animated ones | The refresh period. 33 ms is about 30 frames per second. |
| (bottom) | Show Undeformed Shape | check box | unchecked | **Deformed Shape only** | Overlays a wireframe of the original, undeformed model. |
| (bottom) | Title | edit box | empty | all 3D | Plot title (same as `PLOTTITLE`). |
| buttons | **OK**, **Cancel** | | | | OK applies the settings to the active plot. Cancel discards everything except hide/show requests already applied with the button. |

**Hide/Show list parsing [UI] (inferred grammar, consistent with the text):**
```
entry := int ( WS int )*        # list of ids
       | int '-' int            # inclusive range a..b, a <= b
```
Ignore ids that do not exist, or warn about them. Hiding an element does not delete it; the request is stored in
the model's hide set. The manual does not say whether an element is outside the Display Volume when *any* node
is outside or only when *all* nodes are. I recommend hiding an element only if **all** its nodes lie outside the
box. (inferred, see Section 13)

### 1.5 Shader Options dialog (Options ▸ Shader Options, `SHADEROPTIONS`) [UI]

The dialog title is **"Shader Options"** (PDF p. 195). It sets variables of the OpenGL shaders. These options are
in a separate window so that the Window Options dialog stays similar to the PREP one. **OK** updates the
*current active plot* with the values.

| Field (exact label) | Default | Command argument | Meaning |
|---|---|---|---|
| Node/Bubble Node Size | **10.00** | `[points]` | Maximum point or bubble size in the Node, Bubble and Vector plots (node dots). Units are presumably pixels. |
| Vector/Displacement Scale Factor | **1.00** | `[scale]` | Scale factor for the Vector and Deformed Shape plots. This duplicates "Scale Factor" in Window Options and Load Frame Data, so keep a single variable. (inferred) |
| Element Outline Thickness (% of element) | **0.02** | `[linew]` | Share of each element edge drawn in the outline colour instead of the element colour. Applies to the Element, Cut, Contour and Deformed Shape plots. The text calls it "thinness". |
| Element Shrink (% of element) | **0.06** | `[shrink]` | Amount each face shrinks in the Element Plot when shrink mode is on. |
| buttons | Cancel, OK | | |

Command: `SHADEROPTIONS,[points],[linew],[shrink],[scale]`. All arguments are optional and only the ones given
change. With no arguments the dialog opens.

**Interpretation of the "%" values (inferred):** 0.02 and 0.06 are **fractions**, i.e. 2 % and 6 %, not
0.02 %. Use these formulas:
* Shrink: for each face, `v' = c + (1 − s)(v − c)`, where `c` is the face centroid and `s` is the shrink value.
  Use the element centroid instead if you shrink whole elements.
* Outline: with barycentric or edge-distance shading in the fragment shader, draw a fragment in the outline colour
  when its normalised distance to the nearest edge is less than `t`, the thickness value.

### 1.6 Plot toolbar
The 3D plot toolbar is specified in Section 11.2. Its commands are listed in Section 10.

---

## 2. Element Plot (7.1.2) — `MODELPLOT`

### 2.1 What it shows [UI]
* The **active model**, with all elements drawn as **shaded, filled** faces. Colour depends on the scheme
  selected (Section 2.2). Element outlines are drawn using the outline thickness.
* Purpose: to see where elements sit relative to each other, to check connectivity between elements, and to find
  and fix modelling problems visually.
* Beam and spring elements are drawn as lines. The Element-plot screenshot shows a green stick for a beam on top
  of a red disc of shells or solids. Draw lumped-mass nodes as markers (see `SHOWMASS`). (from screenshot)

### 2.2 Colouring (toolbar G / M / P; `ELECOLOR,<val>`) [UI]
| `val` | Colour elements by |
|---|---|
| 1 | Group number |
| 2 | Material number |
| 3 | Beam/spring property number |

* Elements with no material or property number get a **default colour**.
* Colours come from the **ElemPalette**, which holds up to **128 colours**. It is filled as items are plotted.
  With more than 128 items, `colour index = item number mod 128` (manual: `Elemnum ≡ Colnum (mod 128)`). One
  palette serves groups, materials and properties.
* The rest of the Element plot colours (background, outline, labels and so on) come from the **Element** palette.

### 2.3 Display toggles available on the Element plot [UI]
| Toggle | Command | Argument semantics |
|---|---|---|
| Wireframe mode | `WIREFRAME,<Switch>` | 0 off, 1 on, −1 (default) toggle. Element Plot only. |
| Shrink mode | `SHRINK,[switch]` | 1 on, 0 off, −1 (default) toggle. Element Plot only. |
| Node numbers | `NODENUM,[opt]` | −1 toggle (default), 0 off, 1 on |
| Element numbers | `ELENUM,[opt]` | −1 toggle (default), 0 off, 1 on. **Cannot be on together with GROUPNUM**: turning one on turns the other off. |
| Group numbers | `GROUPNUM,[opt]` | −1 toggle (default), 0 off, 1 on. Cannot be on together with ELENUM. |
| Fixed DOF markers | `SHOWDOF,[label1],…,[label6]` | Labels X, Y, Z, XX, YY, ZZ, DISP (=X,Y,Z), ROT (=XX,YY,ZZ), ALL. With no labels, a check-box window opens (manual Figure 9.7). A marker is drawn on every node that has a fixed DOF in a requested direction. |
| Lumped-mass markers | `SHOWMASS,[opt]` | −1 toggle (default), 0 off, 1 on. **In the element plot the marker shows direction:** a red segment if the node has mass in X, green for Y, blue for Z. Other plots only show that a mass is present. |
| Debug overlay | `DEBUG,[switch]` | 0 off, 1 on, 2 (default) toggle. Shows the view angles (the CNGVIEW arguments) and animation information. |

### 2.4 Data required [CORE]
Nodes (id, x, y, z), elements (id, group, type, connectivity), group → element type, material and property ids
per element, boundary conditions (fixed DOFs per node), lumped masses per node and direction, the hide set, and
the selected nodes.

---

## 3. Node Plot (7.1.3) — `NODEPLOT`

### 3.1 What it shows [UI]
* The active model drawn as **points only**. Only nodes that are **connected to elements** are shown.
* Default colours, all changeable in the **Node** palette or with `COLOR`:

| Node state | Default colour |
|---|---|
| Ordinary (non-interaction) node | **black** dot |
| **Interaction node** | **red** dot |
| Node with a fixed DOF of the kind requested in SHOWDOF | **green border** around the node |
| Node with a lumped mass (SHOWMASS on) | **purple border** |
| Node selected by mouse double-click or `NODESEL` | **blue square border** |

* Node labels are toggled with `NODENUM`. (inferred to apply)

### 3.2 Data required [CORE]
Nodes, element connectivity (to find connected nodes), the **interaction-node set** (the nodes where the
structure meets the excavated soil; it comes from the model's interaction definition), fixed DOFs, masses and
selected nodes.

---

## 4. Cut Plot (7.1.4) — `CUTPLOT,[Cut],[Model]`

### 4.1 What it shows [UI]
* A **wireframe of the whole model**, with the elements that belong to the cut drawn **filled**. The screenshot
  shows red solids inside a black wireframe.
* Purpose: to check that a cut is complete, and to change its contents, before using it to build a **Submodel**
  or a **cross-sectional model** (section-cut).
* **This is the only plot type that can show a model other than the active model.**

### 4.2 Opening dialog "Select Cut to Display" [UI]
| Field | Default (screenshot) | Meaning |
|---|---|---|
| Cut Number | 1 | Number of the cut (a cut is an element set defined elsewhere in the model) |
| Model Number | 1 | Reference number of the model to plot |
| Ok / Cancel | | |

The dialog appears when the menu item is used, or when `CUTPLOT` is typed with no arguments.

### 4.3 Data required [CORE]
Model geometry and the cut's element list.

---

## 5. Two-dimensional line plots (7.2.1–7.2.4)

### 5.1 Line objects [CORE][IO]
* All 2D numeric data is held in **line objects**. A line object is an ordered list of (x, y) pairs plus a
  **name** (a global property, changed with `LINENAME,<Num>,<Name>`) and a **markers** flag (global, set with
  `MARKERS,<Mark>,<Ln1>…<Ln50>`, where 0 is off and 1 is on). Each line object has an integer reference number.
* Line objects last for the session only and are cleared when the UI closes.
* A line's colour depends on its **reference number**. The colours come from the SpecLines palette in spectrum
  plots and from the THLines palette in time-history plots. Both palettes are "populated as lines with different
  line numbers are plotted".
* **Loading into a reference number that is already in use overwrites that line** without warning.
* Ways to load:
  * `READSPEC,<SpecFile>,<numLines>,<Line1> … <LineN>`: reads a spectrum-format file. **The frequency column is
    not counted in numLines.** An `.RS` file with one frequency column and one acceleration column has
    `numLines = 1`. The same applies to `.TFU`, `.TFI` and `.TFD` files. The manual's note reads "TFI, TFI, TFD";
    the first one is presumably TFU. Each further data column becomes one line object, assigned in order to the
    reference numbers Line1…LineN.
  * `READTH,<THFile>,<Pair>,<Num>`: reads a time-history file. `Pair = 0` reads a one-column time history (the
    ACS SASSI output format). `Pair = 1` reads two columns of time/acceleration pairs. The result goes into
    reference number `Num`.
  * The **Line Selection** dialog (Section 5.2).

**Column-selection limitation (by design):** the dialog cannot pick individual columns. To load column k, you must
load **all columns 1…k**. The extra lines get reference numbers counting up from the starting number. (The
command form lets you list reference numbers, but it still reads the columns in order.)

### 5.2 Opening dialog "Line Selection" (7.2.1) [UI]
This dialog opens when Spectrum Plot or Time History Plot is chosen from the menu or toolbar. (The commands take
line numbers directly.)

| Group | Field | Type | Behaviour |
|---|---|---|---|
| (top) | Line list | check list | One row per line object in memory, shown as `<ref>: <name>`, e.g. `1: Test2.rs`. Check the lines to plot. |
| (top right) | Ok, Cancel | buttons | Ok creates the plot with the checked lines. |
| Title _Axis Labels | Title, X-Label, Y-Label | edit boxes | Plot title and axis titles (`PLOTTITLE`, `XTITLE`, `YTITLE`) |
| X Axis Options | Min, Max | real | Axis extent. Blank means use the data extent (min and max over the plotted lines). |
|  | Logarithmic | check | Log X axis |
|  | Show Ticks | check | Show minor ticks (thin grid lines) |
| Y Axis Options | Min, Max, Logarithmic, Show Ticks | same | same, for Y |
| Input Line File | File Name | edit box + `<<` browse button | File to load |
|  | Starting Number | integer | Reference number for the first column loaded |
|  | Lines in file | integer | How many data columns to load, counted from the first and excluding the frequency or time column |
|  | **Add Line(s)** | button | Reads the file and stores the lines in memory. The new lines then appear in the list at the top. |

All of these settings can be changed later with the plotting commands or with the Windows Settings dialog.

### 5.3 "Graph Plot Options" dialog (7.2.2 Spectrum & Time History Options) [UI][CORE]
This is the Windows Settings dialog for a Spectrum or Time History plot. The text calls it the "Line Settings
Window"; its title bar reads **"Graph Plot Options"** (PDF p. 202). It picks which lines in memory appear on the
**active** plot, sets scale and extent, and runs the **line calculations**.

| Group | Field / button | Default (screenshot) | Behaviour |
|---|---|---|---|
| (top) | Line check list | `1: Test.rs` checked | Lines shown on the active plot. **It also selects the inputs for the calculation buttons.** |
| (top right) | Ok, Cancel | | |
| Title _Axis Labels | Title, X-Label, Y-Label | empty | |
| X Axis Options | Min / Max | **0.10 / 100.00** (the data extent of an RS from 0.1 to 100 Hz) | |
|  | Logarithmic, Show Ticks | unchecked | |
| Y Axis Options | Min / Max | **0.02 / 4.12** (data extent) | |
|  | Logarithmic, Show Ticks | unchecked | |
| Line Options | Data Points (highlighted line) | unchecked, greyed out | Shows data-point markers on the highlighted (selected) line. Enabled only when a line is highlighted. (inferred) Same as the global MARKERS property. |
|  | Line Stippling (all lines) | unchecked | Dashed patterns so lines can be told apart in black-and-white printing (same as `STIPPLE`) |
| Spectra Analysis Postprocessing | Spectral Average: **Average** | | Same as `AVERAGE` over the checked lines |
|  | SRSS Combination: **SRSS** | | Same as `SRSS` over the checked lines |
| ⤷ Broadening and Enveloping Spectra | Peak Difference(%) | **0** | `Smooth1`: peak-bridging percentage |
|  | Broaden(%) | **0** | `Smooth2`: peak-broadening percentage |
|  | **Broaden** | | Same as `BROADEN` over the checked lines. The result is named "Envelope". |
| ⤷ Linear Combination and Coefficients | Spectra 1, Spectra 2, Spectra 3 | **1, 1, 1** | Coefficients for the first, second and third checked lines. **The dialog handles at most 3 lines**; the command handles up to 100. |
|  | **Linear Combin** | | Same as `LINECOMBIN` |
| Temporal Analysis Post Processing | **Addition** | | Same as `ADDITION` over the checked lines |
|  | **Subtraction** | | Same as `SUBTRACTION`: first checked line minus the others |
| Post Processing Results | Line number | empty | Reference number that will hold the result. (inferred) If blank, use the next free number. The new line appears **at the end of the list**. |

* **Important difference from PREP:** the calculation buttons **do not write a file**. To save a result, use
  `WRITESPEC` or `WRITETH`.
* The same dialog is used for Transfer Function plots, which are Spectrum plots, and for Time History plots.

### 5.4 Line-object algorithms (the CORE maths behind the buttons and commands)

Notation: the source lines are `L_i = {(x_ij, y_ij)}`, i = 1…n, with n ≤ 100. Each `x_ij` sequence is assumed
strictly increasing; sort it on load if not (inferred).

**Common resampling rule (manual, stated for AVERAGE, SRSS, LINECOMBIN, BROADEN and WRITESPEC):**
1. **Abscissa union:** the result line uses *all* x values of all source lines,
   `X = sort(unique(∪_i {x_ij}))`. Use a tolerance for near-duplicate merging, e.g. relative 1e-9. (inferred)
2. **Interpolation:** inside a line's range, `y_i(x)` is **linear** interpolation in x and y on linear scales,
   even when the axes are logarithmic.
3. **Extrapolation:** outside a line's range, `y_i(x)` is held **constant** at the first or last defined y value.

```
def resample(line, X):           # line: (xs, ys) arrays, xs ascending
    return numpy.interp(X, xs, ys, left=ys[0], right=ys[-1])
```
`numpy.interp` does exactly this.

| Command (default result name) | Syntax | Formula on the union grid X |
|---|---|---|
| `LINECOMBIN` ("Linear Combin.") | `LINECOMBIN,<Dest>,<Line1>,<Coeff1>,…,[Line100],[Coeff100]` | `y(x) = Σ_i c_i · y_i(x)`. Each line needs its own coefficient. |
| `ADDITION` ("Linear Combin.") | `ADDITION,<dest>,<source1>…<source100>` | LINECOMBIN with every `c_i = 1` |
| `SUBTRACTION` ("Linear Combin.") | `SUBTRACTION,<dest>,<source1>…<source100>` | `y = y_1 − Σ_{i≥2} y_i`: LINECOMBIN with c_1 = 1 and all other c_i = −1. Up to 99 lines are subtracted. |
| `AVERAGE` ("Average Line") | `AVERAGE,<dest>,<source1>…<source100>` | `y(x) = (1/n) Σ_i y_i(x)` |
| `SRSS` ("SRSS Line") | `SRSS,<dest>,<source1>…[source100]` | `y(x) = sqrt( Σ_i y_i(x)^2 )` |
| `BROADEN` ("Envelope") | `BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>…<source100>` | Envelope of the sources, then peak bridging (Smooth1 %), then peak broadening (Smooth2 %). See below. |

**Engineering context [CORE]:** these are the standard ISRS post-processing operations. Examples: SRSS of the
X, Y and Z input-direction spectra at one node; a weighted linear combination such as 100-40-40; the average of
several time-history realisations; and the envelope plus broadening of spectra over soil cases (lower bound,
best estimate, upper bound). The manual warns: "A reduced number of frequency steps could affect the accuracy of
the UI spectrum broadening algorithm. Use 301 frequencies as recommended", i.e. 0.1–100 Hz with at least 301
steps (line 6887).

**BROADEN algorithm.** The manual names only the two parameters, "Peak bridging percentage" and "Peak broadening
percentage", and the result name "Envelope". The steps below are a recommended implementation consistent with
ASCE 4 and USNRC RG 1.122 practice. (inferred; see Section 13)
1. **Envelope:** on the union grid X, `E(x) = max_i y_i(x)`, using the resampling rule.
2. **Peak broadening (b = Smooth2/100):** every ordinate at frequency f₀ is spread horizontally over
   `[f₀(1−b), f₀(1+b)]`. Equivalently, `B(f) = max { E(f₀) : f/(1+b) ≤ f₀ ≤ f/(1−b) }`. Evaluate this on a grid
   that includes X plus the shifted points `x(1±b)`, so that plateau edges are kept exactly. With b = 0 nothing
   changes.
3. **Peak bridging (p = Smooth1/100, "Peak Difference(%)"):** find the local maxima of B. For each pair of
   adjacent peaks P_k and P_{k+1}, check the valley between them. If `(min(P_k, P_{k+1}) − valley_min) /
   min(P_k, P_{k+1}) ≤ p`, raise every ordinate between the two peaks to the horizontal line `min(P_k, P_{k+1})`
   wherever it is lower. Repeat until nothing changes. With p = 0 nothing changes.
4. Store the result in `Dest` with the name "Envelope".
The order of steps 2 and 3 is not specified; the order above (broaden first, then bridge) is a proposal. A
single-line BROADEN with both percentages set to 0 returns the line unchanged.

### 5.5 Writing lines [IO][CORE]
* `WRITESPEC,<SpecFile>,<Num1>,…,[Num50]`: writes up to 50 lines into one spectrum-format file. **Every line is
  written on the union of the x values of all the lines**, using the same linear interpolation and constant
  extrapolation as Section 5.4. The layout follows: one frequency column, then one column per line in argument
  order. (inferred)
* `WRITETH,<THFile>,<Num>`: writes a single line in the **one-column** time-history format. The line **must have
  a constant time step**. The time step is taken as `x[1] − x[0]`, the first two points. Recommended: warn if any
  later Δx differs from it by more than a small tolerance. (inferred)
* File ▸ **Export Table** (6.2.3): exports the data of "certain 2D plots" to a **CSV** table for spreadsheets.
  Proposed CSV layout: a header row with x then each line name, then rows on the union grid. (inferred)

### 5.6 Spectrum Plot (7.2.3) — `SPECPLOT,<Line1>,…,<Line50>` [UI]
* Plots lines that use the **spectrum file format**, i.e. x is frequency. The Plot menu calls it **"Spectrum
  TFU-TFI"**.
* **One generic XY axis serves several graph types.** In PREP, the Spectrum, Impedance and TFU-TFI graphs were
  separate plots. In the new UI, the user gets each of these looks by changing the generic spectrum plot with the
  plotting commands: `AXES` for log scales and grids, `XTITLE`, `YTITLE`, `YTITLE2`, `PLOTRANGE`, `PLOTTITLE` and
  `STIPPLE`.
* Command rules:
  * The default extent is the **minimum and maximum over the line objects** being plotted.
  * Reference numbers with no line object are **ignored**.
  * The value **−1 ends the list**; numbers after it are ignored.
  * At most 50 lines.
* Look, from the screenshot (PDF p. 204): a black axis frame; thick horizontal grid lines at the major Y ticks;
  X ticks every 5 Hz on a linear axis; a legend box at the top right with each line's name and a short sample of
  its colour (e.g. `Test2.rs` in red). The axis colours come from the **Spec** palette and the line colours from
  **SpecLines**.
* **Suggested presets (UI convenience, inferred):**
  * RS: log X, frequency in Hz, Y in g.
  * TF amplitude: linear or log X, Y as amplitude ratio.
  * Impedance: X frequency, Y real/imaginary stiffness, with the second Y title for damping.

### 5.7 Time History Plot (7.2.4) — `THPLOT,<Line1>,…,<Line50>` [UI]
* Plots lines in the **time-history file format**, with x as time. It uses the **same axis code** as the
  Spectrum plot, so it looks different from the PREP time-history plot.
* Same command rules as SPECPLOT: data-extent default, unknown numbers ignored, −1 ends the list, at most 50
  lines.
* Look, from the screenshot (PDF p. 205): an acceleration record `Comp.th` drawn in red on a time axis from 0 to
  20 s. Y labels are symmetric (−0.3…0.3). The zero label shows floating-point noise (`-2.7756e-17`), so **round
  tick labels**, e.g. print values smaller than 1e-12 × range as 0. The axis colours come from the **TimeHist**
  palette and the line colours from **THLines**.

### 5.8 Other 2D commands used with these plots (syntax for cross-reference) [UI]
| Command | Syntax | Semantics |
|---|---|---|
| AXES | `AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>` | Each argument is 0 or 1: thick (major) grid lines on X and Y, thin (minor) grid lines on X and Y, log X, log Y |
| PLOTRANGE | `PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>` | Extent of the active 2D plot |
| PLOTTITLE | `PLOTTITLE,<Title>` | Title of 2D and 3D plots. **Ignored by the Soil Layer plot.** |
| XTITLE / YTITLE / YTITLE2 | `XTITLE,<label>` etc. | X-axis title, left Y-axis title, right Y-axis title (the second Y axis is used by the Soil Properties plot) |
| STIPPLE | `STIPPLE,<switch>` | 0 off, 1 on, −1 (default) toggle |
| MARKERS | `MARKERS,<Mark>,<Ln1>…<Ln50>` | Global per-line markers on (1) or off (0) |
| LINENAME | `LINENAME,<Num>,<Name>` | Global line name, used in legends |
| CAPTUREPLOT | `CAPTUREPLOT,<FileName>` | Saves the active plot as **PNG if the file name contains ".png"**, otherwise as **BMP**. Fails if no plot is active. The menu item File ▸ Export Image always writes `.bmp`. |
| CLOSEPLOT | `CLOSEPLOT` | Closes the active plot tab |
| WINDOWSETTINGS | `WINDOWSETTINGS` | Opens the 2D or 3D settings dialog for the active plot |
| LBINCORS | `LBINCORS,<out>,<in>` | Lower-bound incoherent response spectra [ADV]. Specified in the command-reference spec. |

---

## 6. Soil Layer Plot (7.2.5) — `LAYERPLOT`

### 6.1 What it shows [UI][CORE-data]
* A **cross-section of the soil layer column** of the active model, drawn as a vertical stack of rectangles whose
  heights are **proportional to the layer thicknesses**. Layer numbers are printed at the left (1, 2, 3, 4…). The
  fills are shades of grey from the **SoilLayer** palette, set through Options ▸ Colors ▸ "Soil Layers Palette".
* Next to the stack is a **table of layer properties**:

| Column header (screenshot) | Source (`L` command argument) | Shown by default |
|---|---|---|
| Layer | `<nm>` (rows 1…N, then a row labelled **Halfspace**) | always |
| (Thickness) | `<thick>` | **no** (Show Thickness unchecked) |
| Unit Weight | `<weight>` (specific weight) | yes |
| P-Wave Velocity | `<pveloc>` | yes |
| S-Wave Velocity | `<sveloc>` | yes |
| P-Wave Damping Ratio | `<pdamp>` | yes |
| S-Wave Damping Ratio | `<sdamp>` | yes |

* Data source: the **L** command `L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>` and the **TOPL**
  command, which lists the top layers for the SITE module (for embedded models, in the same order as the L
  layers and without repeats).
* `PLOTTITLE` is ignored by this plot.

### 6.2 "Soil Layer Windows Setting" dialog [UI]
The manual's prose and pictures are swapped. The "Controls" paragraph under 7.2.5 describes the soil *properties*
window, and the one under 7.2.6 describes this one. The fields below are transcribed from the screenshot titled
"Soil Layer Windows Setting" (PDF p. 206) together with the 7.2.6 prose that belongs to it.

| Field (exact label) | Default | Meaning |
|---|---|---|
| Start Layer | **1** | First layer to plot. 1 is the minimum and stands for the surface level. |
| EndLayer | **−1** | Last layer to plot. **−1 means find the deepest layer in the model data and use it.** |
| Show Thickness | unchecked | Adds a Thickness column to the table |
| Show Specific Weight | checked | Unit Weight column |
| Show P-Wave Velocity | checked | |
| Show S-Wvae Velocity (sic) | checked | S-Wave Velocity column |
| Show P-Wave Damping Ratio | checked | |
| Show S-Wave Damping Ratio | checked | |
| Ok / Cancel | | |

### 6.3 Implementation notes [UI]
* Rectangle height `h_i = thick_i / Σ_{plotted} thick_j × H_available`.
* The half-space has no thickness. Draw it as a final band of nominal height, or as a hatched base. The
  screenshot shows four layer bands and a table row "Halfspace". (inferred)
* Units follow the model units; nothing is converted. The `L` specific weight is in force per volume units.
  Damping ratios are as entered, i.e. fractions; the screenshot values are dummies.

---

## 7. Soil Properties Plot (7.2.6) — `SOILPROPPLOT,<PropName>`

### 7.1 What it shows [UI][CORE-data]
* For one **dynamic soil property** (a named pair of strain-dependent curves), it plots:
  * **Damping ratio vs. shear strain**: red line with point markers, on the **left Y axis** titled "Damping Ratio";
  * **Shear modulus (G/Gmax, modulus reduction) vs. shear strain**: green line with point markers, on the
    **right Y axis** titled "Shear Modulus", running from 0 to 1;
  * X axis titled **"Shear Strain %"**.
* Legend entries (exact text): `Damping Ratio/Shear Strain` (red) and `Shear Modulus/Shear Strain` (green).
* The axis and line colours come from the **SoilProp** palette.

### 7.2 Opening dialog "Select Dynamic Soil Property" [UI][IO]
The dialog appears when the menu item is used, or when `SOILPROPPLOT` has no name. The name is
**case-sensitive** and must match a property in memory exactly, or no plot appears.

| Element | Behaviour |
|---|---|
| Property list (top left) | Names of the properties defined in the model, e.g. `Clay`, `Rock`, `Sand`. Select one to plot or edit it. |
| **New** | Opens a popup that asks for a new property name. On OK the name is added to the list. |
| **Edit** | Renames the selected property |
| **Delete** | Removes the selected property from the list |
| **Ok** / **Cancel** | Ok plots the selected property |
| Data table (bottom) | Rows 1…n. **Four columns: `Strain` · `Mod. Red.` · `Strain` · `Damp`.** The modulus-reduction curve and the damping curve each have **their own strain column**, so their strain points may differ. Cells can be edited in place. |
| Title | Edit box showing the property name or title, e.g. `Rock` |

Example data from the screenshot for "Rock":

| # | Strain | Mod. Red. | Strain | Damp |
|---|---|---|---|---|
| 1 | 0.0001 | 1 | 0.0001 | 0.4 |
| 2 | 0.001 | 1 | 0.0003 | 0.8 |
| 3 | 0.01 | 0.9875 | 0.001 | 1.5 |
| 4 | 0.1 | 0.9525 | 0.003 | 3 |
| 5 | 1 | 0.90 | 0.01 | 4.6 |
| 6 | 10 | 0.81 | 0.03 | 6 |
| … | (rows continue; the table scrolls) | | | |

The values imply **strain in percent** (the axis says "Shear Strain %") and **damping in percent**, which is the
SHAKE convention. The plot screenshot, however, shows a left axis from 0.2 to 0.9. See Section 13.

**Limits [CORE]**, from manual Section 1.2 for the SOIL module: at most **100 soil material curves** and at most
**11 data points per soil curve**. The table should enforce 11 rows. (inferred mapping)

### 7.3 "Soil Properties Window Settings" dialog [UI]
Transcribed from the screenshot (PDF p. 208). The matching prose is under 7.2.5 "Controls": the user can choose
log axes for X or Y, show minor ticks, and hide or show the 2 lines.

| Group | Field (exact label) | Default |
|---|---|---|
| G,D Axis | Logarithmic | unchecked |
|  | Show Ticks | unchecked |
| Shear Strain Axis | Logarithmic | unchecked |
|  | Show Ticks | unchecked |
| (bottom) | Show Shear Modulus Line | checked |
| (bottom) | Show Daiming Line (sic, "Damping") | checked |
| buttons | Ok, Cancel | |

**Recommendation (inferred):** default the shear-strain axis to log. Strain data spans 1e-4 to 10 %, and on a
linear axis the low-strain curve collapses. Keep the checkbox so users can switch.

---

## 8. Animated plots (7.3)

### 8.1 Data flow overview [CORE][IO]

```
SSI module run with a post-processing "Restart" option
      │  (MOTION → \TFU \ACC \ACCR \RS ; STRESS → \NSTRESS \SOILPRES ; RELDISP → \THD \THDR)
      ▼
ASCII frame files  +  an animation frame LIST file (*.dispani, *.tfiani, *.impani; legacy zpani, thiani, contani, thani)
      │
      ▼  Plot ▸ Process Animation Frame List   (PROCFRAME)      ── or ──  binary DB:  LOAD*DB + ACCDBANI/DISPDBANI/THSDBANI
      ▼
binary GPU-ready frames in a "Frame Storage Dir"  +  an entry in SASSIani.xml
      │
      ▼  Plot ▸ Bubble | Vector | Contour | Deformed Shape   (Load Frame Data dialog, or the command with args)
      ▼
animation in a 3D plot tab (Pause, +/-, Frame Pause, Scale Factor, Colormap range)
```

The SSI modules cannot write the binary layout the graphics card needs. The UI therefore converts the frames
**once** in a separate step that can take a long time; the user should do it only once per data set.

### 8.2 Which result frames feed which animation [CORE]
The table below combines Table 3.2 (PDF p. 50) and the module post-processing option descriptions.

| Producer and option | Sub-directory | Frame naming scheme (Table 3.2) | Content | Intended plot |
|---|---|---|---|---|
| MOTION: Restart for TF | `\TFU` | `TFU_freq_filenum`, e.g. `\TFU\TFU_000.02_00001` (freq = frequency, fnum = frame number) | **Complex** acceleration TF at all active nodal DOFs at one frequency | **Vector** (animated over frequency) |
| MOTION: Restart for ACC | `\ACC` (and `\ACCR` for rotations) | `ACC_time_filenum`, e.g. `\ACC\ACC_00.000_00001` | Acceleration at all nodes at one time step | **Deformed Shape** (animated) or **Bubble** (static) |
| MOTION (always written with the ACC frames) | — | `ACC_max.txt` | Maximum acceleration (ZPA) frame | **Bubble** (static ZPA; BUBBLEPLOT is described as "Plot Bubble (ZPA)") |
| MOTION: Restart for RS | `\RS` | `RS##_freq_filenum`, e.g. `\RS\RS01_000.10_00001` (## = damping number) | RS values at all nodes at one frequency | **Deformed Shape** (animated) or **Bubble** (static) |
| STRESS: Restart for Nodal Stress Contours | `\NSTRESS` | `stress_time_fnum_comp`, e.g. `\NSTRESS\stress_00.000_00001_sig`; maxima `stress_ABS_MAX_comp`, e.g. `stress_ABS_MAX_sig` | **Averaged nodal** stress per component. **sig**: solids normal stress, shells membrane stress. **tau**: solids shear stress, shells membrane shear. **bdsig**: bending stress (shells only). **bdtau**: bending shear (shells only). | **Contour** (static for one time or the maximum, or animated) |
| STRESS: Restart for Soil Pressure Contours | `\SOILPRES` | `press_time_fnum_type`, e.g. `\SOILPRES\pres_00.000_00001_nod`; maxima `pres_ABS_MAX_type`. type = `ele` (element values) or `nod` (nodal values). | Seismic, or total (static + seismic), soil pressure on the walls and mat | **Contour** |
| RELDISP: Restart for Frame Generation | `\THD` (and `\THDR`) | `THD_time_filenum`, e.g. `\THD\THD_00.000_00001` | Relative displacement at all nodes at one time step | **Deformed Shape** (also contour) |
| Binary DBs [ADV] | user dir | — | `Modelname_ACC.bin`, `Modelname_THD.bin` (displacement), `Modelname_STRESS.bin` | ACC DB, RelDisp DB and Stress DB animations |

**Warnings to carry over [CORE]:**
* STRESS frame files contain **average nodal stresses and pressures for plotting only**. They were found by
  setting each node equal to the centre stress of its element (no shape functions) and then averaging. Design
  values are the element-centre values (`ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`, `pres_max_ele`).
  The UI should print this caveat on stress and pressure contour plots, for example in the title bar or status
  line. (inferred UI practice)
* **No frame files are written for TSHELL elements.**
* STRESS frame generation covers only SOLID and SHELL elements in 3D models. In mixed SOLID/SHELL models only the
  membrane averages are framed; bending frames (`bd*`) are written for shells only.
* Before a STRESS or soil-pressure restart, the frames to generate are listed in `Frames.txt`. The format is
  `[# of Frames]`, the frame numbers one per line, `[# of Soil Pressure Groups]`, then the group numbers on one
  line. (See spec `05d`.)

### 8.3 Frame file and list-file conventions [IO]
These conventions are scattered through the manual:
* An **ASCII frame file** has a **header that states the number of rows and columns**. **The second number in the
  header is the number of columns** (MODFRAMES rewrites only that number). Frames written by earlier versions
  have a "legacy header" and must be converted with `MODFRAMES,<cols>,<framelist>` before FRAMECOMBIN can use
  them.
* Proposed concrete layout (inferred): line 1 `nrows ncols`; then `nrows` lines of
  `node_id v_1 … v_(ncols-1)`, or `ncols` values if the node id is not counted as a column. See Section 13.
  Columns hold DOF components, e.g. X, Y, Z, or Re/Im pairs for TF frames. The `Col` argument of BUBBLEPLOT and
  CONTOURPLOT selects one of them.
* An **animation frame list file** has extension `*.dispani`, `*.tfiani` or `*.impani`, plus the legacy
  extensions `zpani`, `thiani`, `contani` and `thani` named in MODFRAMES. Its format is the PREP frame-listing
  format. **The UI ignores the first (header) line.** The remaining lines list the frame files in animation
  order. Inferred layout: one full path per line, or a path relative to the list file's folder.
* `FRAMECOMBIN,<op>,<num>,<InFile1>,…,<InFileX>,<Outfile>` combines frames node by node and column by column:
  op 0 = SRSS, 1 = sum, 2 = average. Use it, for example, to SRSS the X, Y and Z input-direction stress frames
  before plotting.
* `FRAMESEL,<tol>,<Acc>,<Var>` finds the local maxima and minima of an acceleration history that are above `tol`
  % of the global maximum, so you can tell which animation frames are critical.
* `BINFRAMEOUT,<db>,<frame>,<TS>,[Split],<dir>` writes ASCII frames in the UI frame format from a loaded binary
  DB. `frame = -1` with `TS ≤ 0` writes the frame of maximum components. `MAXDBFRAME,<Type>,[dir]` builds a
  one-frame animation of the global maxima.

### 8.4 "Process Animation Frame List" (7.3.1) — dialog "Parse Frame Data", command `PROCFRAME`

| Field (exact label) | Type | Example (screenshot) | Required | Meaning |
|---|---|---|---|---|
| List File Name | edit + `<<` browse | (empty) | **yes** | Path of the animation frame list file (Section 8.3) |
| Frame Storage Dir | edit + `<<` browse | `C:/test/tshells/stressframes/` | **yes** | Directory that will receive the converted binary frames (the "buffer directory") |
| Data Description | edit | `stress animation for the tshell example` | informational; recommended | Label stored in SASSIani.xml and shown in the Load Frame Data list |
| Plot Type | radio: **Bubble**, **Vector**, **Contour**, **Time History**, **Stress DB (Binary)**, **ACC DB (Binary)**, **RelDisp DB (Binary)** | Stress DB (Binary) selected | informational; recommended | Tag stored with the animation. The three "(Binary)" items are for animations built from the SSI response **binary databases** [ADV]. |
| Ok / Cancel | | | | Ok starts processing |

Command: `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]`
* If `AniFile` or `BufferDIR` is blank, the dialog opens.
* `Data`: the label stored in the database.
* `anitype`: an integer tag **for sorting the database only**, which does not change how frames are processed.
  0 = Bubble (default), 1 = Vector, 2 = Contour, 3 = Time History. The binary types have no command code; they are
  created by ACCDBANI, DISPDBANI and THSDBANI instead.

Processing behaviour:
1. Read the list file and skip line 1.
2. For each frame file, parse it and convert it to the binary GPU layout in the Frame Storage Dir. Show a
   **progress bar in the lower right of the window** with the percentage of frames converted. **The UI accepts no
   other input until processing finishes.** Run it as a modal background task.
3. Append an entry to **SASSIani.xml** with the location and information of the animation: description,
   directory, type, number of frames, and default frame range and scale/colour range.
4. Users should not edit the binary frame files, only delete them. Use any internal binary format you like, e.g.
   one float32 array per frame plus an index file. (inferred)

**Binary-database route [ADV]** (Section 9.18 commands): load one DB with `LOADACCDB,<file>`,
`LOADDISPDB,<file>` or `LOADTHSDB,<file>,[sel]` (sel 0 = ACS SASSI format, 1 = ANSYS format). Then create the
animation with `ACCDBANI,<dir>,[label]`, `DISPDBANI,<dir>,[label]` or `THSDBANI,<dir>,[label]`. The label is
stored in the animation database. Only one DB of each kind can be in memory at a time.

### 8.5 Opening animations (7.3.2) — dialog "Load Frame Data"
This dialog appears the first time any animated plot is chosen, and when an animation command has no arguments.

| Group | Element | Example / default | Behaviour |
|---|---|---|---|
| Select From Database | List with columns **Description · Animation Directory · Type · Frames** | `TestAnimation Demo 1 Model · C:\SSI Files\Submodeler … · Displa… · 3110` | Built from SASSIani.xml, which sits in a system-dependent default location so users never have to look for it. The user **must** click a row. |
|  | **Remove Animation** | | Deletes the entry from the list **and deletes the data files** that processing produced |
| Animation Control ▸ Frame Selection | Start | **1** | First frame (1-based) |
|  | End | **number of frames** (3110 in the example) | Last frame |
|  | Stride | **1** | Show every k-th frame |
| Animation Control ▸ Data Scale | Scale Factor | **1.000000** | For Vector and Deformed Shape plots |
|  | (Colormap range Min / Max) | data min and max | For Bubble and Contour plots. Depending on the plot type, the dialog shows either the scale factor or the colour-map limits. |
| buttons | Ok, Cancel | | |

When a row is selected, the controls fill with the defaults stored in SASSIani.xml for that animation.

**Frame sequence [CORE]:** the frames shown are `k = Start, Start+Stride, …` up to and including End (if End is
reached exactly). Count = `floor((End − Start)/Stride) + 1`. Validate `1 ≤ Start ≤ End ≤ Frames` and
`Stride ≥ 1`. (inferred validation)

Command-argument equivalents: `<BufferDir>` is the processed frame storage directory, which identifies the
animation; `<MnF>` = Start; `<MxF>` = End; `<ST>` = Stride; `<MnR>`/`<MxR>` = colour-map minimum and maximum
data values; `<Col>` = data column; `<Scale>` = scale factor.

### 8.6 Animation playback controls (7.3.3 → 7.1.1) [UI]
* The view controls are those of Section 1.
* `PAUSE,[pz]`: −1 toggles (default), 0 starts, 1 stops. The Pause key and the toolbar button do the same.
* `+` / `-` step one frame, for use while paused.
* The frame period is **Frame Pause (ms)**, default 33. Animations loop, wrapping from the last frame back to the
  first. (inferred)
* `DEBUG` shows the animation information: current frame number, frame time or frequency if known, and the view
  parameters. (inferred content)
* Window Options: **Animation Options ▸ Scale Factor** (Vector and Deformed) and **Colormap Value Range**
  (Bubble and Contour) can be changed during playback. **Title** adds a title to the animation plot.

### 8.7 Bubble Plot (7.3.4) — `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>`
* Shows an animation of **all used nodes** of the model as **bubbles** (screen-aligned discs or spheres). Each
  bubble's **size and colour** depend on that node's value in the current frame, read against a **colour bar on
  the right**.
* **Unused nodes are not drawn.** These are:
  1. nodes not connected to any element;
  2. orientation nodes such as beam **K-nodes**;
  3. nodes on elements that were simplified after conversion from an outside program.
* **Colour bar** (from the screenshot): a vertical bar running from dark red at the top (max) through red,
  yellow, green and cyan to blue and dark blue at the bottom (min), i.e. a "jet"-type map. **Five labels** spaced
  evenly from Max to Min, printed with 5 decimals (e.g. 0.43260, 0.20515, −0.02229, −0.24974, −0.47718).
* Value-to-colour mapping [CORE-for-verification]: `t = clamp((v − MnR)/(MxR − MnR), 0, 1)`, colour = jet(t).
  Values outside the range take the end colours (dark blue below, dark red above).
* Value-to-size mapping (inferred): `size = S_max · max(t_min, t)` with `S_max` = Node/Bubble Node Size (default
  10). Alternatively scale by `|v| / max(|MnR|, |MxR|)` so that zero gives the smallest bubble. The manual says
  only that size is "associated with the data". See Section 13.
* The colours of the plot elements come from the **Bubble** palette.
* Typical data: ZPA (`ACC_max.txt`), RS amplitude at one frequency, acceleration at one time. Pick the component
  with `Col`, e.g. X, Y or Z.

### 8.8 Vector Plot (7.3.5) — `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>`
* Shows an animation of the **direction and magnitude** of the response vector at each **used** node (the same
  exclusions as the Bubble plot). The manual says "force"; the data is usually a complex acceleration TF.
* **Each node has 3 vectors: red for the X-direction data, green for Y and blue for Z.** The **Output Direction**
  radio buttons (X, Y, Z, All) choose which ones are drawn.
* **Complex-number drawing rule [CORE-for-verification] (manual):**
  * With no imaginary part, each vector lies along its own axis.
  * With an imaginary part, the vector tilts out of its axis by the same amount in **the other two axes**.
    Manual example: an X vector of 5 + 0.3i is drawn as x = 5, y = 0.3, z = 0.3.
  * Generalisation (inferred): for nodal complex components `u_X = a_X + i b_X`, `u_Y = a_Y + i b_Y` and
    `u_Z = a_Z + i b_Z`, the drawn vectors from node position **p** are
    ```
    V_X = Scale · ( a_X , b_X , b_X )   # red
    V_Y = Scale · ( b_Y , a_Y , b_Y )   # green
    V_Z = Scale · ( b_Z , b_Z , a_Z )   # blue
    ```
    with Scale from the dialog, Window Options or Shader Options (default 1.0). Each vector runs from p to
    p + V. Draw the node as a black dot; the screenshot shows dots with red, green and blue segments.
* The colours of the plot elements come from the **Vector** palette.

### 8.9 Contour Plot (7.3.6) — `CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>`
* Shows **contour data on the actual element faces**. The data is **nodal** (one value per node per frame). The
  colour is **linearly interpolated across each face** from the vertex values (Gouraud-type: per-vertex colour, or
  per-vertex value interpolated and then colour-mapped per fragment). Element outlines are drawn with the outline
  thickness; the screenshot shows a green disc with dark outlines.
* The colour map, colour bar and clamping are the same as the Bubble plot. The range is `MnR`/`MxR` or the
  Colormap Value Range field; the column is `Col`.
* **Recommended (inferred):** interpolate the *scalar* and map it to colour per fragment, not interpolate RGB.
  This gives true linear contours, and true linear contours are what make a contour plot verifiable.
* The colours of the plot elements come from the **Contour** palette.
* Typical data: `\NSTRESS` nodal stress (sig, tau, bdsig, bdtau), static (one time or ABS_MAX) or animated, and
  `\SOILPRES` nodal pressures (`_nod`).

### 8.10 Deformed Shape Plot (7.3.7) — `DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>`
* Animates the model shape from **displacement data at every node** for one time per frame. The tab caption is
  "Displacement Plot".
* Deformed coordinates [CORE-for-verification]: `x' = x + Scale·u_x`, `y' = y + Scale·u_y`,
  `z' = z + Scale·u_z`. Here `Scale` is the scale factor, "Scalar to increase deformation for animation"
  (default 1). The frame supplies `u` per node (relative displacements from RELDISP `\THD`; ACC or RS frames give
  "deformed shapes" of those quantities). Nodes missing from a frame are treated as having zero displacement.
  (inferred)
* **Show Undeformed Shape** (Window Options check box): overlays a **wireframe of the original model**. The
  screenshot shows a grey deformed model with red outlines and a light undeformed wireframe.
* Faces are drawn with the outline thickness. The colours come from the **Deformed** palette.

---

## 9. Data-requirements matrix (what each plot needs) [CORE]

| Plot | Model geometry | Element attributes | BCs / masses | Line objects | Soil data | Frames / DB |
|---|---|---|---|---|---|---|
| Element | nodes, elements | group, type, material, property | fixed DOFs, masses (optional markers) | — | — | — |
| Node | nodes, connectivity | — | interaction nodes, fixed DOFs, masses, selection | — | — | — |
| Cut | nodes, elements of the model given by number | cut element list | — | — | — | — |
| Spectrum | — | — | — | ≥ 1 line (x = frequency) | — | — |
| Time History | — | — | — | ≥ 1 line (x = time) | — | — |
| Soil Layer | — | — | — | — | L layers (thick, weight, Vp, Vs, Dp, Ds), half-space, TOPL | — |
| Soil Properties | — | — | — | — | the named dynamic property (G/Gmax–γ, D–γ) | — |
| Bubble | nodes, connectivity (to find used nodes; K-nodes excluded) | — | — | — | — | scalar per node per frame (column `Col`) |
| Vector | used nodes | — | — | — | — | X, Y, Z components per node, real or complex |
| Contour | nodes, element faces | — | — | — | — | scalar per node per frame |
| Deformed | nodes, element faces/edges | — | — | — | — | u_x, u_y, u_z per node per frame |

---

## 10. Command quick-reference for all plotting commands [IO][UI]

This is a cross-reference. The full argument descriptions are in Section 9.14 of the manual and in the
command-reference spec. Command names are case-insensitive; the minimum abbreviation rules are given in
Chapter 9.

| Command | Syntax | Purpose |
|---|---|---|
| ADDITION | `ADDITION,<dest>,<source1>…<source100>` | Sum of lines |
| AVERAGE | `AVERAGE,<dest>,<source1>…<source100>` | Average of lines |
| AXES | `AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>` | 2D grid lines and log axes |
| BROADEN | `BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>…<source100>` | Envelope with peak bridging and broadening |
| BUBBLEPLOT | `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | Bubble animation |
| CAPTUREPLOT | `CAPTUREPLOT,<FileName>` | Save an image (PNG if ".png" is in the name, else BMP) |
| CLOSEPLOT | `CLOSEPLOT` | Close the active plot |
| CNGCENTER | `CNGCENTER,<X>,<Y>,<Z>` | Set the 3D rotation centre |
| CNGVIEW | `CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>` | Set the 3D view, e.g. to show two plots from the same viewpoint |
| COLOR | `COLOR,<Palette>,<Num>,<R>,<G>,<B>` | Edit a palette colour (0–255) |
| CONTOURPLOT | `CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | Contour animation |
| CUTPLOT | `CUTPLOT,[Cut],[Model]` | Cut plot |
| DEBUG | `DEBUG,[switch]` | 0 off, 1 on, 2 toggle (default) |
| DEFORMPLOT | `DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | Deformed-shape animation |
| ELECOLOR | `ELECOLOR,<val>` | 1 group, 2 material, 3 property |
| ELENUM | `ELENUM,[opt]` | −1 toggle, 0 off, 1 on; turns GROUPNUM off |
| GROUPNUM | `GROUPNUM,[opt]` | −1 toggle, 0 off, 1 on; turns ELENUM off |
| LAYERPLOT | `LAYERPLOT` | Soil layer plot |
| LBINCORS | `LBINCORS,<out>,<in>` | Lower-bound incoherent RS [ADV] |
| LINECOMBIN | `LINECOMBIN,<Dest>,<Line1>,<Coeff1>…[Line100],[Coeff100]` | Linear combination |
| LINENAME | `LINENAME,<Num>,<Name>` | Rename a line (global) |
| MARKERS | `MARKERS,<Mark>,<Ln1>…<Ln50>` | Point markers (global) |
| MODELPLOT | `MODELPLOT` | Element plot |
| NODENUM | `NODENUM,[opt]` | −1 toggle, 0 off, 1 on |
| NODEPLOT | `NODEPLOT` | Node plot |
| NODESEL | `NODESEL,<N1>…<N20>` | Toggle node selection (accumulates) |
| PAUSE | `PAUSE,[pz]` | −1 toggle, 0 start, 1 stop |
| PLOTRANGE | `PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>` | 2D extent |
| PLOTTITLE | `PLOTTITLE,<Title>` | Title (ignored by the soil layer plot) |
| PROCFRAME | `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]` | Process frames (anitype 0 bubble, 1 vector, 2 contour, 3 time history) |
| READSPEC | `READSPEC,<SpecFile>,<numLines>,<Line1>…<LineN>` | Load spectrum columns |
| READTH | `READTH,<THFile>,<Pair>,<Num>` | Load a time history (Pair 0 one column, 1 pairs) |
| RSTCENTER | `RSTCENTER` | Reset the rotation centre to the bounding-box centre |
| RSTVIEW | `RSTVIEW` | Reset the view and zoom |
| SHADEROPTIONS | `SHADEROPTIONS,[points],[linew],[shrink],[scale]` | Shader options |
| SHOWDOF | `SHOWDOF,[label1],…,[label6]` | X, Y, Z, XX, YY, ZZ, DISP, ROT, ALL |
| SHOWMASS | `SHOWMASS,[opt]` | −1 toggle, 0 off, 1 on |
| SHRINK | `SHRINK,[switch]` | 1 on, 0 off, −1 toggle |
| SOILPROPPLOT | `SOILPROPPLOT,<PropName>` | Soil property plot (name is case-sensitive) |
| SPECPLOT | `SPECPLOT,<Line1>,…,<Line50>` | Spectrum plot (−1 ends the list) |
| SRSS | `SRSS,<dest>,<source1>…[source100]` | SRSS of lines |
| STIPPLE | `STIPPLE,<switch>` | 0 off, 1 on, −1 toggle |
| SUBTRACTION | `SUBTRACTION,<dest>,<source1>…<source100>` | First line minus the rest |
| THPLOT | `THPLOT,<Line1>,…,<Line50>` | Time-history plot (−1 ends the list) |
| VECTORPLOT | `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | Vector animation |
| WINDOWSETTINGS | `WINDOWSETTINGS` | Open the settings dialog for the active plot |
| WIREFRAME | `WIREFRAME,<Switch>` | 0 off, 1 on, −1 toggle |
| WRITESPEC | `WRITESPEC,<SpecFile>,<Num1>,…,[Num50]` | Write lines on the union grid |
| WRITETH | `WRITETH,<THFile>,<Num>` | Write a one-column TH; dt = x[1] − x[0] |
| XTITLE / YTITLE / YTITLE2 | `XTITLE,<label>` … | Axis titles (YTITLE2 is the right axis) |

---

## 11. Toolbar reference (Chapter 8)

Both toolbars can be hidden from **View ▸ Toolbars**. They are visible by default, and this setting is **not**
saved in SASSIini.xml (spec `05d`). In the screenshots, the main toolbar is the first row and the 3D plot
toolbar is the second row, with separators grouping the buttons.

### 11.1 Main Toolbar (8.1) [UI]

| # | Icon (visual description) | Description (exact) | Menu option (manual section → menu item) | Command (manual section → command) |
|---|---|---|---|---|
| 1 | blank white page | Create a New Model | 6.1.1 Model ▸ New | — |
| 2 | yellow folder or bin with an up arrow | Open the Model Database | 6.1.2 Model ▸ Open An Existing Model | — (the RESUME command loads binary model data; it is not listed in the table) |
| 3 | yellow folder or bin with a down arrow | Save a Model | 6.1.3 Model ▸ Save | 9.2.29 `SAVE` |
| 4 | funnel/output icon with a small page | Output model to .pre | 6.1.8 Model ▸ Output | 9.2.47 `WRITE,[<file>],[<path>]` |
| 5 | funnel/output icon (ANSYS) | output model to ANSYS | 6.1.9 Model ▸ Export to ANSYS | 9.11.1 `ANSYS,[FileName],[Dir]` |
| 6 | letter **H** | use .hou converter | 6.1.5 Converter for SASSI Model Fixed Format Input | 9.11.4 `CONVERT,SSI,<model>,<filename>` |
| 7 | letter **A** | use ANSYS® converter | 6.1.6 Converter from ANSYS .cdb File | 9.11.4 `CONVERT,ANSYS,<model>,<filename>,<gravity>` |
| 8 | letter **S** | use GT STRUDL converter | 6.1.7 GT-STRUDL Database Converter (**Not Included**) | 9.11.4 `CONVERT,STRUDL,…` (not applicable in this version; disable it) |
| 9 | camera | capture a plot image | 6.2.2 File ▸ Export Image from a Plot | 9.14.6 `CAPTUREPLOT` |
| 10 | yellow 3D cube | Element Plot | 7.1.2 | 9.14.23 `MODELPLOT` |
| 11 | dotted grid of points | Node plot | 7.1.3 | 9.14.25 `NODEPLOT` |
| 12 | orange marker with a cut | Cut Plot | 7.1.4 | 9.14.12 `CUTPLOT` |
| 13 | grey triangle with a rainbow | Spectrum Plot | 7.2.3 | 9.14.40 `SPECPLOT` |
| 14 | red clock | Time History Plot | 7.2.4 | 9.14.44 `THPLOT` |
| 15 | stacked grey layers | Soil Layer Plot | 7.2.5 | 9.14.19 `LAYERPLOT` |
| 16 | brown soil wedge | Soil Properties Plot | 7.2.6 | 9.14.39 `SOILPROPPLOT` |
| 17 | film reel | Process Animation Frames | 7.3.1 | 9.14.30 `PROCFRAME` |
| 18 | red, green and blue dots | Bubble Plot | 7.3.4 | 9.14.5 `BUBBLEPLOT` |
| 19 | black dot with an arrow | Vector Plot | 7.3.5 | 9.14.45 `VECTORPLOT` |
| 20 | concentric red/yellow/blue contour | Contour Plot | 7.3.6 | 9.14.11 `CONTOURPLOT` |
| 21 | small pin with a deflected tip | Deformed Shape Plot | 7.3.7 | 9.14.14 `DEFORMPLOT` |
| 22 | tripod/gear with a yellow mark | Change Analysis Options | "6.5.5" in the manual. The Analysis Options dialog is actually **6.5.4**; 6.5.5 is Windows Settings. | — |
| 23 | tripod/gear | Afwrite current model | — | 9.2.3 `AFWRITE` (runs CHECK first and writes the input files of all requested modules) |

### 11.2 3D Plot Toolbar (8.2) [UI]

| # | Icon | Description (exact) | Command (manual section → command and argument) |
|---|---|---|---|
| 1 | magnifier | Change View of Current Plot | 9.14.9 `CNGVIEW` (opens an input for the six view values; inferred) |
| 2 | magnifier with a blue X | Reset current plot view to default | 9.14.34 `RSTVIEW` |
| 3 | target / cross-hair | Change center of current plot | 9.14.8 `CNGCENTER` |
| 4 | target with a blue X | Reset Center of current plot to default | 9.14.33 `RSTCENTER` |
| 5 | overlapping squares | Wireframe view | 9.14.47 `WIREFRAME` (toggle, −1) |
| 6 | red/blue shrunk blocks | Shrink faces of elements | 9.14.38 `SHRINK` (toggle, −1) |
| 7 | "G" over coloured dots | Colors of elements determined by group | 9.14.15 `ELECOLOR,1` |
| 8 | "M" over coloured dots | Colors of elements determined by material | 9.14.15 `ELECOLOR,2` |
| 9 | "P" over coloured dots | Colors of elements determined by property | 9.14.15 `ELECOLOR,3` |
| 10 | "N1" over a node | Show node labels | 9.14.24 `NODENUM` (toggle) |
| 11 | blue/orange "E" element | Show element labels | 9.14.16 `ELENUM` (toggle) |
| 12 | blue/red "G" element | Show group labels | 9.14.17 `GROUPNUM` (toggle) |
| 13 | support symbol with a check mark | show fixed degrees of freedom | 9.14.36 `SHOWDOF` (opens the DOF check-box window) |
| 14 | node with a star/mass and a check mark | show lumped masses | 9.14.37 `SHOWMASS` (toggle) |
| 15 | red/green play/pause | Pause Start Animation | 9.14.27 `PAUSE` (toggle) |
| 16 | brown bug | Show Debug info | 9.14.13 `DEBUG` (toggle, 2) |

Toolbar rules (inferred): a toggle button shows its pressed/checked state. Buttons that do not apply to the active
plot are disabled, e.g. Shrink and Wireframe are valid "only on the Element Plot" and Pause only on animations.
A toolbar button sends the equivalent command text to the command processor, so the Command History logs it and
it can be replayed in a `.pre` file.

---

## 12. Verification checklist (to make plotting "verifiable") [CORE]

These are suggested acceptance tests. Each has a closed-form expected result.

1. **Resampling:** L1 = {(1,1),(3,3)} and L2 = {(2,10),(4,20)}. `ADDITION,3,1,2` gives the grid {1,2,3,4}.
   L1 at x = 4 is 3 (constant extrapolation) and L2 at x = 1 is 10. Expected y = {11, 12, 18, 23}: at x = 2,
   L1 = 2 and L2 = 10; at x = 3, L1 = 3 and L2 = 15; at x = 4, 3 + 20 = 23.
2. **SRSS:** three identical lines y = 1 give √3 at every x. AVERAGE of y = 1 and y = 3 gives 2.
3. **SUBTRACTION** with three lines: y1 − y2 − y3. **LINECOMBIN** 1.0, 0.4, 0.4 checks the 100-40-40 rule.
4. **BROADEN** of one triangular spike at f₀ = 10 Hz with Smooth2 = 15 gives a plateau from 8.5 to 11.5 Hz at the
   peak value. With Smooth1 = Smooth2 = 0 it returns the envelope of the inputs.
5. **WRITESPEC/READSPEC round trip:** writing and reading back reproduces the values on the union grid exactly.
   **WRITETH** writes a line with x = {0, 0.01, 0.02,…} with dt = 0.01. A non-uniform line triggers a warning.
6. **READSPEC column counting:** a 1 + 3 column file with numLines = 3 creates 3 lines, and the frequency column
   is not counted.
7. **Hide/Show parsing:** `5 9 12` hides 3 elements and `10-20` hides 11. `5 10-20` is rejected, or handled
   according to the policy chosen in Section 13.
8. **Colour map:** a value of MnR − 1 is dark blue, MxR + 1 is dark red, and the midpoint value gives jet(0.5),
   i.e. green. Five colour-bar labels are evenly spaced between Max and Min.
9. **Vector rule:** a node with X component 5 + 0.3i and Scale 1 gives a red vector (5, 0.3, 0.3).
10. **Deformed shape:** with Scale 10 and u = (0.01, 0, 0), the node moves by 0.1 in X.
11. **Frame selection:** Start 1, End 3110, Stride 1 gives 3110 frames. Start 1, End 10, Stride 3 gives frames
    {1, 4, 7, 10}.
12. **ElemPalette:** group 130 uses palette colour 130 mod 128 = 2.
13. **SPECPLOT,1,2,-1,3** plots lines 1 and 2 only. SPECPLOT with an undefined line number ignores it.
14. **Soil layer EndLayer −1** with 7 defined layers plots layers 1–7 plus the half-space row.

---

## 13. Open questions / ambiguities

1. **File formats of the "spectrum" and "time history" line files.** The manual never gives exact layouts.
   * Unknown: header lines (count and contents), comment conventions, delimiters, and whether `.RS`, `.TFU`,
     `.TFI` and `.TFD` files have the same layout.
   * Proposed:
     * spectrum format: optional header lines (skip non-numeric lines), then rows of `x y1 [y2 …]`;
     * TFU/TFI/TFD: frequency followed by Re and Im, and possibly amplitude;
     * one-column TH (READTH Pair = 0, the MOTION `.ACC` output): a header giving the number of points and dt,
       then the values. WRITETH computes a dt, which suggests the dt is stored somewhere.
   * Check these against real MOTION output before freezing them.
2. **What a "line" is for complex TF files.** Reading the columns in order would make Re and Im separate lines.
   It is unclear whether the UI computes the amplitude |TF| itself. Decide whether to add an "amplitude"
   convenience.
3. **Name of a line loaded from a multi-column file.** The screenshot shows `1: Test2.rs`. The name given to
   columns 2…n is unknown. Proposal: `<filename>[col k]`.
4. **"Post Processing Results – Line number"** in Graph Plot Options. It is not described. It is assumed to be
   the destination reference number, with the next free number when blank, since the text says the result
   appears at the end of the list.
5. **BROADEN algorithm details.** The manual names only Smooth1 ("Peak bridging percentage", dialog "Peak
   Difference(%)") and Smooth2 ("Peak broadening percentage", "Broaden(%)"). Unknown: the order of operations,
   the exact bridging criterion, whether operations work on linear or log frequency, and whether every source is
   broadened. The proposal is in Section 5.4.
6. **"Data Points (highlighted line)"**: how a line becomes "highlighted", by selecting it in the list or by
   clicking it on the plot, is not described.
7. **"Show Ticks"**: assumed to be minor ticks (AXES MinTick). It might instead switch all tick labels.
8. **Units and convention of the Soil Properties table.** Strain is labelled "Shear Strain %", but the damping
   example values (0.4…6) look like percent while the plot axis shows 0.2–0.9. Decide on percent damping, as in
   SOIL and SHAKE input, and label the axis accordingly. Also unknown: whether the 11-point limit applies to the
   UI table.
9. **Default 3D view.** RSTVIEW says "top down", but the opening screenshots look oblique. Also unknown: the units
   of CNGVIEW rotations (degrees assumed), whether `zoom` is a scale factor, and the per-key rotation step.
10. **Model Display Volume clipping rule:** does an element count as outside when any of its nodes is outside, or
    only when all are? Is clipping per vertex (GPU clip planes) or per element?
11. **Hide/Show input:** what happens with mixed lists and ranges, or with several groups? The manual says "only
    a list or a single range" and "only one group". Are hide requests saved with the model by SAVE or WRITE
    (`.pre`)? They are "stored in the model data", but this is not stated.
12. **Group-list label format** `Group 1 SHELL - 70`: the meaning of the number is not stated. It is assumed to be
    the element count.
13. **Element Outline Thickness and Element Shrink "%" values** (0.02 and 0.06): fraction or percent? Assumed to
    be a fraction.
14. **Scale factor duplication:** Shader Options "Vector/Displacement Scale Factor", Window Options "Scale
    Factor", Load Frame Data "Scale Factor" and the `<Scale>` command argument may be one variable or several.
    The proposal is a single variable per active plot.
15. **Bubble size mapping:** linear in value, in |value|, or in normalised colour position? Is there a minimum
    size?
16. **Frame file exact layout:** whether the node id is a column, and how ncols counts it; how complex values are
    stored (Re/Im column pairs or amplitude/phase); how the legacy header differs; the exact list-file layout
    after the skipped first line; and how `*.impani` (impedance?) frames are meant to be shown. The manual says
    only that the list format is "the same format frame listing file format from previous PREP".
17. **SASSIani.xml schema and location:** they are only "system dependent default location". Propose a per-user
    application-data directory (e.g. `platformdirs.user_data_dir("ACS_SASSI_EDU")`) and an `<Animation>` element
    with attributes description, directory, type, frames, start, end, stride, scale, cmin, cmax, listfile.
18. **Which column the Vector plot uses** for X, Y and Z, and how the Deformed plot takes u_x, u_y, u_z from a
    frame with more columns (e.g. ACC frames that may include rotations in `\ACCR`). Deformed and Vector have no
    `Col` argument.
19. **The Vector plot's complex rule for Y and Z** components, and for nodes that have all three. The manual gives
    only the X example. The generalisation is in Section 8.8.
20. **Model numbering:** tab captions show "Model 0", but the Cut dialog defaults to Model Number 1. Are
    reference numbers 0-based or 1-based?
21. **The main toolbar's "Change Analysis Options" points to Section 6.5.5 (Windows Settings).** It almost
    certainly opens Options ▸ Analysis (6.5.4). It has no command equivalent.
22. **"Open the Model Database"** has no listed command. Check whether it maps to the RESUME command, which loads
    binary model data saved by SAVE.
23. **NODESEL highlighting in other plots:** "Each plot type has its own different display rules", but only the
    Node plot rule (blue square) is described.
24. **Interaction-node definition** for the Node plot (red). It depends on the model's interaction-node data,
    specified in the model-construction spec.
25. **Soil Layer plot half-space row:** which properties fill it (the layer after the last finite one?), and how
    TOPL affects what is drawn (only top layers, or embedment layers?).
26. **Time History "Plot Type" in Parse Frame Data** (anitype 3): which animated plot it is meant for is not
    stated. Treat it as a database tag only.
27. **Export Table:** "certain 2D plots" are not listed. Propose Spectrum, Time History and Soil Properties.
