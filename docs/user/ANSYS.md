# Working with ANSYS models in SASSI-EDU

This guide is for structural engineers who build their models in ANSYS Mechanical APDL. It covers
three things:

* moving an ANSYS model into SASSI-EDU (`CONVERT,ANSYS`);
* writing a SASSI-EDU model back to ANSYS (`ANSYS`, `ANSYSREFORMAT`);
* checking that the two programs describe the same structure (damping conventions, beam axes,
  verification problem VP-A1).

All commands are described in requirements §3.4.K and decisions D-ANS-01 … D-ANS-08. The code
lives in `sassi/io/ansys_cdb.py` (import), `sassi/io/apdl.py` (export) and
`sassi/prep/commands/conversion.py` (commands).

To use the SSI *results* in ANSYS (the two-step approach of ACS SASSI Option A: SSI motions as
equivalent static or dynamic loads of a refined ANSYS model, module LOADGEN), see
[OPTION_A.md](OPTION_A.md); the APDL model that the `ANSYS` command writes is the structural model of
that second step.

> **Disclaimer (as in ACS SASSI).** The converters handle only a subset of ANSYS elements and
> options, and they have had limited testing. A conversion can finish without errors and still
> describe a different structure. Read every warning and check the converted model before you
> run an analysis (for example, compare fixed-base frequencies, see §6).

---

## 1. Quick start

### ANSYS → SASSI-EDU

In ANSYS (`/PREP7`, after meshing):

```
CDWRITE,DB,mybuilding,cdb          ! writes mybuilding.cdb (blocked format)
```

In SASSI-EDU:

```
* name and directory of the SASSI-EDU model
MDL,bldg,C:/work/bldg
* <model> blank = active model; gravity in model units
CONVERT,ANSYS,,mybuilding.cdb,9.81
* restrain rotations that no element stiffens
FIXROT
* then ETYPE / ETYPEGEN, INT / INTGEN for the SSI model
GROUNDELEV,0
* bldg.pre
WRITE
```

The full syntax is

```
CONVERT,ANSYS,<model>,<filename>,<gravity>,[<prefile>],[<damp>]
```

| Argument | Meaning |
|---|---|
| `<model>` | destination model number; blank = the active model. The model is **replaced** by the converted one. Its MDL name and path are kept, and the active model does not change. |
| `<filename>` | the `.cdb` file (CDWRITE). Plain APDL files with `N`, `EN`/`E`, `MP`, `R` … commands are read too (see §2.1). |
| `<gravity>` | **required.** Gravity in the model's length/time units: 9.81 (m, s), 32.2 (ft, s), 386.4 (in, s). It turns the ANSYS density into the SASSI specific weight (`weight = DENS·g`) and becomes the model `GRAVITY`. |
| `[<prefile>]` | SASSI-EDU extension: also write the converted model as a `.pre` file (the dialog's "Output .pre File Name"). |
| `[<damp>]` | SASSI-EDU extension: damping ratio for materials without `DMPR` (default 0, with a warning). |

The converter also writes an element map, `<model>_cdb.map` in the model directory (or
`<cdbname>_cdb.map` in the working directory when the model has no MDL). Each line holds
`ansys_elem  group  elem` (D-ANS-02), so you can find the SASSI element that came from any ANSYS element.

### SASSI-EDU → ANSYS

```
* optional: one BEAMS group per end-release pattern (from an empty active model)
ANSYSREFORMAT,1,beams.map
* writes <model>.inp in the model directory
ANSYS
```

Read the file into ANSYS with `/INPUT,bldg,inp`. The full syntax is
`ANSYS,[FileName],[Dir],[<dmap>]`, where `<dmap>` sets the damping mapping (§4).

---

## 2. What CONVERT,ANSYS reads

### 2.1 Records

| ANSYS record | Used for |
|---|---|
| `NBLOCK` (any Fortran format, e.g. `(3i9,6e21.13e3)`, `(3i8,6e16.9)`) | nodes; node numbers are kept. Fields are read by **column width**, as the format line says (adjacent `e21.13` fields can touch). Coordinates are global Cartesian; nodal rotation angles THXY, THYZ, THZX are reported, not applied. |
| `EBLOCK` (`SOLID` layout and the short non-solid layout; continuation lines) | elements with TYPE, MAT, REAL, SECNUM |
| `ET`, `KEYOPT` | element types and options |
| `MPTEMP`/`MPDATA` (and `MP`) | EX, PRXY or NUXY (or GXY), DENS, DMPR (or DMPS). The first temperature is used. Every other label is reported (§2.3). |
| `RLBLOCK` / `R`, `RMORE` | real constants |
| `SECTYPE`, `SECDATA`, `SECBLOCK`, `SECOFFSET`, `SECCONTROL` | beam, pipe, link and shell sections |
| `D` with value 0 | fixities (`ALL` = all six DOFs); non-zero values are reported. `D,<component>,…` on a node component (CMBLOCK, also when the CMBLOCK comes later in the file) fixes its nodes. |
| `CMBLOCK` node component `SSI_INT` | interaction nodes (the component that the export writes) |
| `N`, `E`/`EN`/`EMORE`, `TYPE`/`MAT`/`REAL`/`SECNUM` | plain APDL input, e.g. the file written by the `ANSYS` command |

Everything else is reported with a count: `CE`, `CP`, `CERIG`, `RBE3`, `ENDRELEASE`, loads (`F`, `SF`,
`BF`, `ACEL`), `TB` data, global damping (`ALPHAD`, `BETAD`, `DMPRAT`, `DMPSTR`), and unknown commands.

The reader does **not evaluate APDL**. A command that uses a parameter, an expression or a name where a
number is expected (`*SET,tid,4` then `ET,tid,170`; `N,2,L/2,0,0`) is skipped, counted and quoted in
the warnings, and reading goes on. Workbench input files (`ds.dat`) use parameters freely: write the
mesh from Mechanical APDL with `CDWRITE,DB` instead, which writes plain numbers.

Every option of a supported element that the converter neither uses nor knows to be harmless is
reported (spec 04 §11.1): non-zero KEYOPTs, real constants and `SECCONTROL` values that change the
ANSYS stiffness or mass. Output-only KEYOPTs (stress and force printout, layer data storage) are
accepted silently.

### 2.2 Elements

| ANSYS | SASSI-EDU | Notes |
|---|---|---|
| SOLID45, SOLID65, SOLID185 | SOLID (group type 1) | Nodes I…P → 1…8. Prisms, pyramids and tetrahedra keep their repeated nodes. SOLID185 needs KEYOPT(3) = 0 (D-ANS-07). Reduced integration and the mixed u-P formulation (KEYOPT(6)) are reported. SOLID65 cracking and rebar are not converted. |
| SOLID186, SOLID95, SOLID187 | SOLID (corner nodes) | Midside nodes are dropped (warning). A tetrahedron becomes I,J,K,K,L,L,L,L. |
| SHELL63 | SHELL (3) | Thickness = R1 (TK(I)). A tapered thickness is averaged (warning). ADMSUA (R18, added mass per area) is lumped as ADMSUA·A/n at the corners, like the lumped SASSI shell mass (warning). EFS (R5, elastic foundation), RMI (R7 ≠ 1) and KEYOPT(1) ≠ 0 (membrane-only or bending-only) are reported, not converted. |
| SHELL181, SHELL281 (corners) | SHELL (3) | Thickness = sum of the section layers. Multi-layer sections become one homogeneous layer (warning). The Mindlin behaviour of SHELL181/281 is not carried over: SASSI SHELL is a thin Kirchhoff facet. `SECCONTROL` added mass per area is lumped like ADMSUA (warning); KEYOPT(1) = 1 (membrane only) and the transverse shear stiffnesses E11, E22, E12 are reported. |
| BEAM4, BEAM44 | BEAMS (2) | Properties come from the real constants. BEAM44 releases come from KEYOPT(7)/(8), including releases at both ends (§2.4). Tapered BEAM44 sections are averaged (warning); offsets DX…DZ are not represented; BEAM44 real constants beyond R24 are reported. BEAM4 ADDMAS (R12) is lumped as ADDMAS·L/2 at both ends (warning); ISTRN and SPIN are reported. |
| BEAM188, BEAM189, PIPE288 | BEAMS (2) | Properties come from sections: RECT, ASEC, CSOLID, CTUBE and PIPE. BEAM189 I–J (midside K) becomes two elements I–K and K–J. ENDRELEASE is not possible. `SECCONTROL` ADDMAS is lumped as ADDMAS·L/2 (warning). Each element becomes one exact Timoshenko beam, so coarse ANSYS meshes with linear shape functions (KEYOPT(3) = 0) give different results (info). |
| LINK8, LINK180 | BEAMS, axial only | I2 = I3 = 0 and J = 1e-8·A². ROTX/ROTY/ROTZ are fixed only at nodes that no other element stiffens in rotation (nodes joined only by LINKs, SOLIDs, PLANEs, translational springs or masses); a LINK at a shell corner or a beam end leaves the rotations free. LINK180 ADDMAS is lumped (warning); tension-only (KEYOPT(3), TENSKEY) and ISTRN are reported. |
| COMBIN14 | SPRING (7) or GENERAL (9) | KEYOPT(2) = 1…6 gives an SC spring in that DOF. An axial spring (KEYOPT(3) = 0, 1, 2) along a global axis gives an SC spring. An oblique spring gives a 2-node GENERAL element with `K = k e eᵀ` (D-ANS-04). CV1/CV2 viscous damping is not converted. |
| MASS21 | `MT`, `MR` with `MUNITS 0` | Masses in mass units (D-ANS-05). Several MASS21 at one node are summed. KEYOPT(3) = 0, 2, 3, 4 are supported, and KEYOPT(1) = 1 (volume × DENS). |
| MATRIX27 | GENERAL (9) | Stiffness (KEYOPT(3) = 4) → MXR. Mass (KEYOPT(3) = 2) → MXM, in mass units. The real constants C1…C78 are the upper triangle, row by row. Damping matrices and unsymmetric matrices are not supported. |
| PLANE42, PLANE182 | PLANE (4) | 2-D models, see §2.5. Only plane strain is equivalent. |

The converter creates one group per ANSYS element type number (D-ANS-02) and numbers the elements
1…n in file order. Oblique COMBIN14 springs get an extra GENERAL group, numbered after the largest
type number. Converted elements have ETYPE 0. SOLID and PLANE elements below the ground elevation
would therefore be treated as excavated soil (the converter warns). Use `ETYPE`/`ETYPEGEN` to set
the role explicitly.

### 2.3 Materials

```
M,<mat>,EX,nu,DENS*g,beta,beta,1
```

* `nu` comes from PRXY, then NUXY, then `EX/(2 GXY) − 1`. Without any of them, the ANSYS default
  of 0.3 is used (warning).
* Orthotropic data (EY ≠ EX …) is not converted. The isotropic EX is used (warning). Orthotropic
  labels that only repeat the isotropic values (EY = EZ = EX, PRYZ = PRXZ = ν, GXY = GYZ = GXZ =
  EX/(2(1 + ν)) to 1e-6) are accepted silently; a GXY that does not match EX and ν is reported.
* Every other label (ALPX, KXX, C, …) is ignored with a warning (D-ANS-01).
* Damping: `beta = DMPR/2` for both βp and βs (D-ANS-01). `beta = DMPS/2` is used when only DMPS is
  given (both given: DMPR is used, warning). Otherwise the CONVERT `<damp>` argument applies
  (warning). DAMP/BETD/ALPD (Rayleigh) is ignored with a warning.

### 2.4 Beam axes, sections and releases (D-ANS-03)

| | ANSYS (BEAM4/44/188/189) | SASSI-EDU (BEAMS) |
|---|---|---|
| axis 1 / x | I → J | I → J |
| orientation node K | K lies in the element **x–z** plane; z points towards K | K lies in the **1–2** plane; axis 2 points towards K |
| no K node | y parallel to the global X–Y plane (y = global Y for a vertical member), z = x × y; BEAM4 THETA rotates y towards z | a K node is always required |

The converter computes the ANSYS frame of every element. It then uses or creates the SASSI K node so
that **SASSI axis 2 = ANSYS z** and **axis 3 = −ANSYS y**. One created K node is shared along a
straight member. The section properties then map as follows:

| SASSI `R` field | ANSYS RECT (B along y, H along z) | ANSYS ASEC / real constants |
|---|---|---|
| `axial` | B·H | A / AREA |
| `shear2` (along axis 2 = z) | 5A/6 | A/SHEARZ (BEAM4/44), SECCONTROL TXZ/G, else 0 |
| `shear3` | 5A/6 | A/SHEARY, SECCONTROL TXY/G, else 0 |
| `tors` | Roark J of the rectangle | J / IXX (IXX = 0 → IYY + IZZ) |
| `flex2` (about axis 2) | Izz = H·B³/12 | Izz / IZZ |
| `flex3` (about axis 3) | Iyy = B·H³/12 | Iyy / IYY |

An ASEC section with a product of inertia Iyz ≠ 0 is converted in its principal axes. Axis 2 lies
along the principal direction nearest to z, and `Iyz = ∫ y z dA` is assumed.

Releases (BEAM44 KEYOPT(7) at I, KEYOPT(8) at J): KEYOPT is a 6-digit number whose digits are, from
left to right, UX UY UZ ROTX ROTY ROTZ (1 = released). It maps to SASSI KI/KJ as
P1 P2 P3 M1 M2 M3 = UX UZ UY ROTX ROTZ ROTY. For example, KEYOPT(7) = 11 releases ROTY and ROTZ
(a moment hinge), which gives KI = 0,0,0,0,1,1.

SASSI does not allow the same component released at I and J (CHECK Error 10), but ANSYS does: a
pin-ended brace is KEYOPT(7) = KEYOPT(8) = 11. The converter maps such a member onto an equivalent
valid one and warns:

* **M2 or M3 released at both ends.** No end moment in that bending plane means no shear either, so
  the plane has no bending stiffness at all. SASSI represents that exactly with a zero inertia
  (M3 → I3 = 0, M2 → I2 = 0), and the then redundant releases are dropped. The pin-ended brace
  becomes an axial (+ torsion) member; unit test `test_pin_ended_brace_is_an_axial_member_closed_form`
  checks a braced cantilever against the closed form. The consistent mass of that plane follows the
  SASSI beam shape, and the torsional mass ρ(I2 + I3) loses the zeroed inertia.
* **P1, P2, P3 or M1 released at both ends.** One release already removes that force from the whole
  member (a second one only adds a mechanism, singular in ANSYS too), so the release is kept at I.

### 2.5 2-D models

ANSYS 2-D elements lie in the global X–Y plane, and SASSI PLANE elements lie in the X–Z plane. When
the file contains PLANE42 or PLANE182 elements, the whole model is rotated by +90° about X:
(X, Y, Z) → (X, −Z, Y). The DOF labels follow (UY ↔ UZ, ROTY ↔ ROTZ), and HOUSE `<dim>` is set to 1.
Counter-clockwise ANSYS elements stay counter-clockwise in the SASSI (x right, z up) view.

### 2.6 After the conversion

* Run `FIXROT`: SASSI nodes carry six DOFs, so solid-only nodes, shell drilling rotations and
  spring-only nodes need restraints (the converter prints a reminder).
* Set the SSI data: `GROUNDELEV`, `ETYPE`/`ETYPEGEN`, interaction nodes (`INT`, `INTGEN`, or the
  `SSI_INT` component), soil layers, and the analysis options.
* Midside nodes of quadratic elements and other unused nodes are left in the model. AFWRITE fixes
  them (Warning 4); `RMVUNUSED` (tier P1) removes them where your build has it.

---

## 3. What the ANSYS command writes

`ANSYS,[FileName],[Dir],[<dmap>]` writes `/PREP7` input:

* Default name: `<model>.inp` in the model directory.
* Without MDL, the file is `unnamed.inp` in the working directory, with a warning (D-MDL-03).
* `Dir` is a directory, relative to the working directory.

| SASSI-EDU | legacy (default, ANSYS V11–15) | modern (`EDUOPT,ANSYSMODERN,1`) |
|---|---|---|
| SOLID (structure) | SOLID45; KEYOPT(1) = 0 when MOPT `<incomp>` = 0 (extra shapes), else 1 | SOLID185; KEYOPT(2) = 3 with incompatible modes, else 0 |
| PLANE (structure) | PLANE42, KEYOPT(3) = 2 (plane strain) | PLANE182, KEYOPT(3) = 2 |
| SHELL | SHELL63, `R` = THICK | SHELL181 + `SECTYPE,SHELL` |
| TSHELL | SHELL181 | SHELL181 |
| BEAMS | BEAM4 (no releases); BEAM44 with KEYOPT(7)/(8) from KI/KJ | BEAM188 + `SECTYPE,BEAM,ASEC` (+ `SECCONTROL` for the shear areas); releases are **not** written (warning) |
| BEAMS with I2 = I3 = 0 | LINK8 | LINK180 |
| SPRING | COMBIN14, KEYOPT(2) = 1…6, one element per non-zero SC component | same |
| GENERAL | MATRIX27 KEYOPT(3) = 4 (K_R) and KEYOPT(3) = 2 (mass) in global axes. 3-node local input is rotated with `TᵀKT` | same |
| MT / MR | MASS21 KEYOPT(3) = 0, mass units (weight / g when MUNITS = 1) | same |
| M | `MP,EX`, `MP,PRXY`, `MP,DENS` = weight/g, `MP,DMPR` (§4) | same |
| D | `D,n,ALL,0` for fully fixed nodes, otherwise one label per fixed DOF the node has in ANSYS | same |
| interaction nodes | `CMBLOCK,SSI_INT,NODE` | same |
| groups | `CMBLOCK,SSI_G<group>,ELEM` | same |
| excavated soil (ETYPE 2 or resolved 2), soil layers | not exported | not exported |

Beams keep their SASSI K node, so ANSYS z = SASSI axis 2 and IZZ = I2, IYY = I3, SHEARZ = A/As2,
SHEARY = A/As3. TKZ/TKY (stress output only) are the dimensions of the equivalent rectangle.

ANSYS sets beam releases per element **type** (KEYOPT), while SASSI sets them per element. The
exporter writes one ANSYS type per (group, release pattern) and warns when a group has several
patterns. `ANSYSREFORMAT,<Org>,<Map>` instead copies model `<Org>` into the empty active model with
one BEAMS group per pattern. The first pattern keeps the group number; EOUT requests follow the
elements. It also writes `old_group old_elem new_group new_elem` to `<Map>`.

2-D models (PLANE groups) are rotated back to the ANSYS X–Y plane.

Expect small frequency differences that come from the element formulations, not from the converter:

* SASSI SOLID uses ½ lumped + ½ consistent mass, and SHELL uses lumped translational mass. ANSYS
  uses consistent mass unless you set `LUMPM,ON`.
* SHELL181 and BEAM188 are shear-deformable.
* SOLID185 with KEYOPT(2) = 0 uses B-bar integration.

---

## 4. Damping conventions (R1 §1.2, D-ANS-01, D-ANS-06)

SASSI-EDU uses the SHAKE/SASSI complex modulus. Each modulus is multiplied by

```
c(β) = 1 − 2β² + 2iβ√(1 − β²)        (|c| = 1, loss angle δ = 2 arcsin β)
```

ANSYS full-harmonic analysis with structural damping uses `K(1 + i g)`. Two mappings are available:

* **`<dmap>` = 0 (default, D-ANS-06):** `MP,DMPR,mat,2β`. This is the `1 + 2iβ` form, and E is
  unchanged. CONVERT,ANSYS reads it back as β = DMPR/2.
* **`<dmap>` = 1 (exact):** `E_ANSYS = E(1 − 2β²)` and `g = 2β√(1 − β²)/(1 − 2β²)`. This gives
  `E_ANSYS(1 + i g) = E·c(β)` exactly. A model with `CMODFORM,1` uses `c(β) = 1 + 2iβ` instead
  (D-CNV-03); its exact mapping is `E_ANSYS = E`, `g = 2β`, and that is what `<dmap>` = 1 writes for
  it (the file header says which form was used).
* **`<dmap>` = 2:** no damping data.

| β | DMPR (`<dmap>` 0) = 2β | E factor (`<dmap>` 1) = 1 − 2β² | g (`<dmap>` 1) | error of 1 + 2iβ: magnitude | error of 1 + 2iβ: loss angle |
|---|---|---|---|---|---|
| 1 % | 0.0200 | 0.999800 | 0.020003 | +0.02 % | −0.02 % |
| 2 % | 0.0400 | 0.999200 | 0.040024 | +0.08 % | −0.06 % |
| 3 % | 0.0600 | 0.998200 | 0.060081 | +0.18 % | −0.14 % |
| 5 % | 0.1000 | 0.995000 | 0.100377 | +0.50 % | −0.37 % |
| 7 % | 0.1400 | 0.990200 | 0.141039 | +0.98 % | −0.73 % |
| 10 % | 0.2000 | 0.980000 | 0.203059 | +1.98 % | −1.47 % |
| 15 % | 0.3000 | 0.955000 | 0.310582 | +4.40 % | −3.21 % |
| 20 % | 0.4000 | 0.920000 | 0.425998 | +7.70 % | −5.52 % |

For the SSI benchmarks of requirements §6.6 (e.g. VP-20, an ANSYS full-harmonic soil box), use
`<dmap>` = 1.

SASSI allows different βp (on the constrained modulus) and βs (on G), but ANSYS has one damping value
per material. The export uses βs and warns when they differ. SC spring damping and the GENERAL
imaginary stiffness (MXI) have no frequency-independent COMBIN14/MATRIX27 counterpart, so they are
not exported (warning).

> **Check the meaning of DMPR in your ANSYS release.** The decisions above treat `MP,DMPR` as the
> constant structural damping coefficient g (= 2β), as the ACS SASSI manual's workflow implies.
> Some ANSYS releases document DMPR as a constant damping *ratio*: in a full-harmonic analysis
> that release adds (2·DMPR/Ω)K, i.e. g = 2·DMPR. Those releases also offer `MP,DMPS`, the
> structural damping coefficient g. If your release behaves that way, do one of the following:
>
> * replace `MP,DMPR,mat,2β` by `MP,DMPS,mat,2β` in the exported file, or by `MP,DMPR,mat,β`;
> * on import, give the materials `DMPS` instead of `DMPR` (CONVERT reads β = DMPS/2).
>
> To check, run a one-element harmonic test in ANSYS and compare the resonant amplification with
> 1/(2β).

---

## 5. Round trip helper: CONVERT,SSI (D-ANS-08)

`CONVERT,SSI,<model>,<file.hou>,[<prefile>]` rebuilds a model from the HOUSE deck written by AFWRITE.
It also reads `<file>.sit` and `<file>.poi` when they exist. The result is a model whose AFWRITE gives
the same decks again (unit test `test_convert_ssi_round_trip`). Things the decks do not carry:

* local coordinate systems (nodes come back in global coordinates);
* loads and analysis options other than SITE, POINT, HOUSE, HOUSEX, MOPT, CMODFORM, INCOH, WPASS, ME,
  SYMM, AMP, WAVE, TOPL and FREQ;
* the AFWRITE gap nodes and the all-fixed flags of nodes without DOFs, which are removed again.

ETYPE comes back resolved for SOLID/PLANE. Legacy SASSI2000 fixed-format decks are not supported.

---

## 6. Verification

* **VP-A1** (`sassi/verify/problems/vp_ansys.py`, `pytest tests/verification/test_vp_ansys.py`, or
  `VERIFY,VP-A1`). A 20 m steel cantilever, BEAM188 RECT 0.2 (y) × 0.4 (z), is converted from
  `tests/data/ansys/beam188_cantilever.cdb` (K node) and from `beam188_default_orientation.cdb`
  (no K node). Checks:
  * the first 8 fixed-base frequencies equal those of the native SASSI-EDU model to 1e-8;
  * Euler–Bernoulli frequencies in both bending planes to 1 % (the difference is Timoshenko shear,
    ≤ 0.4 %);
  * the beam-axis mapping, for both converted files and both re-imported models: the bending
    frequencies sorted by mode-shape direction (deflection along Y or along Z) equal the native ones
    (1e-8), and the tip deflections under unit Y and Z loads equal `PL³/(3EI) + PL/(G As)` with Izz
    and Iyy respectively (1e-8). A sorted frequency list alone cannot verify the mapping: swapping
    I2 and I3 of a straight prismatic cantilever leaves its spectrum unchanged;
  * export → import round trips (BEAM4 and BEAM188/ASEC) keep the frequencies (1e-8).
* Unit tests: `tests/unit/test_ansys_cdb.py` (reader, every supported element type, unsupported
  items, materials, beam axes, BEAM44 releases checked against 3EI/L³ and 12EI/L³),
  `tests/unit/test_ansys_conversion_rules.py` (LINK nodes at shells and other elements, pin-ended
  BEAM44, reported options, added masses, APDL parameters, CMODFORM),
  `tests/unit/test_ansys_apdl.py` (export text, damping mappings, model → APDL → model with
  identical assembled K and M, 2-D), and `tests/unit/test_ansys_commands.py`.
* Your own ANSYS cross-check (requirements §6.6): export a model with `ANSYS`, then run a modal
  analysis in ANSYS with the same mass assumptions:

  ```
  /INPUT,bldg,inp
  /SOLU
  ANTYPE,MODAL
  MODOPT,LANB,10
  LUMPM,OFF
  SOLVE
  ```

  Compare the frequencies with SASSI-EDU's fixed-base eigenvalues (`fixed_base_frequencies` in
  `sassi/verify/problems/vp_ansys.py`).

## 7. Limitations

The converter does not handle:

* coupling and constraint equations (CE/CP/CERIG/RBE3);
* beam ENDRELEASE and BEAM188 releases;
* section offsets and shell offsets;
* layered or orthotropic materials;
* nonlinear materials (TB);
* viscous spring damping;
* MATRIX27 damping or unsymmetric matrices;
* nodal coordinate systems (rotated nodes);
* loads;
* APDL parameters, expressions and `*DO` loops (the commands are reported and skipped);
* membrane-only or bending-only shells, elastic foundations (SHELL63 EFS), tension-only links,
  mixed u-P solids, initial strains (all reported);
* element types other than those in §2.2.

Option AA (`ANSYSMODELTYPE`, MATRIX50 super-elements, HOUSEFSA/ANALYSFA) is tier P2:
`ANSYSMODELTYPE` is stored only. `GENMATRIXDAMP` is "not usable in this version", as in the manual;
enter MXI rows directly.
