# Option NON: nonlinear structures by equivalent-linear SSI iterations

This guide is for structural engineers who know nonlinear analysis from ANSYS (time integration with
nonlinear materials) and want to understand how ACS SASSI, and SASSI-EDU, include **cracking
reinforced-concrete walls** and **nonlinear springs** (base isolators, pile-soil interfaces, sliding)
in a frequency-domain SSI analysis.

The method is that of the ACS SASSI manual (section 1.5.4, Fig. 1.2, section 6.5.4 "NONLINEAR Module
Options", section 9.17). The requirements are in `docs/spec/00_requirements.md` §3.4.O and §4.15 and in
decisions D-NON-01 … D-NON-12. The code is:

| Part | File |
|---|---|
| hysteresis models, equivalent-linear properties | `sassi/core/hysteresis.py` |
| panel geometry and kinematics, SHEAR, BBCGEN | `sassi/core/panels.py` |
| NONLINEAR module, COMB_XYZ_THD | `sassi/modules/nonlinear.py` |
| commands | `sassi/prep/commands/nonlinear_cmds.py` |
| verification problems VP-45, VP-46, VP-NON1 | `sassi/verify/problems/vp_nonlinear.py` |
| tutorial | `examples/ex07_option_non.pre` |

---

## 1. The idea in one page

The SSI solution of SASSI is linear: it is solved frequency by frequency with complex stiffnesses.
A nonlinear element is therefore replaced by an **equivalent-linear** element whose stiffness and
damping depend on how much it deforms, and the SSI analysis is repeated until the properties no longer
change. This is the SHAKE method applied to the structure.

For every nonlinear element and every iteration (spec 05d §3.5, requirements §4.15):

1. **Deformation history** `x(t)`: the shear strain of a wall panel, computed from the relative
   displacements of its four corner nodes, or the elongation `u_J - u_I` of a spring along its DOF.
2. **Equivalent amplitude** `x_eq = EDF · max|x(t)|`, with the *equivalent-linear displacement factor*
   EDF ≈ 0.8 (0.7 – 0.9; 1.0 gives too soft and 0.6 too stiff elements, manual §6.5.4).
3. **Hysteresis model**: `x(t)` is run through the force-deformation rules of the element (CMS, TAK or
   GMR, §4) to get the force history `F(t)` (files `.thd` / `.ths`: plot one against the other to see
   the loops).
4. **Secant stiffness**: `E_new = E_el · K_sec(x_eq)/K_el`, where `K_sec = F(x_eq)/x_eq` on the
   backbone curve and `K_el = Y1/X1` is the slope of its first segment (D-NON-11). For panels the
   modulus of the panel's material is scaled, with Poisson's ratio kept: shear, axial and bending
   stiffness degrade together (manual §1.5.4). For springs the spring constant of the DOF is scaled.
5. **Damping**: `xi = min(cutoff, scale·xi_h + [xi_el])` (D-NON-04), with the hysteretic damping
   `xi_h = E_D/(4 π E_S)` of the stabilised cyclic loop at `x_eq` (E_D = loop area,
   E_S = x_eq F(x_eq)/2), the *Damping Scale Factor* (0 means 1), the elastic damping `xi_el` added
   when *Include Elastic Damping* is set, and the *Damping Cutoff %* (0 = none; 7 % is the ASCE 4 level
   for cracked concrete).
6. **Ductility and force reduction**: `mu = max|x| / x_cr` (x_cr = BBC point 1, the end of the elastic
   range, not the yield point) and `F_mu = K_el·max|x| / |F(t*)|` at the time of `max|x|` (the
   inelastic absorption factor of ASCE 43-05), in `Panel.fmu` / `Spring.fmu`.
7. **Convergence** (D-NON-06): the iterations stop when the largest relative change of E (or k) is
   below 2 % and the largest change of the damping ratio below 0.5 % (absolute), at most 10 iterations.

The new properties go into a new HOUSE deck, `<model>_new.hou`; the next SSI analysis is a fast
**"New Structure" restart** of ANALYS (`<mode>` 1), which reuses the soil impedance of the first run:
only the structure changed.

---

## 2. The workflow (manual Fig. 1.2)

```
                  elastic SSI analysis (analysis 0)                       iterations k = 1, 2, ...
SITE (X, Y, Z input: FILE1X/Y/Z)                              FCOPY <model>_new.hou -> <model>.hou
POINT, HOUSE                                                  HOUSE
ANALYS  <save> 1 (restart files), <simul> 1 (X/Y/Z cases)     ANALYS  <mode> 1 = New Structure
for X, Y, Z: MOTION (Save Complex TF) + RELDISP               for X, Y, Z: MOTION + RELDISP
             NONLINTHD,<dir>  (keep X_*, Y_*, Z_* THD files)               NONLINTHD,<dir>
COMBXYZTHD  (COMB_XYZ_THD: X + Y + Z per time step)           COMBXYZTHD
RUNNONLINEAR (no PANEL.NON: elastic run, PANEL.NON = 1)       RUNNONLINEAR (PANEL.NON = 1: iteration)
NONLINSAVE,elastic                                            ... until converged (NONLINITER)
```

Why three directions are combined (manual §1.5.4 WARNING): the nonlinear behaviour must be driven by
the **simultaneous** three-component input, so the relative displacements of the X, Y and Z runs are
added time step by time step (COMB_XYZ_THD) before NONLINEAR reads them. The same record is used per
direction as in example 5; real analyses use three statistically independent records.

The NONLINEAR module reads (D-NON-01):

* `<model>.eql` — its input deck, written from the EQL, P, S and BBC data of the model by AFWRITE (AOPT
  NONLINEAR flag) and by RUNNONLINEAR (keyword deck like the other decks, see §10);
* `<model>.hou` — the HOUSE deck of the analysis just made (it holds the elastic or iterated
  properties; NONLINEAR warns when they differ from those it expects);
* the `.THD` files of RELDISP (combined) of the panel corner nodes (`nnnnnTR_X/Y/Z.THD`) and spring
  end nodes;
* the state files `PANEL.NON` / `SPRING.NON` and `Panel_EQL_Matl_Prop.txt` / `SPRING_EQL_Matl_Prop.txt`.

**State machine** (manual §6.5.4): no `.NON` file → *elastic run* (the results are those of the
elastic model; the `.NON` file is created with `1`); `.NON` = 1 → *iteration*: the properties of the
analysis just made are read from `*_EQL_Matl_Prop.txt` (they are compared with the new ones for the
convergence test). To restart from the elastic analysis delete both files (`NONLINRESET`) and rewrite
the elastic HOUSE deck (AFWRITE).

Outputs of one run:

| File | Content |
|---|---|
| `Panel0001.thd`, `SPRING0001.thd` | deformation history x(t) (shear strain, curvature or relative displacement) |
| `Panel0001_AXIAL.thd` | uniform vertical axial strain of the panel (no hysteresis model: informative) |
| `Panel0001.ths` | nonlinear force history F(t) from the hysteresis model |
| `Panel0001.crv` | equivalent-linear curves versus amplitude: E/E_el, xi_h, K_sec, F (written by the elastic run) |
| `Panel_EQL_Matl_Prop.txt` | E_el, E_new, ratio, xi_el, xi_h, xi_new, max\|x\|, x_eq per panel, with the iteration number |
| `Panel.fmu` | max\|x\|, x_cr, ductility mu, elastic and nonlinear forces, F_mu |
| `<model>_new.hou` | the HOUSE deck of the next iteration |
| `NONLINEAR_CONVERGENCE.TXT` | one row per run: analysis, max \|dE/E\| %, max \|d xi\| %, converged |
| `<model>_NONLINEAR.out` | the listing: options, BBCs, panels and springs, results, convergence |

`SPRING` replaces `Panel` for springs (manual §6.5.4).

---

## 3. Defining the nonlinear elements

### 3.1 Global options: EQL

```
EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>
* EDF 0.8, panels, damping cut-off 7 %, scale 1, elastic damping added
EQL,0.8,1,7,0,1
```

Example 7 uses `EQL,0.8,1,0,0,1` (no cut-off); §11 shows what the 7 % cut-off changes there.

`<NonLinOpts>` is a bit mask (D-NON-02): 1 panels, 2 springs, 4 beams (not available). Blank = the
element types for which P / S records exist. `<dampCutoff>` is in **percent**.

### 3.2 Backbone curves (BBC)

A BBC is given **without the origin**; point 1 is the **cracking point** (the end of the elastic
range) and `<yield>` is the number of the yield point. Units are never converted: a shear BBC is
*shear strain (decimal) – shear force* and its first slope must equal `G·A_shear` of the panel (`t·L`);
a spring BBC is *displacement – force* with first slope `k`. RUNNONLINEAR warns when the first slope
differs by more than 1 % (EDU-09).

| Command | Meaning |
|---|---|
| `BBC,<num>,<type>,<points>,<yield>,<file>` | read `<points>` X Y pairs from a file |
| `BBCX,<num>,<points>,<yield>,<X1>,…` / `BBCY,…` | the X (strain) / Y (force) vectors |
| `BBCP,<num>,<point>,<X>,<Y>` | one point (a row of the BBC grid); `<point>` = npoints + 1 appends; skipped rows are (0, 0) until defined (warning) |
| `BBCI,<num>,<yield>,<type>` | yield point number and type (1 CMS, 2 CMB, 3 TAK, 4 GMR; titles only) |
| `DELBBC,<start>,[<end>],[<stride>]` | delete curves |

Beyond the last point the force is held constant (no failure branch); the listing flags elements that
pass the last point. The BBC must not soften (Y non-decreasing); CMS and TAK need strictly increasing
forces (GMR accepts a flat plateau, e.g. an elastic-perfectly-plastic spring).

The manual's warning applies: **use smooth BBCs** — sharp corners (bilinear or trilinear curves)
hurt the convergence of the iterations (see §8).

### 3.3 Wall panels: P, PNLGEN, PLIST

```
* a SHELL group in a vertical plane
P,<num>,<group>,<bbc>,<disp>,<force>
* one panel per vertical shell group (BBC = panel number)
PNLGEN
* group, corners, L, H, t, material, BBC and the checks
PLIST
```

A panel is a group of coplanar shells in a vertical plane, ideally rectangular; it may have any
orientation in plan. Only its four corner nodes are used. They are found by **node-connection
counting** (nodes used by exactly one element of the group); when that does not give four nodes
(triangles, openings) the nodes nearest to the corners of the panel's bounding rectangle are used,
with a warning. The local axes are `e_h` (horizontal, in the panel plane, oriented toward +X, or +Y for
a wall normal to X) and `e_v` = Z.

Panel deformations (requirements §4.15; rigid translations and in-plane rotations give zero):

```
gamma = 1/2 [(u_TL - u_BL) + (u_TR - u_BR)]/H + 1/2 [(w_BR - w_BL) + (w_TR - w_TL)]/L     (Disp. Opt 1)
kappa = [(w_TR - w_TL) - (w_BR - w_BL)] / (L H)                                         (Disp. Opt 2)
eps_v = 1/2 [(w_TL - w_BL) + (w_TR - w_BR)]/H                                           (axial, informative)
```

with `u = d·e_h`, `w = d·e_v`, L the width and H the height of the corner rectangle.

Each panel needs a material of its own (D-NON-12), because NONLINEAR changes E panel by panel:
GROUPMAT (one material per group) or explicit M / MACT. In this version (manual Errors 126 and 128)
panels use Disp. Opt 1 and Force Opt 1 (CMS); `EDUOPT,NONEXT,1` unlocks Takeda (3) and bending
panels (2) as experimental (D-NON-05). CMB (2) is not included (manual).

### 3.4 Nonlinear springs: S

```
* an existing SPRING element, DOF 1/2/3, Force Opt 4 (GMR)
S,<num>,<group>,<elem>,<bbc>,<disp>,<force>
```

A nonlinear spring acts on one translational DOF; split 3-DOF springs into 1D springs (manual
§1.5.4), and give each nonlinear spring an SC property of its own (it is updated). The relative
displacement is `u_J - u_I` of the element's nodes (both must be in the RELDISP requests and must not
be fixed).

### 3.5 Output requests: NONLINMOTDISP

`NONLINMOTDISP` adds the panel corners (X, Y, Z) and the spring end nodes (spring DOF) to the MOTION
requests (NOUT with the transfer-function flag) and the RELDISP requests (RDND), without duplicates.
Run it **before AFWRITE**. MOTION must save the complex TFs (`<cplx>` = 1) and RELDISP uses the
free-field reference (RELFILE blank): the reference cancels in every difference NONLINEAR forms.

---

## 4. The hysteresis models and their sources

| Code | Model | Use | Source of the rules implemented |
|---|---|---|---|
| 1 | Cheng-Mertz Shear (CMS) | panels in shear (low-rise walls) | subroutine HYST04 ("S1 shear hysteresis model") of INRESB-3D-SUP: Cheng, F.Y. and Mertz, G.E. (1989), *A computer program for inelastic analysis of 3-dimensional reinforced-concrete and steel seismic buildings*, Civil Engineering Study 89-31, University of Missouri-Rolla (NSF, NTIS PB90-123225), Appendix C pp. 151-153; model described in Cheng & Mertz (1989) CE Study 89-30 and Cheng, Mertz, Sheu & Ger (1993), J. Struct. Eng. 119(11) |
| 2 | Cheng-Mertz Bending (CMB) | — | not included in this version (manual) |
| 3 | Takeda (TAK) | bending/shear (experimental) | Takeda, Sozen & Nielsen (1970), J. Struct. Div. ASCE 96(ST12):2557-2573, with Otani's (1974, UILU-ENG-74-2029) rules as in HYST06 of INRESB-3D-SUP |
| 4 | General Masing Rule (GMR) | springs | Masing (1926) branches with factor 2; extended memory rules (Pyke 1979; Kramer 1996, §6.4.2) |

**GMR.** Virgin loading follows the BBC. A branch starting at a reversal point (x_r, F_r) is
`F = F_r + 2·F_bb((x − x_r)/2)`. A branch that reaches the previous reversal point continues on the
branch it left (the small loop closes), and a branch that reaches the largest past excursion joins the
backbone. Symmetric cycles give closed loops with `E_D = 8∫₀ˣF dx − 4xF(x)`; for an elastic-perfectly
plastic BBC `xi_h = 2(x − x_y)/(πx)` and `K_sec = F_y/x` (VP-45).

**CMS** (rules numbered as in HYST04 and the manual's Fig. 1.3; PC, DC = cracking point, SI = PC/DC,
DMAX/PMAX = largest past displacement/load in either direction, PM/DM those of the current direction):

* before the first cracking: elastic, `K = SI`; rule 1: loading on the multi-linear backbone;
* rules 2-4, unloading in three force bands measured from PM: above `PA = PM − PC` with
  `S1 = SI·min(1.4675 (DC/DMAX)^0.345, 1)`, down to `PB = PC/2` with `S2 = SI·min(0.7761 (DC/DMAX)^0.5195, 1)`,
  to zero with `S3 = SI·min(0.0707 (DC/DMAX)^1.369, 1)`; the two lower slopes are steepened when needed
  so that the unloading reaches zero no later than `DO`, where the line through the two peaks
  crosses zero; rule 5 unloads toward the reversal point of an open small loop;
* rules 6-7, reloading toward `(D2, P2)`: the point at 0.95·PMAX on the first unloading branch from the
  largest peak (mirrored to the current side), or the previous peak of a direction that is not cracked;
* rules 8-9, **pinching**: after a load reversal, below 0.75·PC, reloading with the slip stiffness
  `SR = SI·min((DC/DMAX)^1.02, 1)` to PC/4, then with the harmonic mean of the slip and the reloading
  stiffness to 0.75·PC — unless the direct line to (D2, P2) is steeper (a high cracking point gives
  unpinched loops, Fig. 1.3);
* rule 10, from (D2, P2) toward `(1.04·DMAX, PMAX)` (degradation factor 1.04 when the current
  direction holds DMAX) until the backbone is met; rule 11, reloading inside small loops toward their
  stored reversal points (up to 10).

HYST04 is load-incremental (it returns the tangent stiffness and the load at which the rule changes);
SASSI-EDU follows the same straight branches exactly in displacement control, changing the rule at its
load limit or at a reversal. Two tests of HYST04 that the original resolves by the small overshoot of
its load steps are made deterministic: a direction counts as cracked once the motion reaches its
cracking displacement, and the rule-10 target on the backbone (`|PX| > |PMAX|`) is tested with a
relative tolerance 1e-9. Beyond the last BBC point HYST04 lets the wall fail (zero stiffness); here the
force is held (a practically flat extension) so that the linearisation can continue, with a listing
note.

**TAK**: trilinear backbone through the cracking point (BBC point 1) and the yield point (BBC yield
point), flat beyond: the manual — "the yield point force is used in the Takeda model for the wall panel
peak capacity". Unloading before yield aims at the cracking point of the opposite direction; after yield
`K_u = K_cy (D_y/D_max)^0.4` with `K_cy = (P_y + P_c)/(D_y + D_c)` (D-NON-03); after a zero crossing the
reloading aims at the peak of the opposite direction (or at its yield point if that line is steeper and
the direction has not yielded); an interrupted reloading leaves its point as the first target of the
next reloading of that side (Otani's points U0, U1 …, kept as a stack). The manual allows TAK only with
the constant elastic damping: SASSI-EDU uses `xi = xi_el` for TAK elements (listing note).

The `.crv` curves use stabilised loops: three symmetric cycles at each amplitude, the last one measured;
E/E_el is the backbone secant (the BBC, or the trilinear Takeda curve).

---

## 5. Shear capacities and backbone generation: SHEAR, BBCGEN

```
* 0 = all panels
SHEAR,<panel>,[fc],[fy],[P],[Nu],[Abe],[fybe]
BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]
```

The panel geometry is taken from the model: `h_W` = vertical extent, `l_W` = horizontal extent,
`t_W` = shell thickness, `A_W = l_W t_W`. Units follow the gravity constant (British: ft, ksi, kips if
g > 20; SI: m, kN/m², kN). The equations use √f'c in psi internally (spec 10 §3.11):

```
ACI 318-08:            V = (alpha_c sqrt(f'c) + rho_H f_y) A_W  <=  10 sqrt(f'c) A_W
Wood 1990:             6 sqrt(f'c) A_W  <=  V = rho_V A_W f_y / 4  <=  10 sqrt(f'c) A_W
Barda et al. 1977:     V = (8 sqrt(f'c) - 2.5 sqrt(f'c) h_W/l_W + N_U/(4 l_W t_W) + rho_V f_y) t_W (0.6 l_W)
Gulec-Whittaker 2009:  V = (1.5 sqrt(f'c) A_W + 0.25 F_VW + 0.20 F_BE + 0.40 N_U) / sqrt(h_W/l_W)
```

`alpha_c` = 3.0 for h_W/l_W ≤ 1.5, 2.0 for ≥ 2.0, linear in between; `F_VW = rho_V A_W f_y`,
`F_BE = A_BE f_y,BE` (arguments 6-7 are A_BE and f_y,BE, D-NON-10; `EDUOPT,SHEARFORCEARGS,1` reads
them as forces). The manual recommends Wood or Gulec-Whittaker for typical nuclear walls and warns
that Barda (barbell walls) and ACI can overestimate the capacity.

BBCGEN (D-NON-09) writes a **22-point** curve of type 1 (CMS) to the panel's BBC number (the panel
number when 0): point 1 = cracking `(V_cr/(G A_W), V_cr)` with `V_cr = 3√f'c A_W` (ASCE 4-17
C.3.3.2) or `CFL·V_u` (CrackingForceLevel 0.10-0.50); points 2-21 equally spaced in strain up to the
yield point `(0.004, V_u)` on `V = V_cr + (V_u − V_cr)[1 − (1 − ξ)²]`; point 22 = failure
`(0.02, 1.02 V_u)`. Its first slope is `G·A_W` by construction.

---

## 6. Building a panel model

Large models are prepared with these commands (spec 11 §2.5; run them on a copy of the model):

| Command | Action |
|---|---|
| `WALLFLR` | delete non-shell elements; each plane with ≥ 5 shells becomes a group (titles: WALL along X/Y, FLOOR, OBLIQUE); the rest go to group 1 |
| `PANELIZE` | split each shell group along its intersections with other groups |
| `EDGE,<group>,[X],[Y],[Z]` | split a wall along the lines of its outer and opening edges (flag 1 = ignore edges parallel to that axis); cells top first, then left to right (Fig. 1.5: `EDGE,1,0,0,1` gives three strips, `EDGE,2` the piers 2, 4, 5) |
| `EDGEMODEL,[x],[y],[z]` | EDGE on every wall group |
| `UNIPNL,<group>` | one group per element (curved walls: one panel per shell) |
| `GROUPMAT` | one material per group (required: D-NON-12) |
| `DGRDFLR,<scale>` | scale E of the floor materials |
| `MERGEPANEL,<model>` | append the panel model's groups and materials to the original model (delete its shell groups first) |
| `PNLGEN` | one panel record per vertical shell group |
| `SOLIDPILE,<group>,[stiff],[soft],[stiff2]` | separate a SOLID pile group by interface springs (X, Y: 1e7, Z side: 10, Z tip: 1e7; 4 % damping; one SC property per spring, targets of `S`) |

---

## 7. Running the iterations

`RUNNONLINEAR` writes `<model>.eql` and runs the module. The loop is written with the SASSI-EDU
driver `NONLINITER` (D-NON-06):

```
VAR,NLHOUSE,"FCOPY,ex07_new.hou,ex07.hou",RUNHOUSE,RUNANALYS
VAR,NLDX,"AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0","EDUOPT,TFFILE,FILE8X","MOTION,...",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,X"
...
VAR,NLNON,"COMBXYZTHD,COMB_XYZ_THD.inp",RUNNONLINEAR
NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10
```

NONLINITER runs the commands of the variables once per pass (`#` = pass number) until NONLINEAR reports
convergence or 10 passes. AFWRITE inside the loop must write only the MOTION and RELDISP decks
(`AOPT`), otherwise it would overwrite the iterated HOUSE deck. The ANALYS deck of the iterations
(`<mode>` 1, New Structure) is written once, before NONLINITER. FCOPY warns in every pass that it
replaces `ex07.hou`: that is the intended step of Fig. 1.2.

`NONLINBAT,1` writes this script for the active model (`<model>_NONLINBAT.pre`) and a
`COMB_XYZ_THD.inp` that lists the THD files (`<file> X_<file> Y_<file> Z_<file>`); `NONLINBAT,0` writes
the single-direction version. The script uses the model's MOTION options for every direction: add a
`MOTION` command to a direction's variable to scale it differently (example 7 scales Z to 0.4 g; for this
symmetric building the vertical input leaves the panel shear strains practically unchanged, and the generic script
reproduces the example's convergence history). `NONLINSAVE,<tag>` keeps copies of the results of a pass under the names
of Fig. 1.2 (`Panel0001_elastic.thd`, `Panel_EQL_Matl_Prop_It3.txt`, `ex07_It3.hou` …);
`NONLINRESET` deletes the state files.

The batch form of the module, as the original executables:

```
python -m sassi.modules.nonlinear < NONLINEAR.inp            # lines: model, model.eql, listing
python -m sassi.modules.nonlinear COMB_XYZ_THD COMB_XYZ_THD.inp
```

---

## 8. Convergence: what to expect

Equivalent-linear iterations converge well when the deformation of the nonlinear element is
*displacement controlled* (a soft element in a long-period system, e.g. base isolators) and slowly when
it is *force controlled* (a stiff wall whose shear is set by the inertia forces): then the fixed point is
where the backbone force equals the demand, and each iteration reduces the distance to it only by the
factor `1 − (x F'/F)` — close to 1 just after cracking, where the backbone is nearly flat. If the
softening also moves the structure toward the peak of the input spectrum, the iterations can drift away
(this is the manual's warning about sharp BBC corners: "oscillation or divergence"). Practical advice:

* use smooth BBCs and check them with PLIST and the `.crv` curves;
* use enough SSI frequencies around the structural frequencies, which fall as the walls crack (manual
  §1.5.4);
* look at the convergence history (`NONLINEAR_CONVERGENCE.TXT`); a change that grows from one
  iteration to the next means divergence, not slow convergence;
* the 2 % / 0.5 % test compares two successive iterations, so a slowly converging sequence passes it
  before reaching the fixed point; compare the last properties with the `.crv` curve at the last x_eq.

---

## 9. Verification

| VP | What is checked | Result |
|---|---|---|
| VP-45 | elastic-perfectly-plastic BBC under GMR: closed loops, `xi_h = 2(x − x_y)/(πx)`, `K_sec = F_y/x` at 1.25 … 50 x_y (1e-6); the module on a harmonic history; the state machine (elastic run creates SPRING.NON; .NON = 1 iterates; the properties of SPRING_EQL_Matl_Prop.txt are used; deleting both restarts) | pass: closed-form errors ≤ 2e-16; the state-machine check 5.9e-8 (the convergence file prints 5 decimals) |
| VP-46 | SHEAR equations of the spec 10 §3.11 example (British, 1e-6) and BBCGEN (22 points, yield 21, failure point) | pass: British and BBCGEN to 1e-6; the three SI references are printed to 0.1 kN and are compared within half their last digit, 0.05 kN (lead decision D-W3-08; observed ≤ 0.031 kN), and the SI computation path agrees with the British results × 4.44822162 kN/kip to 1e-6 |
| VP-NON1 | SDOF with a nonlinear spring on a rigid base, harmonic input, iterated through HOUSE → ANALYS restart → MOTION → RELDISP → NONLINEAR: every response equals the closed-form steady state, the iteration follows an independent closed-form fixed-point iteration (1e-6), the converged stiffness and damping equal the BBC secant and the Masing loop damping at the converged amplitude, convergence in 5 iterations | pass: errors ≤ 1e-8; converged k/k_el = 0.7708, ξ = 0.1611 |

Unit tests: `tests/unit/test_hysteresis.py`, `test_nonlinear_panels.py`, `test_nonlinear_cmds.py`,
`test_nonlinear_module.py`, `test_nonlinear_example.py`.

---

## 10. Interpretations and limitations of this build

* **The `.eql` deck** is a keyword deck (D-NON-01). NONLINEAR is an ordinary module (lead decision
  D-W3-09): AFWRITE writes the `.eql` when the NONLINEAR flag of AOPT is set, and RUNNONLINEAR writes it
  again before it runs the module. CHECK applies the manual's Errors 121-128 (Error 121 when the model
  has neither panels nor springs; Errors 126 and 128 are lifted by `EDUOPT,NONEXT,1`) and reports the
  further findings of the deck builder as EDU-44 (error) and EDU-45 (warning).
* The CMS and TAK rules come from the INRESB-3D-SUP listing and Otani's description; ACS SASSI's own
  implementation is not published, so loops may differ in details (pinching levels, inner loops).
* No axial hysteresis model (manual): the axial strain is written for information only; following
  ASCE 4-17 and the manual, treat the structure as nonlinear under the horizontal components.
* Nonlinear beams (B, NonLinOpts 4) are not usable (manual).
* The relative displacements of fixed nodes are not available from RELDISP: nonlinear springs must
  connect free nodes.
* Gravity (axial) forces in the walls must be included by the analyst in the BBC (SHEAR/BBCGEN `Nu`).

For ANSYS users: the converged `E` of each panel and `k` of each spring are the cracked/secant
properties to use in an ANSYS model of the same building (e.g. a harmonic or modal analysis with the
cracked stiffness); the damping ratios map to ANSYS material damping as for the elastic model.

---

## 11. Tutorial: example 7

`examples/ex07_option_non.pre` is a complete Option NON analysis of a two-storey reinforced-concrete
shear-wall box (40 ft x 40 ft, storeys of 13 ft, 1 ft walls, 875 kips of equipment on each slab) on a
surface mat, 32.5 ft of soil (Vs 1,300 ft/s) on rock, under the RG 1.60 record of the examples scaled to
0.6 g in X and Y and 0.4 g in Z (Demo 9 of the manual is an RC shear-wall building at 0.60 g). Units are
ft, kip, s (`HOUSE,32.2`). It runs in 30–120 s and leaves about 105 MB in `examples/ex07/`: 34 MB of
text results and 72 MB of binary inter-module files (FILE2, the ANALYS restart files `COOTKnnn`, …),
which can be deleted after the run.

| Step | Commands | Remarks |
|---|---|---|
| model | `L`, `TOPL` (25 sublayers of 1.3 ft), `N`/`FILL`/`NGEN`, `M,11` … `M,18`, `GROUP,2` … `GROUP,9` (SHELL, one per wall and storey), `GROUP,10` (slabs), `MT`, `INT`, `POINT`, `HOUSE` | one material per panel (D-NON-12) |
| panels | `PNLGEN`, `PLIST` | panels 1–4: storey-1 walls S, E, N, W; 5–8: storey 2 |
| capacities | `SHEAR,0,4,60,0.005,0` | British (g = 32.2: ft, ksi, kips): f'c 4 ksi, f_y 60 ksi, ρ 0.5 %; per wall ACI 2,821 kips, Wood 2,186 kips (its lower bound 6√f'c A_W), Barda 2,608 kips, Gulec-Whittaker 1,716 kips |
| backbones | `BBCGEN,0,1,4,60,0.005,0,0,0,0.3` | ACI 318-08 V_u, cracking at 0.3 V_u (V_cr = 846 kips, γ_cr = 8.82e-5) |
| options | `EQL,0.8,1,0,0,1`, `NONLINMOTDISP` | EDF 0.8; ξ = 4 % + ξ_h; the 12 corner nodes in NOUT/RDND |
| elastic SSI | three `SITE` runs (FILE1X/Y/Z), `ANALYS,0,0,0,1,1,…,1`, `RUNPOINT`, `RUNHOUSE`, `RUNANALYS` | restart files saved, X/Y/Z cases |
| elastic NONLINEAR | `NONLINBAT,1`, `NONLINRESET`, `VAR` NLDX/NLDY/NLDZ/NLNON, `FOREACH` …, `NONLINSAVE,elastic` | analysis 0 |
| iterations | `ANALYS,0,0,1,…` (New Structure), `AFWRITE`, `VAR,NLHOUSE,…`, `NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10` | analyses 1, 2, … |
| results | `MOTION` spectra of nodes 13, 63, 113 (mat, floor, roof centres) | converged model |

Convergence history (`NONLINEAR_CONVERGENCE.TXT`; analysis 0 = the elastic model):

| analysis | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| max \|dE/E\| % | 25.5 | 20.9 | 16.1 | 11.1 | 7.4 | 4.0 | 2.6 | **1.7** |
| max \|dξ\| % | 6.86 | 3.02 | 1.54 | 0.70 | 0.33 | 0.07 | 0.02 | **0.01** |

Converged properties (`Panel_EQL_Matl_Prop.txt`, `Panel.fmu`):

| Panels | E/E_el | ξ_h | ξ | max\|γ\| | γ_eq | μ | F_μ |
|---|---|---|---|---|---|---|---|
| 1–4 (storey 1) | 0.374 | 12.5 % | 16.5 % | 3.64e-4 | 2.91e-4 | 4.13 | 3.14 |
| 5–8 (storey 2) | 1.000 | 0 | 4.0 % | 7.3e-5 | 5.8e-5 | 0.82 | 1.00 |

The four storey-1 walls behave alike because the building is symmetric and the same record drives X and
Y. They are cracked: 4.1 times the cracking strain, but far from the yield strain of 0.4 %. Their shear
at the peak strain, 1,114 kips, is 39 % of V_u; the initial stiffness at the same strain would give
K_el·max\|γ\| = 3,493 kips, hence F_μ = 3.1. The storey-2
walls stay below the cracking strain (μ = 0.82), so they keep E_el and the elastic 4 %. The 16.5 %
damping is above the 7 % that ASCE 4 accepts for cracked concrete. The method gives it, but a design
analysis would cap it (`<dampCutoff>` 7).

The same model with one change (iteration counts of this build):

| Variant | Change | Result |
|---|---|---|
| the example | — | converged after 7 iterations |
| damping cut-off 7 % | `EQL,0.8,1,7,0,1` | not converged after 10 passes, although the changes still fall (25.5 … 3.9, 3.2 %). With less damping the walls deform more: E/E_el = 0.21, ξ = 7 %, max\|γ\| = 9.0e-4 (μ = 10.3) |
| Gulec-Whittaker, cracking 3√f'c A_W | `BBCGEN,0,4,4,60,0.005,0,0,0,0` | converged after 5 iterations (max \|dE/E\| 7.3, 6.1, 6.0, 5.6, 5.3, 1.7 %; the damping changed by 3.6 % at analysis 4): E/E_el = 0.72, ξ = 13.5 %, μ = 1.8 |

These variants show the points of §8. In force control the converged state depends strongly on the
damping, and the convergence rate depends on the shape of the backbone near the solution. The
convergence test compares two successive analyses, not the solution itself.

Tests: `tests/unit/test_nonlinear_example.py` runs the example and checks the clean run, the backbone
curves, the convergence, the cracked first storey, the next HOUSE deck and the WRITE → INP round trip.
