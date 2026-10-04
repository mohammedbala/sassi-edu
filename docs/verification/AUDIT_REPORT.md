# SASSI-EDU final audit report

**Package:** `final_audit` (docs/internal/wave3_packages.md). **Date:** 2026-10-02. **Code version:** 0.1.0.

**Auditor:** an independent reviewer who did not write the code. The review covered cross-module physics
and numerics, the verification suite, and the features added last:

- incoherency, wave passage and multiple excitation;
- 2D plane strain and SYMM planes;
- nonlinear-soil SSI iterations;
- TSHELL and the node optimizer;
- Option NON;
- Option A LOADGEN;
- SOIL-NON;
- water modelling;
- AFWRBAT.

**Bottom line:**
- No defect was found that makes a verified result wrong in the paths the verification problems (VPs) exercise.
- Conventions are consistent from SITE through ANALYS, MOTION, RELDISP, STRESS and LOADGEN. This covers units, gravity, complex modulus, damping assignment, phase and time delay, control point and ground elevation.
- The independent checks of section 2 confirm this to round-off.
- The audit found and fixed **four** program defects (section 3):
  - an inconsistent EXCSTRCHK rule between CHECK and HOUSE, with a misleading claim of exactness;
  - a silent wrong-node trap in MOTION/RELDISP after the node optimizer;
  - a gap in CHECK's EDU-06 that let a singular model pass CHECK;
  - Option NON reading the wrong nodes' histories after the node optimizer.
- The audit also corrected **eleven** weaknesses in the verification suite (section 4). These were tautological or blind references, a loop bug, tolerances looser than achieved, NaN-hiding comparisons and stale texts. The remaining points are reported with recommendations.
- The full test suite passes (section 7).

---

## Contents

1. [Scope and method](#1-scope-and-method)
2. [Independent end-to-end checks](#2-independent-end-to-end-checks)
3. [Program defects found and fixed](#3-program-defects-found-and-fixed)
4. [Audit of the verification suite](#4-audit-of-the-verification-suite)
5. [Lead decisions: observations](#5-lead-decisions-observations)
6. [Files changed](#6-files-changed)
7. [Test-suite results](#7-test-suite-results)
8. [Residual limitations](#8-residual-limitations)

---

## 1. Scope and method

**Reading.** I read the following:
- `docs/ARCHITECTURE.md`;
- requirements §§0.3, 2.6, 4.4, 6.3 and 7.16–7.18;
- manual §1.5.1 rules 3–11, the application guidelines 13–14 and §§9.7–9.8;
- the code paths of a full seismic run, from AFWRITE and SITE through POINT, HOUSE, ANALYS (`ssi_solver`, `flexibility`, `freefield`), MOTION, RELDISP, STRESS and LOADGEN;
- the new-feature modules where they meet that chain:
  - `incoherency` (FILE77, FFL/FFM);
  - `symmetry` / `ssi2d`;
  - `nlsoil` restarts;
  - the `renumber` / optimizer numbering;
  - the LOADGEN unit conversions.

**Independent checks.** I designed checks that use as little of the program's own helper code as possible:
- closed forms written out in the test files;
- exact invariances: unit scaling, translation, time-shift;
- an independent oscillator for the response spectra;
- mutation tests.

All of them run through the module decks or through the command language: Interpreter, CHECK, AFWRITE and RUN, which is the GUI path (L17). They are permanent regression tests in `tests/unit/test_audit_*.py` (35 tests).

**Verification-suite audit.** Four independent read-only reviews covered all 74 registered VPs and their pytest wrappers. They:
- classified every `check`, `require` and `inform` as independent, own-code (tautology risk), mis-transcribed or mislabelled;
- recomputed the closed forms;
- ran mutation experiments.

I verified each finding I acted on myself before changing anything.

**Rules followed:**
- minimal fixes, each covered by a regression test;
- no reference value or tolerance weakened. One tolerance was *tightened* (VP-48), and several criteria were made stronger or independent (VP-TS1, VP-H2, VP-04, VP-46, VP-I1; section 4);
- lead decisions respected; observations are in section 5;
- `docs/` (except this file) and `sassi/ui/` not edited.

---

## 2. Independent end-to-end checks

All results are from runs on 2026-10-02. "rel" is the maximum relative difference.

| # | Check (test) | Modules on the path | Reference (independent of the code) | Result |
|---|---|---|---|---|
| E1 | **SI vs British units.** The same model (2-layer site, surface SHELL mat, 2-mass stick) runs in m–kN–t and in ft–kip with g = 9.80665/0.3048 (`test_audit_endtoend.py::test_units_*`) | SOIL (3 iterations), SITE, POINT, HOUSE, ANALYS, MOTION, STRESS, RELDISP, LOADGEN (static and dynamic) | Exact unit scaling. Dimensionless results must be identical; dimensional results must scale exactly. A hard-coded 9.81, or a module that uses another gravity, would break it | FILE8 translations 1.7e-11; rotations ×0.3048 exact to 2e-11. `.ACC` / `.RS` (g) 1e-11. `.THD` ×1/0.3048 2.5e-11. Beam forces and moments ×1/4.448 and ×1/(4.448·0.3048) ≤ 3e-8. SOIL `ACC`/`SN` identical; `SS` ×0.0209. LOADGEN `F` 2e-12, `D` 4e-12, ACEL/D tables 6e-12 |
| E2 | **Wave-passage delay in the time domain** (`test_wave_passage_delays_the_motion_downstream`) | HOUSE WPASS → FILE77 → ANALYS FFM → MOTION | With e^{+iωt}, a later arrival is e^{−iωx/V}. Nodes 40 m and 100 m downstream (V = 400 m/s) must be delayed by exactly 20 and 50 samples | H = e^{−iωx/V} to 1.2e-12. Histories equal the shifted reference history to 7.5e-15. The opposite sign gives 100 % error |
| E3 | **Inclined SV apparent velocity** (`test_inclined_sv_arrives_later_downstream`) | SITE (30°, inclined body wave) → ANALYS → MOTION | Delay x·sinθ/V_hs: 4 and 10 samples | 9e-14 (FILE8) and 1.5e-15 (histories) |
| E4 | **In-structure response spectra** (`test_isrs_equals_an_independent_frequency_domain_oscillator`) | MOTION RS (Nigam–Jennings, sub-stepping) | Exact absolute-acceleration transfer function of each 5 % oscillator applied to the zero-padded `.ACC` history; ZPA = PGA | ≤ 0.47 % over 0.1–100 Hz; SA(100 Hz) within 0.4 % of the PGA |
| E5 | **Translation invariance** (`test_translation_of_model_ground_and_control_point_leaves_the_solution_unchanged`) | AFWRITE, SITE, POINT, HOUSE, ANALYS | An embedded FV excavation with mat and stick under inclined SV (k ≠ 0), moved by (100, −50, 7) m together with `gelev` and the control point, must give the same FILE8 | 4.5e-14 |
| E6 | **Damping assignment and complex-modulus form**, β_p = 0.03 ≠ β_s = 0.08, CMODFORM 0 and 1 (`test_audit_materials.py`) | SITE (SV, P), SOIL (SH, vertical P), HOUSE (L-table excavated soil, M type 1 E–ν, M type 3 Vp–Vs), AFWRITE | Within-layer cos(k* z); SHAKE amplification 1/cos(k* H); u^T K u = G* γ² V and M* ε² V for affine fields; r^T M r = ρV, with c(β) written out | SITE 7.5e-10 (the mixed mass removes the O((kh)²) error of vertical waves), against 5e-3 to 2e-2 for swapped dampings. SOIL 2.4e-7 (limited by the 6-decimal frequencies of the SAF file). HOUSE 1e-15. Both forms pass |
| E7 | **Zero-SSI identity through the command path** (`test_zero_ssi_identity_runs_through_the_command_path`) | Interpreter, CHECK, AFWRITE, SITE, POINT, HOUSE, ANALYS | U = U′_f at all 27 nodes, compared with raw FILE1 data, not the code's free-field function | < 1e-8. Before the audit this model could not reach HOUSE: CHECK reported an error. See F-01 |
| E8 | **RELDISP vs double integration** (scratch; covered by VP-34) | MOTION, RELDISP | Frequency-domain double integration of the full-period acceleration difference | 2.3e-14. Double-integrating the truncated 24 s `.ACC` output instead gives a peak 16 % too large (maximum difference 23 % of the peak), because truncation introduces drift. Use RELDISP, as the manual says |

Already covered by existing tests, so not duplicated:
- **Reciprocity** of the complex-symmetric SSI system through FORCE + ANALYS: `tests/unit/test_review_analys.py`.
- **SDOF on a rigid disk** against an exact 3-DOF model with exact-kernel BEM impedances: VP-17.
- **SOIL vs SITE consistency** on a strain-softened profile: `tests/integration/test_examples.py::test_ex04_site_reproduces_the_soil_amplification`.
- **WRITE → INP model-state identity** for every example: `same_state` on the full canonical model.

---

## 3. Program defects found and fixed

### F-01 (high): EXCSTRCHK. CHECK and HOUSE disagreed, and HOUSE called an inaccurate model "exact"

**Finding.**
- **CHECK** reported every excavation-interior node shared with a structural element as `Error EXCSTRCHK`. That error blocks the HOUSE deck.
- **HOUSE** accepted every such node that was an interaction node, with the warning *"accepted: the flexible-volume general assembly is exact there"*.
- **Consequences:**
  - The zero-SSI identity (VP-16, VP-T3, required by §6.3 with "same mesh") and near-field soil on the excavation mesh ran from decks, but were refused through the command language.
  - A flexible structure on interior excavation nodes was accepted by HOUSE and called exact.
  - VP-52 (CHECK must report an error) and `test_house.py` (HOUSE must only warn) enforced the contradiction.

**Manual.**

| Source | What it says |
|---|---|
| §9.8.1 | The condition is "incorrect from a SASSI modeling point of view", and "Any node reported … may cause incorrect SSI analysis results" |
| Rule 4 | The excavated soil connects with the structure only at the foundation–soil interface, or only with the *near-field soil* elements |
| Rule 11 | "the FV method provides often close results for SSI models with … unique mesh in the basement … limited to … stiff walls and floors. For flexible structural systems placed in the basement … can affect significantly the accuracy", and FI methods give "very poor results" |
| Application guideline 13a | "none of the internal excavation volume nodes shall be connected to a structural node" |

**Numerical evidence** (scratch runs). A 6 m × 6 m box, 2 layers, 4 × 4 mesh. "Shared" means the structure uses the interior excavation nodes; "separate" means it has its own nodes there. The FI shared case was forced through HOUSE by a test patch.

| Configuration | Shared vs separate, max relative difference of the ATF |
|---|---|
| FV, stiff SOLID block on the excavation mesh | 0.46 % |
| FV, soft (soil-like) SOLID block | 2.5 % |
| FI-EVBN, stiff SOLID block | 5.0 % |
| FI-EVBN, soft SOLID block | 32 % at the top; 478 % at mid-depth (\|H\| 1.085 instead of 0.244 at 20 Hz) |
| FV, flexible intermediate floor (SHELL, t = 0.1 m) at mid-depth, vertical P | separate mesh: resonance peak 9.74 at 15.0 Hz. Shared: no resonance (max \|H\| = 1.0); the floor is tied to the soil through the residual of X_ff − (K*_e − ω²M_e) |
| Same, t = 0.05 m | separate 15.0 at 7.6 Hz; shared 1.85 at 17 Hz |
| Same, t = 0.4 m | max deviation 32 % of the peak |

So the HOUSE claim was wrong for flexible subsystems. The CHECK blanket error was stricter than the manual: the manual's EXCSTRCHK lists nodes, and FV tolerates soil-like elements on the mesh. It also contradicted the zero-SSI requirement.

**Decision.** One rule, in `sassi/core/house_lib.py::excstrchk_kind`, used by both CHECK and HOUSE:

- **Error** when the shared interior node is not an interaction node. In FI the structure is tied to −(K*_e − ω²M_e) without the impedance.
- **Error** when any structural element at the node is a BEAMS, SHELL, TSHELL, SPRING or GENERAL element, a K node, or a SOLID/PLANE element not built on the excavation mesh. These are flexible subsystems.
- **Warning** when the node is an interaction node and *every* structural element there is a SOLID/PLANE whose nodes all belong to the excavation. This covers near-field or backfill soil, a solid basement on the excavation mesh, and the zero-SSI identity, for which FV is close or exact.
- The non-linear-soil exemption of CHECK (no message for PINGRP groups at interaction nodes) is kept.
- The new HOUSE warning text says the result is *close*, not exact, and that separate nodes are the correct model.
- The EXCSTRCHK command now labels each listed node `[CHECK error]` or `[CHECK warning]`.

This refines G-15 / requirements §4.4 item 3 ("EXCSTRCHK condition (error)") for the soil-replacement case; section 5 asks the lead to ratify it.

**Tests:** `tests/unit/test_audit_excstrchk.py`, 7 tests:
- the rule table;
- CHECK and HOUSE classify five configurations identically (SOLID FV/FI, pile FV/FI, floor FV);
- the zero-SSI identity through the command path.

`tests/unit/test_house.py::test_excstrchk_and_edu22` was updated: a pile on an interior FV node is now an error in HOUSE, as in CHECK. VP-52 is unchanged and passes.

### F-02 (medium): node optimizer. Post-processing requests silently named other nodes

**Finding.**
- With Optimize Model (`HOUSEX,1`), FILE4 and FILE8 use the new numbers. The manual has the user select MOTION/RELDISP nodes from the `.hounew`.
- ANALYS and LOADGEN translate model numbers automatically, but MOTION and RELDISP do not.
- Requirement §2.6 asks for a warning ("HOUSE optimizer used: MOTION/RELDISP node requests are in the new numbering"), but none existed.
- **Evidence:** example 1 with `HOUSEX,1` reported the mat nodes 37 and 78–81 for the requested 41 and 82–85, which are the stick floors. The ISRS of a "top floor" was a mat corner. The run ended "MOTION finished with status OK; 0 warning(s)".

**Fix.**
- `motion.optimizer_numbering_note` (also called by RELDISP) runs when `<model>.map` exists. HOUSE renames a stale map `.prev`.
- It lists every requested node whose model number differs, as `new = model`, and every request absent from the map.
- The manual's semantics (new numbers) are unchanged.

**Tests:** `tests/unit/test_audit_optimizer.py` (MOTION and RELDISP list the mapping; no note without the optimizer).

### F-03 (medium): EDU-06. CHECK missed free drilling rotations of shells on SOLID elements

**Finding.**
- CHECK's EDU-06 considered only nodes with shells and springs.
- Consider a basemat or wall SHELL on excavated SOLID nodes (the usual embedded model). Its drilling rotation is unrestrained, because the SOLID has no rotations.
- CHECK passed. HOUSE warned with "use FIXROT/FIXSHLROT", but FIXROT (manual §9.7) treats shell-only nodes only. ANALYS then stopped on a singular system.

**Fix.**
- `check.drilling_unrestrained` now also covers shell nodes shared with SOLID/PLANE elements. Its message asks for D, or for a rotational spring when the shell is oblique.
- The HOUSE hint now mentions D for such nodes.
- FIXROT itself is unchanged; it follows the manual.

**Test:** `test_audit_endtoend.py::test_check_reports_free_drilling_rotations_of_shells_on_solids`.

### F-04 (medium): Option NON after the node optimizer read other nodes' histories

**Finding.**
- NONLINEAR reads the relative-displacement histories `nnnnnTR_X.THD` of its panel corners and spring ends by their **model** numbers. Those numbers come from the `.eql` deck, and the `<model>_new.hou` it writes uses them too.
- After a renumbering by the optimizer, MOTION and RELDISP write the files under the **new** numbers. NONLINMOTDISP requests the model numbers, which MOTION reads as new numbers.
- A file named after a model number therefore holds another node's history, and NONLINEAR would compute panel strains and equivalent-linear properties from the wrong nodes without any message.
- A run of example 7 with `HOUSEX,1` completed with 0 warnings. That optimizer kept the identity numbering, so it was harmless there; any model the optimizer actually renumbers is exposed.

**Fix.** NONLINEAR stops with a clear message when `<model>.map` is not the identity: run HOUSE without the optimizer for Option NON. An identity map is accepted.

**Test:** `tests/unit/test_audit_optimizer.py::test_option_non_refuses_a_renumbered_model`.

---

## 4. Audit of the verification suite

**Result:**
- No VP passes without testing its main claim. All 74 pass.
- **No reference constant is mis-transcribed.** The reviews recomputed every closed form and compared every table with R2. Two examples:
  - VP-42 matches Jakub & Roesset (1977), Tables 1–2.
  - The VP-04 Table B-2, PGA profile and transfer-function values match R2 H.1 exactly.
- No pytest wrapper weakens its VP.

The weaknesses below are of three kinds: references computed with the code under test, criteria blind to the defect class their label names, and stale or overstated texts that the verification manual prints.

### 4.1 Corrected (each with a regression test in `tests/unit/test_audit_vpsuite.py`)

| VP | Finding | Correction |
|---|---|---|
| VP-TS1 | **Loop bug.** The MITC4 parallelogram equilibrium patch sat inside the triangle loop and used the leaked loop variable `eint` (= 1). "EINT 1" was recorded three times and EINT 0 was never tested | Moved into the `for eint in (0, 1)` loop. Both pass (9.3e-16) |
| VP-TS1 | The Navier reference took κ from `tshell.SHEAR_FACTOR`, the element's own constant | Literal 5/6 |
| VP-TS1 / VP-54T | The reference resultants come from `tshell.plate_rigidities`, the constants the element uses. The membrane rigidity had no independent check | New independent test: the element's resultants under constant membrane, shear and transverse-shear states of a distorted quad, against E t/(1−ν²), G t, (5/6) G t written out |
| VP-H2 | "G* = ρVs² c(β_s)" compared \|G*\| of the VP's *own* material object with ρVs². Since \|c(β)\| = 1 it could not see the damping. A conjugated or dropped damping passed 11/11 | G* and M* are now recovered from HOUSE's Ke with affine fields and compared in complex form with c(β) written out (1e-15). The mutation test fails as it should |
| VP-04 | The damping reference used `shake.interp_log_strain`, the function SOIL uses | Independent `np.interp` in log10 strain. A linear-in-strain mutation of SOIL is now detected |
| VP-04 | The D-W1-06 spectral criterion rows were labelled "diagnostic", and the notes said "VP-04 SA checks fail" | Labels "(D-W1-06)" and corrected notes. Wrapper updated; a dead filter and an unused xfail text were removed |
| VP-46 | The BBCGEN cracking-force reference was the program's own `ShearCapacities.cracking` | 3√f′c·A_W written out (1832.82 kips) |
| VP-I1 | The title called the reference "independent". It is built with `flexibility_matrix` and `free_field_at_nodes`, the functions ANALYS uses. Its f→0 and monotonicity criteria were applied to the VP's own reference instead of the ANALYS result. An info label ("EPRI: FFL slightly conservative") contradicted its own numbers | Title and docstring state what is compared. The criteria apply to the ANALYS result. The label is neutral |
| VP-E1, VP-28, VP-37, VP-LA1/LA2 | Worst-error accumulations with `e > worst` or the built-in `max` drop NaN. The framework's `worse` exists for exactly this | NaN-propagating |
| VP-48 | The brute force evaluates every vertex and both window ends, so it is exact, yet the tolerance was 2e-3·max(y) (observed 1.2e-14) and its note was wrong | **Tightened** to 1e-9·max(y); note corrected |
| Texts printed by the manual | VP-26 docstring ("the X and Y criteria fail"), VP-06 (0.5 % cut-offs, convergence from above), VP-H1 ("independent element-library assembly"; it is an identity of the deck translation), VP-38T module text ("VP-38T fails", "strict xfail"), `NAFEMS_8X8_NOTE` ("lead decision requested"), VP-51 ("48 hexes … failing"), VP-09 note | Rewritten to the lead decisions D-W3-05, D-W1-03, D-W3-01, D-W2-05 and D-W1-04 that the code implements |

### 4.2 Reported, not changed (recommendations)

**Medium:**
- **VP-36:** the SRP 3.7.1 criteria use the same evaluator EQUAKE uses to select its motion. A reviewer's fully independent re-evaluation reproduced every verdict:
  - min SA/target 0.952 / 0.934 / 0.930;
  - max 1.132–1.182;
  - longest run below target 6 / 3 / 7;
  - PSD ≥ 0.962;
  - ρ = −0.130 / −0.061 / +0.034.

  So nothing is wrong today, but the VP cannot detect an evaluator defect. Recommendation: build an independent checker into the VP.
- **VP-30:** it lists MOTION, EQUAKE and SOIL but exercises only `core.spectra`. E4 above now covers MOTION's RS path. EQUAKE `.rso` and SOIL `RS` are still covered only through their use of the same routine.
- **VP-43:** the Lamb ν = 1/4 values of the requirement row are not checked (the title says so). A reviewer confirmed the R2 values. Recommendation: implement the check, or record a decision deferring it.
- **VP-28 (d):** the option-6 ISRS criterion uses 152 frequencies against 58 for options 0–5, because the spline needs a denser grid (spec 01 check C11 recommends at least 200 frequencies with option 6). In a reviewer's rerun with 58 frequencies, option 6 errs 11.4 % against the 2 % criterion. The label states the frequency count, but no lead decision records this.
- **VP-16, VP-T3 (by design):** the identity compares with the same `free_field_at_nodes` ANALYS uses. It cannot detect a HOUSE-vs-SITE material error common to the structure and the excavated soil. That error class is now covered by E6 and the new VP-H2 checks; VP-S2 already detects density and damping-sign errors.
- **VP-26, VP-I1:** the X and Y factor sets are identical (isotropic horizontal coherency), so an X/Y swap of factor sets in ANALYS would go unnoticed.

**Low:**
- Loose tolerances relative to the achieved accuracy, but equal to the requirement table:
  - VP-15: 1 % against 5e-10;
  - VP-S1: 1 % against an exact factor known to 1e-4;
  - VP-26 commutator: 1e-6 against 3.7e-10;
  - VP-02b: compared at 5/3 Hz instead of R2's printed 1.6667 Hz.
- Rows that cannot fail on a program error:
  - VP-39b: 5.033 Hz recomputed from the VP's own k and m;
  - VP-07: symmetry after explicit symmetrisation;
  - VP-27: trace rows;
  - VP-53: its own synthetic CRITFREQ data;
  - VP-30: 1/(2ζ) against 10.000;
  - VP-36: grid and SD unit rows;
  - VP-54W: LISTPOOLINTER "0 of 0" would pass.
- Own-code references:
  - VP-37: SOLID/PLANE patch references use `isotropic_D` / `plane_strain_D`;
  - VP-SN2: secant reference `mkz_stress`;
  - VP-35: material and beam frames;
  - VP-48 T-B3: `Line.at`.
- **VP-NON1:** the "time domain: stabilised loop damping" row re-runs `HY.cyclic_loop` rather than reading NONLINEAR's time-domain output.
- **VP-LA2:** a D table whose RELDISP history is identically zero is not compared.
- **VP-25:** uses synthetic FILE81/FILE82; the real ANALYS → COMBIN path is covered by a unit test.
- **VP-E1:** HOUSE `imp`/`incomp` and SITE `hs`/`hslaw` are not compared; FILE8 equality catches them only where they change results.

---

## 5. Lead decisions: observations

None of these changes a criterion. Each needs the lead's attention.

1. **EXCSTRCHK (F-01).**
   - Requirements §4.4 item 3 and G-15 call the condition an error.
   - The audit's rule makes it a warning only for interaction nodes used solely by SOLID/PLANE elements on the excavation mesh. That case is required by the VP-16 "same mesh" zero-SSI identity (§6.3) and by near-field soil.
   - Proposed decision text: *"EXCSTRCHK severity follows `house_lib.excstrchk_kind`: error unless the node interacts and every structural element at it is a SOLID/PLANE on the excavation mesh (warning)"*, with the evidence table of F-01.
2. **D-W3-07 (water).**
   - The decision says VP-W1 gives the impulsive base shear "within 1.5 % of exact potential flow, mesh-convergent".
   - Observed (rerun by the auditor): +1.46 % (2 Hz), +1.54 % (4 Hz), +1.86 % (8 Hz), for both the SHELL and the SOLID tank. The VP criterion is 5 % (`W1_RTOL`, also in the §6.3 row), while the pytest wrapper adds a one-sided 3 % bound. VERIFY and pytest therefore apply different criteria.
   - VP-W1 has no mesh-convergence row; the unit test `test_water_physics.py::test_mesh_refinement_converges_to_the_exact_series` has it.
   - Recommendation: state the observed band, align `W1_RTOL` with the wrapper (3 %), and cite the convergence test.
3. **D-W1-05 (VP-38).** These numbers come from a reviewer's independent, converged Legendre–Ritz computation:
   - FV16: λ₁ = 3.4710, so NAFEMS 0.421 Hz is **+0.73 %** high, not "≈ 0.6 %".
   - FV12: NAFEMS mode 4 (4.233 Hz, ×2) is +1.02 % above the exact value (Leissa), and mode 5 is +0.82 %. The SHELL 32 × 32 mode 4 is −0.37 % from the exact value but −1.38 % from NAFEMS, which leaves only 0.62 % of the 2 % tolerance.
   - No change is needed, but the decision text should give these numbers.
4. **Verification manual.** `docs/verification/VERIFICATION_MANUAL.md` (generated 17:45 UTC) is stale:
   - it contains 56 of the 74 VPs;
   - it reports VP-38T as failing, but VP-38T passes under D-W3-01.

   It must be regenerated (`python -m sassi.verify.report`). The corrected VP texts of section 4.1 then flow into it.
5. **R2 D.5 typo.** The Love dispersion equation prints `√(1/c²−1/β1²)`; it should be `√(1/β1²−1/c²)`. The numbers in R2 are right.
6. **Strict xfail.** `tests/unit/test_review_site_point.py::test_vp06_first_love_cutoff_frequency` still says VP-06 "does not check" the 11.547 Hz cut-off. VP-06 checks it with the provisional 10 % of D-W1-03 (observed 12.40 Hz). The xfail correctly documents that 0.5 % is not met; only its reason text is outdated.

---

## 6. Files changed

**Program:**
- `sassi/core/house_lib.py`: `excstrchk_kind`, the shared EXCSTRCHK severity rule (F-01).
- `sassi/modules/house.py`:
  - EXCSTRCHK uses `excstrchk_kind`, with a new error message for interaction nodes and a corrected warning (F-01);
  - the EDU-06 hint mentions D (F-03).
- `sassi/prep/check.py`:
  - EXCSTRCHK severity through `excstrchk_severity` / `excstrchk_kind`, and the module docstring (F-01);
  - `drilling_unrestrained` covers shell nodes on SOLID/PLANE elements (F-03).
- `sassi/prep/commands/checks.py`: the EXCSTRCHK command labels each node with its CHECK severity (F-01).
- `sassi/modules/motion.py`: `optimizer_numbering_note` and its call (F-02).
- `sassi/modules/reldisp.py`: optimizer note for RDND / RELFILE nodes (F-02).
- `sassi/modules/nonlinear.py`: `_refuse_renumbered_model` (F-04).

**Verification problems (section 4.1):**
- `sassi/verify/problems/vp_tshell.py`
- `vp_house.py`
- `vp_soil.py`
- `vp_nonlinear.py`
- `vp_incoherency.py`
- `vp_site_point.py`
- `vp_generation.py`
- `vp_examples.py`
- `vp_motion.py`
- `vp_elements.py`
- `vp_loadgen.py`
- `vp_linemath.py`

**Existing tests updated:**
- `tests/unit/test_house.py`: the EXCSTRCHK assertion (F-01).
- `tests/verification/test_vp_soil.py`: VP-04 label check; dead filter and unused xfail text removed.
- `tests/verification/test_vp_generation.py`: comments.

**New:**
- `tests/unit/test_audit_excstrchk.py` (7)
- `test_audit_endtoend.py` (8)
- `test_audit_optimizer.py` (3)
- `test_audit_materials.py` (6)
- `test_audit_vpsuite.py` (11)
- this report.

**Concurrent changes by other agents.** These were not made by this audit and were not reviewed in it:
- `docs/` (manuals, requirements, command reference);
- `sassi/ui/`;
- `sassi/modules/loadgen.py`, `sassi/prep/registry.py` and `sassi/prep/commands/soilnon_cmds.py`, modified at 17:32–17:33 local time.

The final test run below includes the state of all files at the time it ran.

---

## 7. Test-suite results

| Run | Command | Result |
|---|---|---|
| Baseline, before any change | `.venv/bin/python -m pytest -q` | 2690 passed, 1 xfailed, 256 warnings in 526 s |
| Final, after all fixes | `.venv/bin/python -m pytest -q` | **2776 passed, 1 xfailed**, 256 warnings in 514 s |
| Audit tests alone, final versions | `.venv/bin/python -m pytest -q tests/unit/test_audit_*.py` | 35 passed |

The xfail is the strict xfail of item 6 in section 5. All 74 VPs pass, including VP-38T.

The count grew by 86 tests. This audit added 35 of them; the other agents working at the same time added the rest.

---

## 8. Residual limitations

These are known limits of the verification, stated plainly:

- **POINT with β_p ≠ β_s.** There is no closed form for the point-load flexibility. E6 covers SITE, SOIL and HOUSE; POINT inherits the FILE2 layer moduli of SITE, and the FILE2 meta carries CMODFORM.
- **Restart reuse.** A New-Structure restart (ANALYS mode 1) reuses X_ff when the interaction nodes, layer table, gravity and CMODFORM are unchanged (FILE90). A change of the POINT central-zone radius or of the SITE half-space model (`nl`, `hslaw`) is **not** detected. This follows G-23 and the manual, and remains the user's responsibility.
- **FIXROT.** Following manual §9.7, FIXROT treats only shell-only and solid-only nodes. Shells on excavated SOLIDs need D. CHECK now says so (F-03).
- **Optimizer numbering.** Following the manual, MOTION and RELDISP requests are new numbers. F-02 adds a warning listing the mapping but does not translate, and Option NON refuses renumbered models (F-04) rather than translating. ANALYS (FORCE loads) and LOADGEN translate automatically.
- **EXCSTRCHK.** The warning class (FV, SOLID/PLANE on the excavation mesh) is "close", not exact, unless the elements reproduce the excavated soil. The audit measured 0.5 % for a stiff block and 2.5 % for a soft block. Users should still model separate nodes (manual rule 11).
- **Not re-derived by the auditor:**
  - the Abrahamson coherency coefficients of `sassi/data/coherency/*.json` (models 2–6; not checked against the primary EPRI sources in this audit);
  - the Cheng–Mertz and Takeda hysteresis rules (7.15 item 5);
  - SOIL-NON beyond VP-SN1/SN2;
  - the full TSHELL formulation beyond patch, closed-form and NAFEMS checks.
- **LOADGEN APDL.** The APDL output is verified by a parser and an in-code replay (VP-LA2); it has **not** been run in ANSYS (D-W3-11).
- **Not checked by CHECK:**
  - G-14 (structural interaction nodes below grade outside the excavation): HOUSE warns;
  - SOIL's `grav` against the HOUSE gravity. A consistent mismatch changes only SOIL stresses and FILE88 G values, not the wave solution.
- **Strain-compatible soil in embedded models.** With SITEX soil mode 1, SITE uses the FILE88 properties of SOIL, but the excavated soil of HOUSE keeps the L table. HOUSE warns that the L layers must be updated, but no command copies FILE88 into the L table: it is a manual step.
- **File precision.** The SOIL SSAF `.TFU` prints frequencies with 6 decimals, which limits text-file comparisons to about 1e-6.
- **Incoherency.** In the strongly incoherent limit the deterministic AS (algebraic-sum) factors depend on the eigenvector basis of the (nearly) degenerate coherency matrix. D-W3-06 makes them reproducible; for Σ = I and 4 nodes, s = 1.37, 1.03, 0.51, −0.90 with Σ s² = N. This is a property of the AS approximation, which the manual restricts to rigid foundations (EDU-13). It is not a defect.
- **Nonlinear options.** Nonlinear soil SSI and Option NON are equivalent-linear iterations. SOIL-NON is a 1D time-domain column. None is a general nonlinear time-domain SSI solver; the manual's methods are the same.
- **Status.** SASSI-EDU is educational and not qualified for licensing or design work.
