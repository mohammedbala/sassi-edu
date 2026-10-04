# Foundation impedance and inertial SSI: verification study (VP-10, VP-11, VP-13, VP-14, VP-17)

This report goes with `sassi/verify/problems/vp_impedance.py` and `tests/verification/test_vp_impedance.py`.
It records how the SITE + POINT + HOUSE + ANALYS chain was set up for the foundation-impedance
problems of requirements §6.3, how sensitive the results are to the mesh, the central-zone radius R0
and the soil layering, and why some checks against the R2 references fail.

Numbers are from SASSI-EDU on 2026-10-02. All stiffnesses are normalised by G R (translations) and
G R³ (rotations) for the disk of radius R, or by G B and G B³ for the square of half-width B. Unless
stated otherwise the soil is uniform with Vs = 100 m/s, ρ = 2 t/m³ and ν = 1/3 (Vp = 200 m/s,
G = 20 000 kPa). It is near-elastic, with β = 1e-4.

## 1. Summary

| VP | Result | Failing checks | Cause (evidence below) |
|---|---|---|---|
| VP-10 disk statics | **pass** | — | Extrapolated welded values are within 0.6–2.3 % of the relaxed closed forms and within 0.3–1.4 % of an independent welded BEM. |
| VP-11 disk dynamics | **xfail** | c_r against Veletsos–Verbic at a0 = 0.5 / 1.0 / 1.5 (+95 / +56 / +29 %) | The Veletsos–Verbic rocking damping is too low (by 53 / 41 / 28 %) compared with a rigorous relaxed solution and with the exact low-frequency limit. SASSI-EDU agrees with the rigorous solution to within 1.2 %. All other Veletsos–Verbic checks and all high-frequency dashpots pass. |
| VP-13 Tajirian–Tabatabaie | **xfail** | Vertical resonances 1 and 2 at β = 5 % | The vertical compliance peaks do not sit at the 1-D column frequencies (2n−1)Vp/(4H). For the first, the layer has a zero-group-velocity onset at A0 = 0.509 (continuum P-SV dispersion; SITE 0.511), and damping shifts the peak further down to 0.485. The second is hidden by 5 % damping (it is at 1.585, +0.9 %, when β = 1 %). The horizontal resonances, the third vertical resonance and the statics all pass. |
| VP-14 square statics | **xfail** | K_xx = K_yy against Pais–Kausel (+6.9 %) | Pais & Kausel's rocking fit is 7 % below the welded exact-kernel BEM value (6.457) for L = B. It is also below the relaxed BEM value (6.237). SASSI-EDU is 0.7 % below the BEM. K_z, K_x, K_y and K_zz pass. |
| VP-17 inertial SSI | **xfail** | Peak \|u_t/u_g\| against 6.71 (−17 %) | Two causes. First, R2's 6.71 comes from Veletsos–Verbic impedances, whose rocking damping is too low; rigorous relaxed impedances give 6.19. Second, the welded mat adds sliding–rocking coupling (K_xθ = −0.075 K_x R; the BEM gives −0.081), which lowers the peak to 5.56. F.2 fed with the rigorous welded impedance (exact-kernel BEM) gives 5.563, as SASSI-EDU does. The F.2 identity holds to 8e-10. Statics, T̃/T and the resonance frequency pass. |

No tolerance was changed. Each failing check keeps its R2 reference and tolerance. Each test asserts
that every other check passes and reports the known failures as xfail with the observed numbers.

**Regression guards of the known failures.** An xfail must not hide a regression of the very
quantity that fails. Two layers of protection:

1. Every failing quantity has a guard check *inside the VP* against a derived reference that passes
   today (these are ordinary checks: if one fails, the test fails).
2. The pytest wrapper (`tests/verification/test_vp_impedance.py`) stores the documented value of
   every known failure with a band; a known failure whose computed value leaves the band fails the
   test instead of reporting xfail.

| known failure | documented value | band | guard check inside the VP (derived reference) | observed |
|---|---|---|---|---|
| VP-11 welded c_r at a0 = 0.5 / 1 / 1.5 | 0.0458 / 0.1246 / 0.1859 | ±3 % | relaxed c_r vs the relaxed exact-kernel BEM (10 %) | −1.2 … +0.6 % |
| VP-13 vertical resonance 1 (β = 5 %) | 0.4848 | ±1 % | peak below the continuum ZGV onset 0.5094, and within 6 % of it | −4.8 % |
| VP-13 vertical resonance 2 (β = 5 %) | none (NaN) | — | peak of the β = 1 % sweep vs 3Vp/(4H) (2 %) | 1.5855, +0.9 % |
| VP-14 K_xx = K_yy | 6.4136 | ±1 % | vs the welded exact-kernel BEM (5 %) | −0.7 % |
| VP-17 peak \|u_t/u_g\| | 5.563 | ±2 % | peak and resonance vs F.2 fed with the welded exact-kernel BEM impedance (5 %) | +0.004 %, +0.2 % |

## 2. Model rules

* **Interaction mesh.** Interaction nodes sit on the ground surface: a ring mesh for the disk (ring k
  carries 6k nodes, element size h = R/n) or a square grid (h = 2B/n). Stiff massless BEAMS connect
  them to the centre node (`builders.surface_rigid_mat`). K_G = Tᵀ X_ff T depends on X_ff only, so the
  mat does not change K_G.
* **Central zone (R1 §3.4).** R0 = 0.85 h for ring meshes and 0.90 h for square meshes.
* **Soil column.** The top sublayer thickness is h. That is the cube-like near-surface element that the
  R0 rules are calibrated for. Sublayers then grow by 12 % per sublayer, up to min(λ_min/8, 4R), down
  to 30 R (statics) or 20 R (dynamics). Below that come 20 generated half-space sublayers (UNIFORM law,
  D-W1-01) and base dashpots. The high-frequency band (a0 = 6–8) uses sublayers ≤ λ_s/8 down to 3 R.
* **Welded and relaxed contact.**
  * ANALYS `<impe>` = 2 gives K_G = Tᵀ X_ff T (D-ANL-07). Here every translation of the interaction
    nodes follows the rigid body, so K_G is a **welded**-contact impedance.
  * **Consistency checks of K_G** (`kg_checks`, in every VP: VP-10/11/14 at every frequency, VP-13 at
    three frequencies of each sweep and in the statics, VP-17 for the seismic run). K_G is compared with
    Tᵀ F_ff⁻¹ T formed by an independent LU solve (`numpy.linalg.solve`), not with the
    `ssi_solver.impedance_matrix` / `global_impedance` path that ANALYS itself uses; and K_G must be
    complex symmetric (reciprocity). Observed: 1e-16 to 2e-15, and symmetry to 1e-16. The flexibility
    F_ff comes from the same `flexibility_matrix`, so these checks cover ANALYS's inversion and
    condensation, FILE11 and the node ordering, not POINT; POINT is covered by the physics references.
  * The **relaxed**-contact impedance (no shear traction under vertical and rocking motion, as assumed
    by Boussinesq and by Veletsos–Verbic) comes from the same F_ff. Vertical and rocking use the
    vertical block, Kᵣ = T_zᵀ F_zz⁻¹ T_z. Horizontal and torsion use the horizontal block.
* **Mesh extrapolation.** Three self-similar meshes are used: R0/h and t0/h are fixed, so the
  horizontal and vertical discretisations refine together. Richardson extrapolation with the formal
  order p = 1 uses the two finest meshes. The observed order from all three meshes is reported. It is
  1.1–2.3 wherever the sequence is monotone.
* **FILE11.** ANALYS runs with FILE11 holding K_G and T only (`file11_without_matrices`).
  Otherwise X_ff of every frequency would be written (up to 100 MB per frequency here). The switch
  is `analys.file11_x_cap(0)`, a context variable: only ANALYS runs of the calling thread are
  affected (an ANALYS run in another thread of the same process still stores X_ff, as D-ANL-07
  requires). It no longer patches the module global.
* **Disk use.** Once K_G, F_ff and the transfer functions are in memory, each run deletes its binary
  intermediate files (FILE1/2/3, FILE8, FILE11, HOUSE matrices; `purge_work_files`) and keeps only
  the decks and listings. A run used to leave 7–75 MB per VP; it now leaves 0.1–8 MB (VP-13: the
  300-frequency ANALYS listings). `KEEP_WORK_FILES = True` keeps everything for debugging.

## 3. Independent references built for this study (derived, type D)

The R2 closed forms are approximate for several quantities: welded versus relaxed contact, and
curve fits. To decide whether a failed check is a SASSI-EDU error or a reference error, two
independent rigorous solutions were coded in `vp_impedance.py`. They are unit-tested in
`tests/unit/test_impedance_helpers.py` (and probed independently in `tests/unit/test_review_impedance_vps.py`).

1. **Static rigid-foundation BEM** (`static_bem`). It uses the exact Boussinesq/Cerruti surface
   kernels, piecewise-constant tractions on quadrilateral cells and collocation at the centroids. The
   singular self-cell integrals are evaluated in polar coordinates. Error is O(h); the two finest
   meshes are extrapolated.
   * Disk, ν = 1/3, n = 16/24 rings:
     * relaxed: K_v 6.0011 (exact 6.000) and K_r 4.0023 (exact 4.000);
     * welded: K_v 6.1310 (Mossakovskii exact 4GR ln(3−4ν)/(1−2ν) = 6.1299), K_h 4.8419,
       K_r 4.1479, K_t 5.3365 (Reissner–Sagoci 5.3333), K_xθ −0.394.
   * Square, n = 24/32 cells per side, welded: K_z 7.0570, K_x 5.5860, K_zz 8.5883 and K_xx 6.4572.
     The relaxed K_xx is 6.2380.
2. **Dynamic BEM with exact Lamb kernels** (`lamb_kernels`, `dynamic_bem_impedance`,
   `relaxed_bem_impedance`, `welded_bem_impedance`).
   * The surface Green's functions of the half-space come from Hankel transforms of the Lamb
     symbols. The integration runs on a contour lifted into the upper half plane, so it passes
     above the branch points and the Rayleigh pole (time factor e^{+iωt}). The kernels reproduce the
     Wong (1975) / Apsel (1979) table of R2 D.3 to 0.001, including the vertical–horizontal
     coupling kernel G_xz (column U_r0, all 11 values) that welded contact needs.
   * The smooth dynamic remainders are tabulated and splined in r. The table now runs to an upper
     bound of the largest source–observation distance (`max_pair_distance`: 2R for a disk,
     2√2 B for a square); it used to stop at 2.2 max|V|, too short for a square. For disks the
     change moves the results by 5e-9.
   * Relaxed contact: the disk results, extrapolated from n = 16/24, are stored in `VP11_BEM`
     (regenerating the a0 = 1.0 entry reproduces it to 4 digits).
   * Welded contact (`welded_disk_impedance`, n = 8/12 extrapolated; 12/16 agrees within 0.3 %): the
     rigorous counterpart of K_G, with the dynamic sliding–rocking coupling (at a0 = 1.07:
     K_xθ = (−0.527 + 0.074i) G R², −0.112 Re K_x R). It is the reference of the VP-17 peak guard.
3. **Continuum dispersion of a layer on a rigid base** (`layer_dispersion_det`, `layer_zgv_onset`).
   The P-SV equations of a homogeneous layer are written as a first-order system y′ = A y; the
   propagator expm(A H) gives the secular function (free top, fixed base) without the thin-layer
   discretisation of SITE. It reproduces the k = 0 column frequencies (2n−1)πV/(2H) exactly. For
   the VP-13 layer the second P-SV mode is born at a ZGV point A0 = 0.5094 (k R = 0.170); SITE with
   24 sublayers gives 0.511 (+0.3 %).
4. **Exact low-frequency radiation damping** (`lowfreq_damping`). This is a first-order perturbation
   with the static punch tractions and the imaginary part of the Lamb kernel (body waves plus the
   Rayleigh pole). For the relaxed disk at ν = 1/3 it gives:
   * c_h(0) = 0.580;
   * c_v(0) = 0.783;
   * c_r = 0.240 a0².

## 4. VP-10: rigid disk statics on a half-space

R = 5 m at a0 = 0.05. The meshes have 8, 12 and 16 rings (217, 469 and 817 nodes).

| | n = 8 | n = 12 | n = 16 | extrapolated | R2 B.1 (relaxed) | err | welded BEM | err |
|---|---|---|---|---|---|---|---|---|
| K_v / GR | 6.0991 | 6.0977 | 6.1000 | **6.1069** | 6.0000 | +1.8 % | 6.1310 | −0.4 % |
| K_h / GR | 4.8892 | 4.8659 | 4.8562 | **4.8272** | 4.8000 | +0.6 % | 4.8419 | −0.3 % |
| K_r / GR³ | 4.1590 | 4.1296 | 4.1200 | **4.0913** | 4.0000 | +2.3 % | 4.1479 | −1.4 % |
| K_t / GR³ | 5.5348 | 5.4487 | 5.4093 | **5.2912** | 5.3333 | −0.8 % | 5.3365 | −0.8 % |

The relaxed-contact values from the same F_ff, extrapolated, are K_v 5.976 and K_r 3.948 (exact 6
and 4). The welded-minus-relaxed gap, 2–4 %, is the one R2 B.1 anticipates ("welded vs relaxed
contact").

**Sensitivity to R0 and to the top sublayer t0.** Each cell shows the n = 8 / 12 / 16 values, then
the extrapolated value.

| t0, R0 | K_v | K_h | K_r | K_t |
|---|---|---|---|---|
| h, 0.75 h | 5.976 / 6.006 / 6.026 → 6.087 | 4.823 / 4.817 / 4.817 → 4.817 | 3.977 / 3.990 / 4.005 → 4.052 | 5.414 / 5.356 / 5.334 → 5.266 |
| **h, 0.85 h** | 6.099 / 6.097 / 6.100 → **6.106** | 4.889 / 4.866 / 4.856 → **4.827** | 4.159 / 4.130 / 4.120 → **4.091** | 5.535 / 5.449 / 5.409 → **5.291** |
| h, 0.95 h | 6.226 / 6.191 / 6.175 → 6.126 | 4.955 / 4.915 / 4.895 → 4.837 | 4.353 / 4.276 / 4.240 → 4.129 | 5.652 / 5.538 / 5.482 → 5.316 |
| h/2, 0.85 h | 6.015 / 6.036 / 6.051 → 6.095 | 4.791 / 4.796 / 4.802 → 4.818 | 4.029 / 4.033 / 4.042 → 4.070 | 5.339 / 5.308 / 5.298 → 5.269 |
| h/4, 0.85 h | 5.975 / 6.006 / 6.027 → 6.088 | 4.725 / 4.748 / 4.764 → 4.809 | 3.968 / 3.987 / 4.005 → 4.058 | 5.214 / 5.215 / 5.223 → 5.247 |

* **Effect of R0 and t0 on a single mesh.** On one mesh, R0 and t0 change the stiffness by up to
  ±3 % (K_r at n = 16 runs from 4.005 to 4.240). They calibrate the near-singular diagonal of F_ff,
  and that diagonal is an O(h) contribution.
* **After extrapolation.** All variants fall within ±1 % of one another and within 2.5 % of the
  welded BEM.
* **Order of convergence.** It is clean (p ≈ 1.1) for R0 = 0.95 h. For smaller R0 it is irregular,
  because the tributary-area error and the central-zone error have opposite signs.
* **Column depth and growth rate.** The 0.85 h ring rule, with t0 = h, is the best overall choice.
  Doubling the depth of the column (60 R instead of 30 R), or using a 6 % instead of a 12 % growth,
  changes the stiffness by less than 0.1 %.

## 5. VP-14: square footing statics

B = 5 m at a0 = ωB/Vs = 0.05. The grids have 16, 24 and 32 cells per side.

| | n = 16 | n = 24 | n = 32 | extrapolated | Pais & Kausel | err | welded BEM | err | Gazetas |
|---|---|---|---|---|---|---|---|---|---|
| K_z / GB | 7.1546 | 7.1162 | 7.0990 | **7.0475** | 7.05 | −0.04 % | 7.0570 | −0.1 % | 6.81 |
| K_x = K_y / GB | 5.7200 | 5.6707 | 5.6477 | **5.5790** | 5.52 | +1.1 % | 5.5860 | −0.1 % | 5.40 |
| K_zz / GB³ | 9.2236 | 8.9902 | 8.8787 | **8.5441** | 8.31 | +2.8 % | 8.5883 | −0.5 % | 8.35 |
| K_xx = K_yy / GB³ | 6.7710 | 6.6446 | 6.5869 | **6.4136** | 6.00 | **+6.9 %** | 6.4572 | −0.7 % | 5.40 / 5.58 |

* **Symmetry.** K_x = K_y and K_xx = K_yy hold to 1e-16.
* **Sensitivity.** The extrapolated rocking K_xx is 6.357, 6.414 and 6.472 for R0 = 0.8 h, 0.9 h and
  1.0 h with t0 = h, and 6.393 for t0 = h/2. The order is p ≈ 1.0–2.3. The +6–8 % deviation from
  Pais & Kausel therefore does not depend on the discretisation choices.
* **Evidence on the reference.** The same BEM reproduces the exact disk values to 0.1 %, and
  SASSI-EDU's square agrees with the BEM to within 0.7 %.
  * At L/B = 1 the Pais–Kausel rocking coefficient, 4.0/(1−ν) = 6.0, is 7 % below the welded value
    and 4 % below the relaxed value.
  * The two published fits disagree with each other by 10 % for this quantity: Gazetas gives 5.40
    and 5.58.
  * R2 B.2 itself states that the fits were "not re-verified here".
  * **Recommendation:** replace the VP-14 rocking reference with the welded BEM value 6.46 GB³ (type D),
    or widen it to the fit's real accuracy (≈ ±8 %). This is a decision for the lead.

## 6. VP-11: rigid disk dynamic impedance

Disk R = 5 m, 12 rings. The coefficients are k = Re K / K0 and c = Im K / (a0 K0), where K0 is the
stiffness of the same model at a0 = 0.05.

| a0 | | k_h | c_h | k_r | c_r | k_v | c_v | k_t | c_t |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | SASSI-EDU welded (K_G) | 0.994 | 0.586 | 0.944 | 0.0458 | 0.982 | 0.789 | 0.956 | 0.0321 |
| | SASSI-EDU relaxed | 0.991 | 0.590 | 0.940 | 0.0507 | 0.971 | 0.791 | 0.956 | 0.0321 |
| | relaxed BEM, exact kernels | 0.990 | 0.585 | 0.938 | 0.0501 | 0.969 | 0.793 | 0.956 | 0.0310 |
| | Veletsos–Verbic | 1 | 0.65 | 0.953 | **0.0235** | 0.952 | 0.789 | — | — |
| 1.0 | SASSI-EDU welded | 0.974 | 0.592 | 0.840 | 0.1246 | 0.921 | 0.805 | 0.866 | 0.0931 |
| | SASSI-EDU relaxed | 0.967 | 0.603 | 0.829 | 0.1359 | 0.893 | 0.824 | 0.866 | 0.0931 |
| | relaxed BEM | 0.965 | 0.599 | 0.825 | 0.1353 | 0.883 | 0.824 | 0.867 | 0.0913 |
| | Veletsos–Verbic | 1 | 0.65 | 0.840 | **0.0800** | 0.863 | 0.859 | — | — |
| 1.5 | SASSI-EDU welded | 0.949 | 0.606 | 0.747 | 0.1859 | 0.823 | 0.838 | 0.783 | 0.1465 |
| | SASSI-EDU relaxed | 0.946 | 0.621 | 0.731 | 0.2018 | 0.788 | 0.877 | 0.782 | 0.1465 |
| | relaxed BEM | 0.941 | 0.617 | 0.724 | 0.2005 | 0.765 | 0.872 | 0.783 | 0.1440 |
| | Veletsos–Verbic | 1 | 0.65 | 0.712 | **0.1440** | 0.793 | 0.915 | — | — |

**Agreement with the rigorous solution.** SASSI-EDU (relaxed contact) matches the rigorous relaxed
BEM to within 3 % for every coefficient (c_r within 1.2 %).

**Where Veletsos–Verbic fall short.**

* **Rocking damping.** The approximations follow the rigorous rocking stiffness k_r to within 2 %,
  but their rocking damping c_r is 53 / 41 / 28 % too low. At low frequency the exact limit is
  c_r = 0.240 a0², against 0.1 a0² for Veletsos–Verbic. SASSI-EDU's relaxed c_r at a0 = 0.25 is
  0.0148 − 2β/a0 = 0.0140, against 0.0150 exact.
* **Horizontal damping.** α1 = 0.65 is a band average. The exact low-frequency value is c_h(0) =
  0.580, and the rigorous c_h rises from 0.585 to 0.617 over 0.5 ≤ a0 ≤ 1.5. R2 C.2 already notes
  that k_h is "truly mildly frequency-dependent". The welded c_h passes at 10 % (−9.8 / −8.9 / −6.8 %).

**Mesh sensitivity.** In an earlier run with t0 = h/2, the welded c_r at a0 = 0.5 was 0.0449–0.0450
for 8, 12 and 16 rings.

**High-frequency dashpots** (Im K / ω divided by ρVsA, ρVpA, ρVpI, ρVsJ; R2 C.1). The values below are
for 12 rings with sublayers ≤ λ_s/8, down to 3 R.

| a0 | horizontal | vertical | rocking | torsion |
|---|---|---|---|---|
| 6 | 0.979 | 0.995 | 1.032 | 0.973 |
| 7 | 0.973 | 0.975 | 1.029 | 0.980 |
| 8 | 0.956 | 0.963 | 0.999 | 0.979 |

All values are within the 5 % tolerance.

* **Mesh and depth sensitivity** (t0 = h/2 runs). With 16 rings instead of 12, or a 6 R column, or
  λ/12 sublayers, every value moves by at most 2 %.
* **Finite a0.** The exact solution itself approaches the plane-wave limit from below at finite a0,
  and the cone model predicts about −2 % for rocking at a0 = 6. Part of the 2–4 % shortfall is
  therefore physical.

## 7. VP-13: Tajirian & Tabatabaie (1985)

* **Model.** r = 0.5, H = 3, Vs = 1, Vp = 2, ρ = 1, on a rigid base. There are 24 sublayers of 0.125 and
  4 rings, so h = 0.125 = the sublayer thickness and R0 = 0.106. The frequency grid is A0 = 0.01–3.0
  in steps of 0.01.
* **Resonance definition.** A resonance is a local maximum of |C_ii|, where C = K_G⁻¹, within ±10 %
  of the expected value. Its position is refined with a parabola.

| resonance | expected A0 | β = 5 % | err | β = 1 % (diagnostic) | β = 15 % |
|---|---|---|---|---|---|
| horizontal 1 | 0.262 | 0.2627 | +0.4 % | 0.261 | 0.262 |
| horizontal 2 | 0.785 | 0.7911 | +0.7 % | 0.790 | none |
| horizontal 3 | 1.309 | 1.2917 | −1.3 % | 1.291 | none |
| vertical 1 | 0.524 | **0.4848** | **−7.4 %** | 0.492 and 0.508 | none |
| vertical 2 | 1.571 | **none** | — | 1.585 (+0.9 %) | none |
| vertical 3 | 2.618 | 2.6486 | +1.2 % | 2.651 | none |

* **Robustness of the vertical peak.** The first vertical peak at β = 5 % stays at A0 = 0.4845–0.4848
  for 4 or 8 rings, 24 or 48 sublayers, and R0 = 0.75 h or 0.85 h. It is a property of the
  layer–disk system, not of the discretisation.
* **What the layer's modes show.** The second Rayleigh mode, the thickness-stretch family, starts
  to propagate at A0 = 0.5094 with k R = 0.170 > 0 (continuum P-SV dispersion, section 3; SITE Mode 1
  with β = 1e-4 gives A0 = 0.511, checked within 1 %). That is a zero-group-velocity (backward-wave)
  point below the k = 0 cut-off Vp/(4H), as is usual for layers with ν < ~0.45.
* **Consequence for the reference.** Radiation starts at that point, and the damped compliance peak
  lies below it. The (2n−1)Vp/(4H) values of R2 C.7 are 1-D column frequencies; they are not
  compliance peaks of the 3-D problem. The second vertical resonance exists (1.585 at β = 1 %), but
  5 % damping erases it.
* **Guards.** The VP checks that the β = 5 % peak lies below the continuum ZGV onset and within 6 % of
  it (observed −4.8 %), and that the second vertical peak of the β = 1 % sweep is at 3Vp/(4H) within
  2 % (1.5855, +0.9 %). All K_G of the sweeps pass the consistency checks of section 2.
* **Recommendation:** check the vertical resonances on the lightly damped run, with a ZGV-aware
  reference for the first one, or against the digitised figures (requirements §6.7). This is a
  decision for the lead.

**Statics.** These use |K_G| at A0 = 0.01; |c(β)| = 1 in the SASSI damping form, so this is the
elastic stiffness. The meshes are 4/24, 8/48 and 12/72 (rings/sublayers), with h equal to the
sublayer thickness.

| | 4/24 | 8/48 | 12/72 | extrapolated | reference | err |
|---|---|---|---|---|---|---|
| K_h / GR | 5.3723 | 5.2751 | 5.2492 | **5.1975** | 8GR/(2−ν)(1+R/2H) = 5.2000 | −0.05 % |
| K_r / GR³ | 4.2925 | 4.1612 | 4.1330 | **4.0767** | 8GR³/(3(1−ν))(1+R/6H) = 4.1111 | −0.8 % |
| K_v / GR | 7.1836 | 7.1427 | 7.1428 | 7.1431 | [VERIFY] 4GR/(1−ν)(1+1.28R/H) = 7.28 | −1.9 % |
| K_t / GR³ | 5.8224 | 5.5333 | 5.4480 | 5.2774 | [VERIFY] 16GR³/3 = 5.333 | −1.0 % |

## 8. VP-17: SDOF on a rigid massless disk

* **Model.** R = 10 m with 9 rings. The mat links are E = 1e6 G BEAMS. The SDOF has h = 10 m,
  m = 942.48 t and T = 0.5 s (k = 148 830 kN/m), with β = 5 % hysteretic and CMODFORM 1, so
  k* = k(1 + 2iβ) as in R2 F.2.
* **Identity.** FILE8 of the seismic run equals the generalised F.2 solution fed with the code's own
  impedance S at the centre node, to 7.8e-10. Here S is the inverse of the 2×2 compliance from two
  FORCE load cases on the mat alone. F.2 is generalised to coupled S_xθ.
* **Effect of the links.** S differs from the ideal rigid-mat K_G by 1.4e-3, which is the flexibility
  of the links.

| quantity | SASSI-EDU | reference (R2 F.1) | err |
|---|---|---|---|
| static K_x (kN/m) | 9.765e5 | 9.60e5 | +1.7 % (welded) |
| static K_θ (kN m) | 8.299e7 | 8.00e7 | +3.7 % (welded; welded BEM: +3.7 %) |
| Veletsos–Meek T̃/T with the code's springs | 1.1540 | 1.15805 | −0.3 % |
| resonance of \|u_t/u_g\| | 1.693 Hz | 1.717 Hz | −1.4 % |
| peak \|u_t/u_g\| | **5.56** | 6.71 | **−17 %** |

**Decomposition of the peak.** The R2 values (1.7175 Hz, 6.709) are reproduced exactly from F.2 with
Veletsos–Verbic impedances.

| impedance fed to F.2 | resonance (Hz) | peak |
|---|---|---|
| Veletsos–Verbic (R2) | 1.717 | 6.709 |
| rigorous relaxed BEM (VP11_BEM) | 1.709 | 6.191 |
| SASSI-EDU relaxed contact | 1.712 | 6.176 |
| SASSI-EDU welded, coupling removed | 1.719 | 6.334 |
| rigorous welded BEM (`welded_disk_impedance`), coupling removed | 1.717 | 6.335 |
| rigorous welded BEM (`welded_disk_impedance`) | 1.689 | 5.563 |
| SASSI-EDU welded (the actual model) | 1.693 | 5.563 |

* **Rocking damping.** It accounts for −8 %: Veletsos–Verbic underestimate it (VP-11).
* **Sliding–rocking coupling of the welded mat.** It accounts for a further −12 %. The static
  K_xθ = −0.075 K_x R is confirmed by the welded BEM (−0.081). The coupling raises the foundation
  flexibility seen by a horizontal force at height h by about 8 %.
* **Guard.** The VP checks the SASSI peak and resonance against F.2 fed with the rigorous welded
  impedance (5 %): 5.5633 vs 5.5631 and 1.693 vs 1.689 Hz. Agreement this close is partly
  coincidental; term by term the 9-ring K_G is within 2–4 % of the extrapolated BEM.
* **Recommendation:** derive the VP-17 worked-example reference from rigorous welded impedances
  (the guard value 5.563 above), or restrict the comparison to the resonance frequency. This is a
  decision for the lead.

## 9. Runtimes and disk use (Mac mini, one process)

| VP | runtime | marked slow | files left in the work directory |
|---|---|---|---|
| VP-10 | 8 s | no | < 0.5 MB |
| VP-11 | 10 s | no | < 0.5 MB |
| VP-13 | 40 s | no | 7.6 MB (mostly the 300-frequency listings) |
| VP-14 | 9 s | no | < 0.5 MB |
| VP-17 | 92 s | yes (`slow=True`) | 0.4 MB |

The static BEM references take about 5 s in total (VP-10 and VP-14); the welded dynamic BEM of the
VP-17 guard about 9 s, the continuum ZGV onset of VP-13 about 3 s. The real-axis part of the
wavenumber integrals now uses real-argument Bessel functions, which made the kernel tables about
ten times faster.
