# SASSI-EDU Theory Manual

This manual explains the methodology that SASSI-EDU implements, with the equations and enough of the
derivations for a structural engineer to follow them. It is the equivalent of the theory basis of a
commercial SSI code. Every section names the code that implements the equations and the
verification problems (VPs) that check them; the [Verification Manual](../verification/VERIFICATION_MANUAL.md)
gives the numbers. The normative sources are the SASSI theory note
[`docs/spec/R1_sassi_theory.md`](../spec/R1_sassi_theory.md) (R1), the requirements
[`docs/spec/00_requirements.md`](../spec/00_requirements.md) (§4 methodology, §7 decisions D-xxx) and
the benchmark collection [`docs/spec/R2_benchmarks.md`](../spec/R2_benchmarks.md) (R2).

> SASSI-EDU is an educational re-implementation, not affiliated with or endorsed by the vendor of
> ACS SASSI. The original SASSI reports (UCB/GT/81-02, SASSI2000 theory manual) were not accessible
> to this project; the formulation is rebuilt from the open literature by the SASSI authors and
> their collaborators (R1 §11). Where a step had to be reconstructed this manual says so.

**Notation.** Matrices are bold in the text where useful and plain in the formula blocks; `*` marks
complex (damped) quantities; `ᵀ` is the transpose (never the conjugate transpose: the SSI matrices
are complex *symmetric*); `i = √-1`; `ω = 2πf`.

## Contents

1. [Conventions](#1-conventions)
2. [Complex-modulus damping](#2-complex-modulus-damping)
3. [Flexible-volume substructuring](#3-flexible-volume-substructuring)
4. [The thin-layer method (SITE Mode 1)](#4-the-thin-layer-method-site-mode-1)
5. [Half-space: variable depth and viscous boundary](#5-half-space-variable-depth-and-viscous-boundary)
6. [Free-field motion and control-point normalisation (SITE Mode 2)](#6-free-field-motion-and-control-point-normalisation-site-mode-2)
7. [Point-load solutions: central zone and transmitting boundary (POINT)](#7-point-load-solutions-central-zone-and-transmitting-boundary-point)
8. [Flexibility and impedance assembly](#8-flexibility-and-impedance-assembly)
9. [Solution of the SSI equation (ANALYS)](#9-solution-of-the-ssi-equation-analys)
10. [Transfer-function interpolation](#10-transfer-function-interpolation)
11. [Convolution, response spectra, relative displacements and stresses](#11-convolution-response-spectra-relative-displacements-and-stresses)
12. [The SHAKE equivalent-linear method (SOIL)](#12-the-shake-equivalent-linear-method-soil)
13. [Element formulations (HOUSE)](#13-element-formulations-house): the element library, the TSHELL thick shell, the node optimizer
14. [Incoherency, wave passage and multiple excitation](#14-incoherency-wave-passage-and-multiple-excitation)
15. [Two-dimensional SSI and symmetry planes](#15-two-dimensional-ssi-and-symmetry-planes)
16. [Nonlinear soil SSI iterations](#16-nonlinear-soil-ssi-iterations)
17. [Option NON: nonlinear structures by equivalent linearisation](#17-option-non-nonlinear-structures-by-equivalent-linearisation)
18. [SOIL-NON: nonlinear time-domain site response](#18-soil-non-nonlinear-time-domain-site-response)
19. [Water: hydrodynamic mass by soft solids](#19-water-hydrodynamic-mass-by-soft-solids)
20. [Option A: SSI loads for a second-step ANSYS analysis](#20-option-a-ssi-loads-for-a-second-step-ansys-analysis)
21. [Spectrum-compatible motions (EQUAKE)](#21-spectrum-compatible-motions-equake)
22. [Code and verification map](#22-code-and-verification-map)
23. [References](#23-references)

---

## 1. Conventions

| Item | Convention | Source |
|---|---|---|
| Harmonic time factor | `exp(+iωt)` | R1 §0, D-CNV-01 |
| Horizontal propagation | `exp(-ikx)`; outgoing waves have Im k < 0, or Re k > 0 when k is real | R1 §0 |
| Cylindrical outgoing waves | Hankel functions of the second kind `H(2)μ(kρ)` | R1 §3.2 |
| Vertical axis | **z up**; layer interfaces numbered 1 (ground surface) downwards | R1 §0 |
| In-plane thin-layer variables | `{u_x, i·u_z}` ("i-scaled vertical", Kausel form): all layer matrices symmetric; physical `u_z = -i·φ_z` | R1 §2.2 |
| FFT | forward `A_k = Σ a_n exp(-2πikn/N)` (`numpy.fft.rfft`), inverse with 1/N (`irfft`) | D-CNV-02 |
| Frequency grid | `Δf = 1/(Δt·NFFT)` (or SITE `<fstep>`), SSI frequencies `f = n·Δf` with integer n | requirements §4.0.1 |
| Transfer functions | seismic: `H = U/U_cp` (total motion per unit control motion, dimensionless, the same for displacement and acceleration); vibration: displacement per unit load factor | D-CNV-06 |

With `exp(+iωt)` a time delay τ is the factor `exp(-iωτ)`, a viscous dashpot c contributes `+iωc` to
the dynamic stiffness and a hysteretic material has a modulus with a *positive* imaginary part. The
transfer functions of ANALYS can therefore be multiplied directly by `numpy.fft.rfft` of the input.

---

## 2. Complex-modulus damping

**Code:** `sassi/conventions.py` (`cfactor`, `complex_lame`), `sassi/elements/base.py`
(`material_from_M`, `material_from_layer`), `sassi/core/tlm.py` (`layer_moduli`),
`sassi/core/shake.py` (`complex_velocity`). **Verified by:** VP-01, VP-02a, VP-04, VP-39, VP-39b.

### 2.1 The SASSI form

Soil and structural damping are represented as **hysteretic** (rate-independent) damping, through a
complex modulus. The Berkeley codes (SHAKE after the 1973 Udaka-Lysmer revision, SHAKE91, TLUSH,
SASSI2000) use

```
G* = G · c(β),     c(β) = 1 - 2β² + 2iβ·√(1-β²)                                  (2.1)
```

`c(β)` has modulus 1: writing `c = exp(iδ)` gives `cos δ = 1 - 2β²`, `sin δ = 2β√(1-β²)`, i.e.
`δ = 2·arcsin β`. Damping rotates the modulus in the complex plane without changing its magnitude.
The complex wave velocity is

```
V* = √(G*/ρ) = V·(√(1-β²) + iβ),      |V*| = V                                    (2.2)
```

and a vertically propagating wave `exp(-ik*z)` with `k* = ω/V* = (ω/V)(√(1-β²) - iβ)` decays as
`exp(-βωz/V)`: per radian of phase the amplitude decreases by the factor `exp(-β/√(1-β²)) ≈ exp(-β)`.

Evidence that SASSI uses this form (R1 §1.1): for the hysteretic SDOF transfer function the recovered
damping `1/(2|H|max)` is 4.99/9.95/14.83/19.60 % for β = 5/10/15/20 % with Eq. (2.1), against
4.98/9.81/14.37/18.57 % for `1 + 2iβ`; Ostadan et al. (2004) published 5.0/9.9/14.7/19.5 % for
SASSI2000. Only Eq. (2.1) reproduces the 19.5 % at β = 20 %.

The simpler form `c = 1 + 2iβ` (original SHAKE 1972) is available with `CMODFORM,1` for benchmarks;
its modulus `√(1+4β²)` slightly stiffens the material. Damping ratios β ≥ 0.5 are rejected (EDU-04):
the real part `1 - 2β²` of Eq. (2.1) vanishes at β = 1/√2, and such damping ratios have no physical
meaning for soils or structures long before that.

### 2.2 Which modulus carries which damping

The shear modulus carries the S-wave damping and the **constrained** modulus `M = λ + 2G = ρVp²`
carries the P-wave damping (D-CNV-04):

```
G* = G·c(βs)      M* = M·c(βp)      λ* = M* - 2G*                                   (2.3)
```

| Element / layer | Complex stiffness |
|---|---|
| soil layers (SITE, POINT), excavated soil, SOLID, PLANE | Lamé constants `(λ*, G*)` from Eq. (2.3) |
| BEAMS | `E*, G*` from `(M*, G*)` (`ν* = (M*-2G*)/(2(M*-G*))`, `E* = 2G*(1+ν*)`); `E* = E·c(β)` when βp = βs |
| SHELL, TSHELL | `E* = E·c(β)`, ν real (CHECK warns when βp ≠ βs, EDU-27) |
| SPRING | `k* = k·c(damp)` |
| GENERAL | `K* = MXR + i·MXI` as entered |

Masses are always real. Because every element stiffness is linear in the moduli, an element whose
βp = βs has `K* = c(β)·K₀`.

### 2.3 The hysteretic SDOF and the meaning of β

A mass m on a spring `k* = k·c(β)` whose base moves with unit amplitude has the total-motion transfer
function

```
H(ω) = k* / (k* - mω²)                                                             (2.4)
|H|² = 1 / [ (1 - 2β² - r²)² + 4β²(1-β²) ],      r² = mω²/k
```

The minimum of the denominator is at `r² = 1 - 2β²`, which gives the peak

```
|H|max = 1 / (2β·√(1-β²))     at  ω = ω₀·√(1-2β²)                                 (2.5)
```

(25.0050, 10.0125, 5.0252, 3.3715, 2.5516 for β = 2, 5, 10, 15, 20 %). VP-01 runs exactly this system
through SITE, POINT, HOUSE, ANALYS and MOTION and checks Eq. (2.5) to 1e-6; VP-39/39b do the same
with SPRING and GENERAL elements.

**Equivalent ANSYS input** (R1 §1.2): ANSYS structural damping multiplies the stiffness by `(1 + ig)`.
`E_ANSYS = E(1 - 2β²)` and `g = 2β√(1-β²)/(1 - 2β²)` give the same complex modulus.

---

## 3. Flexible-volume substructuring

**Code:** `sassi/core/ssi_solver.py`, `sassi/modules/analys.py`, `sassi/modules/house.py`.
**Verified by:** VP-16 (zero-SSI identity), VP-15, VP-E1, VP-S2, the impedance problems VP-10 to
VP-17, and example 2 (FV against FI).

### 3.1 The idea

The soil-structure system is split into three parts (Lysmer et al. 1981; Tabatabaie 1982):

```
   SSI system  =  free field  +  structure  -  excavated soil
```

* the **free field**: the horizontally layered site *without* excavation, which SITE and POINT solve
  semi-analytically (thin-layer method); it is represented by its dynamic stiffness `X_ff` at the
  *interaction nodes* and by its motion `U'_f` there;
* the **structure**: the finite-element model of the building *including* its basement and any
  near-field soil (HOUSE matrices `K*_s`, `M_s`);
* the **excavated soil**: a finite-element model of the soil volume that the basement replaces, with
  the free-field layer properties (HOUSE matrices `K*_e`, `M_e`). It is subtracted, because the free
  field already contains that soil.

### 3.2 Derivation

Let `C(ω) = K* - ω²M` denote a dynamic stiffness. Consider the body "structure minus excavated
soil": its dynamic stiffness is `C^s - C^e`, assembled on one set of nodes (the structure and the
excavation share the nodes of the basement-soil interface). The rest of the world is the free field
with the excavated volume still in place but "cancelled" by the negative stiffness. The free field
has the free-field motion `U'` when nothing else acts on it; imposing a different displacement `U`
on its interaction nodes requires the forces `X_ff (U - U')`, where `X_ff` is the free-field dynamic
stiffness (impedance) condensed to those nodes. By action and reaction the free field exerts
`-X_ff (U - U')` on the body, and the body's equilibrium (no external force for a seismic analysis)
is

```
(C^s - C^e) U = -X_ff (U - U')     ⇒     [C^s - C^e + X_ff] U = X_ff U'              (3.1)
```

`U` are **total** displacement amplitudes. With the degree-of-freedom sets of the manual, `s`
structure only, `i` shared by structure and excavation (basement interface), `w` excavation only,
and `f` interaction DOFs, Eq. (3.1) is manual Eq. 2.1 (R1 §4.4):

```
| C^s_ss    C^s_si                    0              | | U_s |   | F_s                         |
| C^s_is    C^s_ii - C^e_ii + X_ii    -C^e_iw + X_iw | | U_i | = | X_ii U'_i + X_iw U'_w + F_i |   (3.2)
| 0         -C^e_wi + X_wi            -C^e_ww + X_ww | | U_w |   | X_wi U'_i + X_ww U'_w       |
```

where `X_ab` is the block of `X_ff` (zero when a or b is not an interaction set) and `F` are external
forces (FORCE; `U' = 0` for a pure vibration analysis). SASSI-EDU assembles it in the general form
(D-ANL-01)

```
C(ω) = A_sᵀ(K*_s - ω²M_s)A_s - A_eᵀ(K*_e - ω²M_e)A_e + A_fᵀ X_ff A_f
b    = A_fᵀ X_ff U'_f   (seismic)      or      b = P (vibration, FILE9)                  (3.3)
```

with Boolean gather matrices `A`; the structure and the excavated soil are on one global DOF map
(HOUSE), and `A_f` picks the translations of the interaction nodes.

**Zero-SSI identity** (the first built-in check). If the "structure" is identical to the excavated
soil, `C^s = C^e` and Eq. (3.1) gives `X_ff (U - U') = 0`, i.e. `U = U'` at every node: the free
field is recovered exactly. VP-16 builds such a model (FV, SV, SH and P input) and finds `U = U'` to
round-off (about 2e-15 relative, tolerance 1e-8).

### 3.3 Method variants: the choice of the interaction set

Eq. (3.1) is exact when **every** node of the excavated soil is an interaction node: the
**flexible-volume (FV)** method. The variants use fewer interaction nodes, which makes `X_ff` (a
dense matrix) much smaller:

| Method | Interaction set f | Consequence |
|---|---|---|
| FV (direct) | all excavated-soil nodes | exact (reference) |
| FI-FSIN, subtraction method (SM) | soil-foundation interface: lateral and bottom faces | the `w` rows keep only `-C^e_ww`; the system is near-singular where `det C^e_ww(ω) ≈ 0`, i.e. at the natural frequencies of the excavated volume with the interface nodes fixed: spurious peaks at and above the first excavated-volume frequency (DOE OE-3 2011-02; Mertz et al. 2011) |
| FI-EVBN, modified subtraction (MSM) | FSIN plus the ground-surface face of the excavation | the surface constraint raises the excavated-volume frequencies; much more accurate, anomalies possible above the new frequency |
| FFV (Ghiocel 2013) | EVBN plus internal horizontal levels every `skip` levels | the sub-volumes between interaction levels are small and stiff; close to FV at a fraction of the cost |

In SASSI-EDU the variants differ **only** by the interaction set (INT/INTGEN); the equations and the
code are the same. Example 2 shows the SM anomaly at 6 Hz on a 10 × 10 × 5 m box and its removal by
FI-EVBN. Screening frequencies used in practice (DOE/STP 2011): the soil-layer frequency `Vs/(4H)`
for embedment H and the excavated-volume frequency `f_EV`.

---

## 4. The thin-layer method (SITE Mode 1)

**Code:** `sassi/core/tlm.py` (`column_matrices`, `rayleigh_modes`, `love_modes`,
`modal_flexibility`), `sassi/modules/site.py`. **Verified by:** VP-05 (Rayleigh phase velocity),
VP-06 (Love dispersion), VP-43 (monotone convergence), VP-02b, and through POINT VP-07 to VP-09.

### 4.1 Discretisation

The site is a stack of horizontal layers on a rigid base or a visco-elastic half-space. Each layer is
divided into thin **sublayers** in which the displacement varies **linearly** with depth, while the
horizontal dependence stays analytic, `exp(iωt - ikx)` (Lysmer and Waas 1972; Kausel 1981). This is a
finite-element discretisation in z only, so the accuracy depends on the sublayer thickness relative
to the wavelength: the 1/5-wavelength rule `h ≤ Vs/(5 f_cut)` holds for the **mixed mass**
`M = ½ M_consistent + ½ M_lumped` used by SASSI-EDU in SITE, POINT and the excavated soil (D-CNV-05).
At `h = λ/5` the mixed mass keeps the 1-D amplification within about 9 %, whereas the consistent mass
errs by up to 17 % (R1 §2.1); the mixed mass makes the shear-wave dispersion error of fourth order
(R2 A.3, VP-03).

### 4.2 Layer matrices (Kausel form)

For a sublayer of thickness h with Lamé constants `λ*`, `G*` (Eq. 2.3) and density ρ, in the DOF order
`{u_x1, u_y1, ũ_z1, u_x2, u_y2, ũ_z2}` (1 top, 2 bottom interface, `ũ_z = i·u_z`, z up), Kausel's
Table 1 gives (R1 §2.2; dots are zeros):

```
        h | 2(λ+2G)    .     .    λ+2G      .     .  |          1 |  .      .    λ-G     .      .   -(λ+G) |
A_m =  ---|    .      2G     .     .        G     .  |   B_m = ---|  .      .     .      .      .     .    |
        6 |    .       .    2G     .        .     G  |          2 | λ-G     .     .     λ+G     .     .    |
          |  λ+2G      .     .   2(λ+2G)    .     .  |            |  .      .    λ+G     .      .  -(λ-G)  |
          |    .       G     .     .       2G     .  |            |  .      .     .      .      .     .    |
          |    .       .     G     .        .    2G  |            |-(λ+G)   .     .   -(λ-G)    .     .    |

        1 |  G     .      .     -G     .      .    |
G_m =  ---|  .     G      .      .    -G      .    |        M_m^cons = (ρh/6) [2I 1I; 1I 2I]
        h |  .     .    λ+2G     .     .   -(λ+2G) |        M_m^lump = (ρh/2) [ I  0; 0  I]
          | -G     .      .      G     .      .    |        M_m      = ½ M_m^cons + ½ M_m^lump
          |  .    -G      .      .     G      .    |                 = ρh [5/12 I  1/12 I; 1/12 I  5/12 I]
          |  .     .  -(λ+2G)    .     .    λ+2G   |
```

The sublayer matrices are assembled by overlapping their blocks at shared interfaces, exactly as in
finite-element assembly, and give the interface load-displacement relation in the
frequency-wavenumber domain

```
P̄ = (A k² + B k + G - ω² M) Ū                                                     (4.1)
```

The direction sub-blocks used below are `A_x` (with λ+2G), `A_z = A_y` (with G), `G_x = G_y` (with G),
`G_z` (with λ+2G), the x-row/z-column coupling block `B_xz`, and `C_• = G_• - ω²M`. The
i-scaled vertical variable makes every matrix complex symmetric; the Waas/SASSI form (vertical
displacement down, no i, a skew-symmetric B entering as `iBk`) is algebraically identical
(`T⁻¹ L_Waas T = L_Kausel` with `T = diag(1, i, ...)`, R1 §2.3).

A rigid base deletes the bottom-interface DOFs; a half-space adds the dashpots of section 5.

### 4.3 Eigenproblems

Setting `P̄ = 0` gives the wavenumbers k of the free waves of the discrete layered medium at a given
frequency.

**Rayleigh (P-SV, in-plane).** The quadratic eigenproblem `(A k² + B k + C) φ = 0` is linearised
*without* doubling its size by the substitution `Z = {φ_x; k φ_z}` (Kausel 1981, Eq. 16):

```
( k² | A_x     0  |  +  | C_x   B_xz | ) | φ_x   |
     | B_xzᵀ  A_z |     |  0    C_z  |   | k φ_z |  = 0                              (4.2)
```

a generalised linear eigenproblem in k² of size 2N_f (N_f free interfaces). SASSI-EDU solves it as
the standard eigenproblem of `-Ā⁻¹ C̄` (Ā is block lower triangular, so the inverse is two
symmetric solves). Each root `k_j = √(k_j²)` is taken on the branch with **Im k_j < 0**, or Re k_j > 0
when k_j is real (decaying or outgoing in +x), and `φ_z = Z_lower/k_j`.

**Love (SH, anti-plane).** `(A_y k² + C_y) φ_y = 0`, N_f modes.

**Normalisation** (Kausel 1981, Eq. 22a; the same as Waas'):

```
φ_xᵀ A_x φ_x + φ_zᵀ A_z φ_z + φ_zᵀ B_xzᵀ φ_x / k_j = 1        (Rayleigh)            (4.3)
φ_yᵀ A_y φ_y = 1                                              (Love)
```

With this normalisation the modes are orthogonal (`Y_iᵀ Ā Z_j = 0` for i ≠ j) and the wavenumber-domain
flexibility has the modal expansion (Kausel 1981, Eqs. 49-52)

```
F_xx = Φ_x D Φ_xᵀ,   F_xz = k Φ_x K⁻¹ D Φ_zᵀ,   F_zz = Φ_z D Φ_zᵀ,   F_yy = Φ_y D_L Φ_yᵀ,
D = diag(1/(k² - k_j²)),  K = diag(k_j)                                             (4.4)
```

which SASSI-EDU checks against the direct inverse of Eq. (4.1) to round-off (unit tests of
`tlm.modal_flexibility`, R1 V1). The modes are sorted by decreasing Re k (shortest wavelength
first). FILE2 stores, per frequency, the generated sublayers and `{k_j, φ_x, φ_z}` (2N_f Rayleigh
modes) and `{k_l, φ_y}` (N_f Love modes); one FILE2 serves POINT2 and POINT3, because the
cylindrical and the plane-strain eigenproblems are identical.

**What VP-05 and VP-06 show.** The fundamental Rayleigh velocity of a deep uniform stratum converges
at O(h²) to the exact `c_R/Vs` (0.9194 at ν = 0.25) within 0.5 % at 10 sublayers per wavelength for
ν ≤ 1/3; for ν ≥ 0.45 linear thin layers lock volumetrically and need about 40 sublayers per
wavelength (lead decision D-W1-02; CHECK warns for ν > 0.47). The Love dispersion of a layer over a
half-space is reproduced within 0.5 %, but the variable-depth half-space cannot represent the
grazing field near a Love cut-off frequency (D-W1-03).

---

## 5. Half-space: variable depth and viscous boundary

**Code:** `sassi/core/tlm.py` (`halfspace_depth`, `halfspace_sublayers`, `build_column`,
`vertical_response`, `outcrop_response`). **Verified by:** VP-02b, VP-08, VP-10, VP-11 (radiation
damping), UT-19.

A visco-elastic half-space under the user layers is simulated with two techniques together
(manual §4.1.2 item 20, R1 §2.4):

1. **Variable depth.** SITE adds `nl` sublayers (SITE `<nl>`, 4-20, 0 = rigid base) with the
   half-space properties and a total thickness of one and a half shear wavelengths of the half-space,

   ```
   H_hs(f) = 1.5 λ_s = 1.5 Vs_hs / f                                               (5.1)
   ```

   because the fundamental Rayleigh mode has decayed at that depth. The buffer therefore becomes
   deeper at low frequency and thinner at high frequency.
2. **Viscous boundary.** Lysmer-Kuhlemeyer dashpots per unit area are attached at the bottom of the
   buffer: `c_s = ρV*_s` on the horizontal DOFs and `c_p = ρV*_p` on the vertical one, with the
   *complex* velocities of Eq. (2.2), so that for vertically propagating waves the dashpot is the
   exact impedance of the damped half-space. With `exp(+iωt)` a dashpot adds `+iωc` to the diagonal
   of G at the base interface (in the i-scaled variables the vertical entry is unchanged because the
   i factors cancel).

**Sublayer law.** The manual says only that the thicknesses increase with depth. SASSI-EDU offers
`EDUOPT,HSLAW,UNIFORM` (`h_i = H_hs/nl`, the default), `GEOMETRIC` (`h_i = h_1 q^(i-1)`,
`h_1 = min(h_last·Vs_hs/Vs_last, H_hs/nl)`, q found by bisection so that Σh_i = H_hs) and `LINEAR`;
the non-uniform laws fall back to uniform when a sublayer would exceed λ_s/8. VP-02b compared the
outcrop transfer function of a uniform layer on a half-space with the exact solution: the geometric
law missed the 1 % criterion (1.10 % at nl = 20) and the uniform law passed (0.20 % / 0.95 % at
nl = 20 / 10), so the default is UNIFORM (lead decision D-W1-01, the fallback written into decision
D-SIT-02). R1 V4 had found ≤ 0.7 % (nl = 20) and ≤ 1.5 % (nl = 10) with uniform sublayers, but up to
20 % with nl = 5.

**Within and outcrop motions.** The *within* transfer function between two interfaces above the
half-space does not depend on the half-space model; the half-space model matters for the motion
*below* the control point and for the *outcrop* motion. SASSI-EDU computes the discrete outcrop
motion at the top of the half-space as the surface motion of the half-space-only column (same
generated sublayers, dashpot and incident wave), R1 §2.7a.

---

## 6. Free-field motion and control-point normalisation (SITE Mode 2)

**Code:** `sassi/core/freefield.py`, `sassi/core/tlm.py` (`vertical_response`, `choose_mode`),
`sassi/modules/site.py`. **Verified by:** VP-02b, VP-15, VP-16, VP-24, VP-S2, example 4.

The seismic environment is a superposition of plane **wave fields** w: body waves (P, SV, SH) at
vertical or inclined incidence, and surface waves (Rayleigh, Love). For each wave and frequency SITE
computes the field at every user interface in the SITE axes x'y'z' (z' up, x' in the vertical plane
of propagation, y' = z' × x').

### 6.1 Vertically propagating body waves

With k = 0 the layered column decouples into 1-D columns: SV (x') and SH (y') use the shear matrices
(`G_s`, built with G*), P (z') the compression matrices (`G_p`, built with M* = λ* + 2G*). For an
incident (up-going) wave of displacement amplitude `E_b` in the half-space, the traction of the
half-space on the column base is `t = iωρV*(u_I - u_R) = iωρV*(2u_I - u)`, which gives (R1 §2.7a)

```
(G_• - ω²M + iωc e_N e_Nᵀ) u = 2iωc E_b e_N            (dashpot base, E_b = 1)          (6.1)
u_N = 1                                              (rigid base)
```

### 6.2 Normalisation to the control point

The control motion is specified at the **top of the control layer** `<cl>` in the direction `<cm>`
(x', y' or z'), as the **within** (in-column) motion (D-SIT-07). Each wave field is divided by its own
motion there:

```
Û_w(z) = u_w(z) / u_w,cm(z_cp)                                                       (6.2)
```

so the control point moves with unit amplitude; this ratio of displacements is also the ratio of
accelerations. With several waves the total field is

```
U(z, x') = Σ_w r_w(f) Û_w(z) exp(-i k_w (x' - x'_cp))                                 (6.3)
```

with participation ratios `r_w` given at two frequency numbers (Frequency 1 and 2), interpolated
linearly in between, held constant outside, each in (0, 1] and summing to 1 (Errors 53-55). Because
each field is unit-normalised and the ratios sum to one, the control point has unit motion. When the
design motion is an *outcrop* motion, the user runs SOIL first and uses the computed within motion
(R1 §2.7a; examples 4 and 6).

### 6.3 Inclined body waves and surface waves

* **Inclined P, SV, SH** (D-SIT-08): apparent wavenumber `k = ω sin θ / V_hs` (V = Vs for SV/SH, Vp
  for P; θ from the vertical). The generated buffer is not used: the user column sits on the
  **exact** half-space stiffness at wavenumber k (from P and SV potentials, Kausel and Roesset 1981),
  `K_out` for the radiating down-going field and `K_in` for the incident up-going field, and
  ```
  (A k² + B k + G - ω²M + K_out e_N e_Nᵀ) U = (K_out - K_in) u_I                        (6.4)
  ```
  with `u_I` the displacement of the pure incident wave at the top of the half-space. For θ = 0,
  `K_out = -K_in = iωρV*` and Eq. (6.4) reduces to Eq. (6.1).
* **Surface waves**: one mode of FILE2, chosen by WAVE `<opt>`: 1 "shortest wavelength" = the
  propagating mode with the largest Re k, 2 "least decay" = the mode with the smallest |Im k|
  (ties within 1e-8 max|k| broken by the largest Re k). A mode is *propagating* when
  `|arg k| ≤ atan(0.5) + atan(loss)`, the elastic sector widened by the material loss angle
  `tan(δ/2)` of the column (D-SIT-09, D-W1-08). The field is `φ_j(z) exp(-ik_j(x' - x'_cp))` with the
  physical vertical component `u_z = -iφ_z`, normalised by Eq. (6.2).

The manual's definitions of "shortest wavelength" and "least decay" and the treatment of the
half-space for inclined waves are reconstructions (R1 §10 items 1 and 3).

### 6.4 Evaluation at the interaction nodes

FILE1 holds, per wave and frequency, the unit-normalised field at the user interfaces, the
wavenumber and the ratio. ANALYS evaluates it at an interaction node on interface m at (x, y):

```
x' = (x - x_c) cos a + (y - y_c) sin a                                               (6.5)
U'(node) = R_z(a) Σ_w r_w U_w(m) exp(-i k_w x')
```

where a is the coordinate transformation angle (ANALYS `<ang>`, from x' to the global x axis,
counter-clockwise about z) and (x_c, y_c) the control point, which has x' = 0. Vertical incidence has
k = 0: every node of an interface moves in phase. VP-24 checks the rotation on a rotationally
symmetric model (`ATF_x = cos a · ATF(a=0)`, `ATF_y = sin a · ATF(a=0)`).

---

## 7. Point-load solutions: central zone and transmitting boundary (POINT)

**Code:** `sassi/core/axisym.py` (POINT3), `sassi/core/strip2d.py` (POINT2), `sassi/modules/point.py`,
`sassi/core/greens_tlm.py` (exact reference). **Verified by:** VP-07 (Boussinesq/Cerruti statics and
reciprocity), VP-08 (Wong 1975 dynamic surface Green functions), VP-09 (far field against the exact
thin-layer series), VP-43, and the impedance problems VP-10 to VP-14.

### 7.1 What POINT computes

"Each column of the flexibility matrix is formed by applying a unit point load at the interaction
degree of freedom associated with that column and by computing the resulting displacements at all
the interaction nodes" (Tajirian and Tabatabaie 1985). Because the site is horizontally layered, the
response to a load at depth n depends only on the horizontal distance to the observation point:
POINT solves one load per **load interface** n = 1 ... L+1 (L = POINT `<layer>`, the number of
embedment layers) on the axis of an axisymmetric model, and ANALYS translates and rotates the
solution to every pair of interaction nodes (section 8).

The load is expanded in Fourier harmonics about the vertical axis: a vertical load excites the
harmonic μ = 0, a horizontal load μ = 1, with the displacement pattern (symmetric case, R1 §3.2)

```
u_ρ = ũ_ρ cos μθ,     u_θ = -ũ_θ sin μθ,     u_z = ũ_z cos μθ                         (7.1)
```

A true point load gives an infinite displacement under it. SASSI therefore spreads it over a small
**central zone** of radius R0 (POINT `<rad>`), a finite-element core whose size calibrates the
diagonal (self) terms of the flexibility to the mesh of interaction nodes (section 7.5).

### 7.2 The core

One radial axisymmetric finite element (nodes on the axis and at ρ = R0) per SITE sublayer,
including the generated half-space sublayers and the base dashpots (decision D-PNT-03; the
discretisation of the original "special cylindrical axisymmetric elements" is a reconstruction,
R1 §3.3). Bilinear interpolation `L0 = 1 - ρ/R0`, `L1 = ρ/R0` in ρ and the thin-layer linear
interpolation in z; strains

```
ε = [ ∂ũ_ρ/∂ρ,  (ũ_ρ - μũ_θ)/ρ,  ∂ũ_z/∂z,  ∂ũ_ρ/∂z + ∂ũ_z/∂ρ,
      μũ_ρ/ρ + ∂ũ_θ/∂ρ - ũ_θ/ρ,  ∂ũ_θ/∂z + μũ_z/ρ ]ᵀ                                  (7.2)
```

with the isotropic `D(λ*, G*)`, `K_c = ∬ BᵀDB ρ dρ dz`, 3 × 3 Gauss points (ρ = 0 is never evaluated),
the mixed mass in z (as in SITE) and the consistent mass in ρ. Axis constraints: μ = 0 → ũ_ρ = ũ_θ = 0;
μ = 1 → ũ_ρ = ũ_θ (one tied unknown), ũ_z = 0. A viscous base adds `iωc ∫ L_a L_b ρ dρ` on the bottom
interface.

### 7.3 The exterior and the transmitting boundary

For ρ ≥ R0 the field is a superposition of outgoing Rayleigh and Love modes of FILE2:

```
ũ(ρ) = Ψ_μ(ρ) α,     α = [a_1 ... a_2Nf ; b_1 ... b_Nf]

         | Φ_x diag H'_μ(k^R ρ)           Φ_y diag μH_μ(k^L ρ)/(k^L ρ) |
Ψ_μ(ρ) = | Φ_x diag μH_μ(k^R ρ)/(k^R ρ)   Φ_y diag H'_μ(k^L ρ)         |               (7.3)
         | -Φ_z diag H_μ(k^R ρ)           0                            |
```

with `H_μ = H(2)μ` and `H'_μ = dH_μ(x)/dx` (the columns of Kausel's C_μ operator with J → H(2)). The
consistent generalised forces that the exterior exerts on the core at ρ = R0 (per radian) are

```
f_ρ = R0 [ E_(λ+2G) ∂ũ_ρ/∂ρ + E_λ (ũ_ρ - μũ_θ)/R0 + Q_λ ũ_z ]
f_θ = R0 [ E_G (μũ_ρ/R0 + ∂ũ_θ/∂ρ - ũ_θ/R0) ]                                         (7.4)
f_z = R0 [ Q_G ũ_ρ + E_G ∂ũ_z/∂ρ ]
```

with the thin-layer boundary integrals `E_c = ∫ c N Nᵀ dz` (per sublayer `c h/6 [2 1; 1 2]`) and
`Q_c = ∫ c N N'ᵀ dz` (per sublayer `c/2 [1 -1; 1 -1]`). Collecting them per mode column as `T_μ`, the
**transmitting boundary** (consistent boundary of Waas, Lysmer and Kausel) is

```
R_μ = -T_μ Ψ_μ(R0)⁻¹          (a full 3N_f × 3N_f matrix)                             (7.5)
```

Numerically, the Hankel functions are evaluated with the exponentially scaled `hankel2e` and the
basis is scaled to R0 (column j multiplied by `exp(ik_j R0)`), which avoids overflow for ρ ≥ R0 and
underflow of evanescent modes; R_μ does not depend on that scaling.

### 7.4 Solution

For each load interface n (unit load on the axis):

```
[ K_c - ω²M_c + R_μ ] Ũ = F̃,      F̃ = P/π (μ = 1, on the tied axis unknown),  F̃ = P/(2π) (μ = 0)   (7.6)
α = Ψ_μ(R0)⁻¹ Ũ_boundary,     ũ(ρ) = Ψ_μ(ρ) α   for ρ ≥ R0
```

(π and 2π are the angular integrals of cos², sin² and 1.) FILE3 stores, per frequency, load
interface and μ, the amplitudes α, the core-axis displacements at interfaces 1 ... L+1 and the mode
rows at the observation interfaces, so that ANALYS evaluates the field **exactly at any distance**
without tabulation (D-PNT-02).

**Exact reference.** For a point load in a thin-layer medium Kausel (1981, Eqs. 88-89) gives the
Green functions as Hankel series over the same modes, e.g. for a vertical unit load
`u_z = (1/4i) Σ_l φ_z^ml φ_z^nl H_0(k_l^R ρ)`. They are implemented in `sassi/core/greens_tlm.py` and
are the reference of VP-09: the POINT3 far field agrees within 2.35 % at 3.3 R0 (lead decision
D-W1-04), 0.8 % at 6.7 R0 and 0.25 % beyond 13 R0 on the dominant components (R1 V8). In the static
limit on a deep graded stratum the Green functions match Boussinesq and Cerruti within 0.1-1.7 %
(VP-07, R1 V3).

### 7.5 The central-zone radius

With one radial element, every pair of interaction nodes at a distance r ≥ R0 lies in the exact far
field; R0 only calibrates the near-singular diagonal terms. The rules of the manual are
R0 = 0.90 h for square meshes, 0.85 h for triangular (circular) meshes and R0 = h in 2D (h = mesh
size). The self-flexibility of the core is much softer than that of a uniform disk of radius R0
(equivalent disk radius 0.25-0.29 R0, R1 V9), so a "tributary disk" model would not reproduce
SASSI's diagonal terms. For non-uniform meshes `RADIUS` computes `r_e = Scale·√A_plan` per excavated
element and the average (D-PNT-05).

### 7.6 POINT2 (plane strain)

The 2D counterpart: a strip |x| ≤ R0 of two plane-strain elements per sublayer, with Waas-Lysmer
transmitting boundaries at x = ±R0 built from the same modes (outgoing mode j:
`u_x = φ_xj exp(-isk_j(x - x_b))`, `u_z = -isφ_zj exp(...)`, s = ±1), `R = -T Ψ⁻¹`; the far field is
`Σ α_j ψ_j exp(-ik_j(|x| - R0))`, the vertical response to a horizontal load being odd in x. The
reference is Kausel's line-load Green function (R1 V6, V7); VP-T1 checks the FILE3 far field against
an independent modal-sum implementation of it: within 0.33 % at |x| = 2R0 and 0.3 % farther out (lead
decision D-W3-12). HOUSE `<dim>` = 1 selects POINT2, and ANALYS solves 2D (plane-strain) models with
it (section 15).

---

## 8. Flexibility and impedance assembly

**Code:** `sassi/core/flexibility.py` (`flexibility_matrix`, `cylindrical_components`),
`sassi/core/ssi_solver.py` (`impedance_matrix`, `rigid_body_transform`, `global_impedance`).
**Verified by:** VP-07 (reciprocity), VP-10, VP-11, VP-13, VP-14, VP-17.

### 8.1 The 3 × 3 blocks

For interaction node i on interface m at (x_i, y_i) and node j on interface n, let
`r = |x_i - x_j|` and `θ = atan2(y_i - y_j, x_i - x_j)` (from the load point j to the observation point
i). From POINT (load at n, observed at m) take `ũ, ṽ, w̃ = (ũ_ρ, ũ_θ, ũ_z)` for μ = 1 and
`p̃, q̃ = (ũ_ρ, ũ_z)` for μ = 0. Rotating the cos θ / -sin θ pattern of Eq. (7.1) to Cartesian axes gives
(rows u_x, u_y, u_z at i; columns P_x, P_y, P_z at j; R1 §4.1):

```
        | ũc² + ṽs²    (ũ - ṽ)sc     p̃c |
F_ij =  | (ũ - ṽ)sc    ũs² + ṽc²     p̃s |      c = cos θ,  s = sin θ                    (8.1)
        | w̃c           w̃s            q̃  |
```

* r = 0 (same vertical line): `diag(ũ_axis, ũ_axis, q̃_axis)` from the core-axis values;
* 0 < r < R0 (mesh finer than the central zone): linear in r between the axis value and the value at
  R0, consistent with the core shape functions;
* r ≥ R0: the exact expansion `Ψ_μ(r) α`. When a model has more than 6,000 distinct distances, the
  expansion is tabulated on a dense grid (geometric near R0, spacing at most λ_min/48) and
  interpolated by cubic splines of `r·u(r)`.

**Reciprocity.** `F_ij = F_jiᵀ` holds exactly for the exact Green functions and to about 1e-4 with
the finite-element core (R1 V10, VP-07 checks `‖F - Fᵀ‖/‖F‖ < 1e-3`). SASSI-EDU symmetrises
`F ← ½(F + Fᵀ)` (D-ANL-02), because SASSI treats the impedance as symmetric.

**2D.** For POINT2 data the blocks are the 2 × 2 P-SV and 1 × 1 SH blocks with the odd coupling terms
multiplied by `sgn(x_i - x_j)`; ANALYS inverts the in-plane (UX, UZ) part only (section 15.1). For a
half or quarter model (SYMM) the blocks of the image nodes are added before the inversion, Eq. (15.2).

### 8.2 Impedance

```
X_ff(ω) = F_ff(ω)⁻¹          (manual Eq. 4.3: K + iωD = (f + ig)⁻¹)                     (8.2)
```

a full, complex symmetric `3n_f × 3n_f` matrix, computed by a dense LU factorisation and inversion in
place (LAPACK getrf/getri; never Cholesky). Its size dominates the memory of large models: 16 bytes
per entry, `(3n_f)² × 16` bytes per matrix (14.4 GB for 10,000 interaction nodes).

### 8.3 Global (unconstrained) foundation impedance

With ANALYS `<impe>` > 0 the impedance is condensed to the six rigid-body DOFs of the interaction
nodes about the control point (D-ANL-07):

```
K_G(ω) = Tᵀ X_ff T,   u_j = u0 + θ × d_j,   d_j = r_j - (x_c, y_c, z_c)

       | 1  0  0    0    dz  -dy |
T_j =  | 0  1  0   -dz   0    dx |                                                    (8.3)
       | 0  0  1    dy  -dx   0  |
```

and written as FOUNSTIF (Re K_G), FOUNDASH (Im K_G/ω), FOUNDAMP (Im K_G/(2|Re K_G|)) and FOUNIMPD
(|K_G|). This is the impedance of the soil at the interaction nodes when they are forced to move as a
rigid body; for a surface foundation it is the rigid-foundation impedance. VP-10, VP-11, VP-13 and
VP-14 compare it with the classical solutions for rigid disks and squares (Veletsos and Verbic 1973,
Pais and Kausel 1988, Tajirian and Tabatabaie 1985) and with independent boundary-element solutions
(welded and relaxed contact, [impedance_study.md](../verification/impedance_study.md)).

---

## 9. Solution of the SSI equation (ANALYS)

**Code:** `sassi/core/ssi_solver.py`, `sassi/modules/analys.py`, `sassi/modules/force.py`.
**Verified by:** VP-01, VP-15, VP-16, VP-17, VP-22 (restarts), VP-23 (simultaneous cases), VP-24,
VP-40 (moving load), VP-41 (low-frequency check), VP-E1.

### 9.1 Per frequency

For every SSI frequency `f_q = n_q Δf`:

1. **Frequency survey**: the frequency numbers must exist in FILE1 (or FILE9) and in FILE3, otherwise
   the run stops before any solve.
2. **Flexibility** `F_ff` from FILE3 (section 8.1) and **impedance** `X_ff = F_ff⁻¹`.
3. **Dynamic matrix** `C = (K*_s - ω²M_s) - (K*_e - ω²M_e)` (sparse, HOUSE matrices) plus
   `A_fᵀ X_ff A_f`.
4. **Load**: seismic `b_f = X_ff U'_f` with the free field of Eq. (6.5) at each interaction node
   (incoherent motion, wave passage and multiple excitation multiply the load or the motion by the HOUSE
   factors of FILE77, Eq. 14.7); vibration `b` = the FORCE vector, `P_k(ω) = a_k exp(-iωt_k)` for a load
   factor a_k and arrival time t_k on the reference load history (a delay is a phase lag; moving loads
   are loads with increasing arrival times, VP-40).
5. **Solve** (section 9.2) and write the transfer function of every equation to FILE8 (fixed DOFs:
   H = 0).

### 9.2 Schur complement

`X_ff` is dense on the interaction DOFs f, everything else is sparse. With the remaining DOFs n:

```
| C_nn   C_nf        | | U_n |   | b_n |
| C_fn   C_ff + X_ff | | U_f | = | b_f |

S   = C_ff - C_fn C_nn⁻¹ C_nf                    (dense Schur complement)
U_f = (S + X_ff)⁻¹ (b_f - C_fn C_nn⁻¹ b_n)                                            (9.1)
U_n = C_nn⁻¹ (b_n - C_nf U_f)
```

`C_nn` is factorised by a sparse LU (SuperLU) and `S + X_ff` by a dense LU, both complex symmetric
but solved with LU (never Cholesky or a Hermitian solver). Every right-hand side (simultaneous cases,
restarts) reuses both factorisations. If `C_nn` alone is singular or ill-conditioned (estimated
reciprocal condition number below 1e-12, for example an undamped structure with clamped interaction
DOFs at one of its fixed-base frequencies) the full system with the dense block inserted is
factorised by a sparse LU instead, with a warning. The cost per frequency is dominated by the
O(n_f³) dense operations.

### 9.3 Restarts

| ANALYS `<mode>` | Name | Stored by the initiation run | Re-used |
|---|---|---|---|
| 1 | New Structure | `COOXqqq` = X_ff per frequency | X_ff; the structure is re-assembled and re-factorised |
| 2 | New Seismic Environment | `COOTKqqq` = the factorised system (sparse LU of C_nn, Schur data, dense LU of S + X_ff) | X_ff and the factorisation; only b changes |
| 3 | New Dynamic Loading | same as 2 | same as 2 |

Each record carries the frequency number and value and the FILE90 hashes of the interaction nodes,
the layering and (COOTK) the structure matrices; a restart refuses records that do not match, and
`COOXI`/`COOTKI` index the records by frequency number so that any subset of the saved frequencies
can be re-run. A restart gives the same FILE8 as a fresh initiation run to round-off (VP-22,
tolerance 1e-10).

### 9.4 Simultaneous cases and the low-frequency check

With `<simul>` = 1 the three coherent input directions (FILE1X/Y/Z, angle 0) are three right-hand
sides of one factorisation and give FILE8X/Y/Z; with `<simul>` = Nl ≥ 2 in a vibration analysis the
Nl load cases FILE9001 ... give FILE8001 ... VP-23 checks that each case equals a single run to
1e-12. An incoherent stochastic simulation with Ns samples (`<simul>` = Ns ≤ 50) solves the 3 Ns
right-hand sides FILE1X/Y/Z × FILE77001 ... FILE77Ns and writes `FILE8{3(s-1)+d}` (sample s,
direction d = 1, 2, 3; section 14.4). In 2D, `<simul>` = 1 analyses the in-plane X and Z cases only.

As the frequency tends to zero the whole system moves rigidly with the free field, so every
transfer function in the input direction tends to 1. ANALYS reports the deviation at the first SSI
frequency (G-19); VP-41 checks it on several models.

---

## 10. Transfer-function interpolation

**Code:** `sassi/core/interp.py` (`interpolate_tf`, `smooth_tf`, `phase_adjust`), used by MOTION and
STRESS. **Verified by:** VP-28 (exactness, ISRS from interpolated TFs), VP-32 (SRSS TF), VP-33 (phase
adjustment), VP-53 (CRITFREQ).

ANALYS solves at a few tens to a few hundred SSI frequencies; the convolution needs the transfer
function at every Fourier frequency `f_k = kΔf` up to the cut-off.

### 10.1 The SASSI (Tajirian) interpolant

Over a short band, every SSI transfer function behaves like that of a **two-degree-of-freedom system
with hysteretic damping**, whose transfer function is a ratio of quadratics in ω² (no odd powers,
because hysteretic damping does not depend on frequency; Tajirian 1981, TLUSH 1981):

```
H(ω) = (C1 ω⁴ + C2 ω² + C3) / (ω⁴ + C4 ω² + C5)                                        (10.1)
```

The five complex constants are fitted to five computed values `H_p = H(ω_p)` of a **window** of five
consecutive SSI frequencies:

```
[ ω_p⁴   ω_p²   1   -ω_p² H_p   -H_p ] · {C1 C2 C3 C4 C5}ᵀ = ω_p⁴ H_p,     p = 1..5          (10.2)
```

(the MHI transcription prints +H_p in the last column; moving the denominator of Eq. 10.1 across
gives -H_p, and only that sign reproduces an exact 2-DOF transfer function, R1 §5.1). SASSI-EDU uses
the scaled variable `x = (ω/ω_max)²` of the window for conditioning. When the window holds data of
fewer than two modes (an SDOF-like or flat transfer function) the 5 × 5 system has rank 4: every
member of the solution family interpolates, and the minimum-norm least-squares (SVD, relative
singular-value cut 1e-10) solution is used. The fit reproduces an exact 2-DOF hysteretic transfer
function to about 1e-14 (R1 V5; VP-28 checks 1e-10).

**Guard against spurious poles** (D-MOT-01 as amended by D-W1-10). When the fitted denominator nearly
vanishes on the window span (`min|D| < 1e-3 max|D|`, `D(x) = x² + C4 x + C5`) the window *may* hold a
spurious near-real pole. The pole is accepted when it is cancelled by a numerator zero, when it is a
physical damped mode (for `exp(+iωt)` a hysteretic mode has its pole at `x_p = (ω0/ω_max)² c(β)`,
i.e. `0 < arg x_p < π` with implied damping `sin(arg x_p/2) ≥ 1e-4`), or when it is remote from the
span; otherwise the complex cubic through the four nearest window points replaces the rational form
in that window. (The literal ratio test alone would destroy exact fits on lightly damped poles.)

### 10.2 Window schemes (options 0-6)

For the interval [f_j, f_j+1] the candidate windows are those of five consecutive points that contain
it: start index m with `max(1, j-3) ≤ m ≤ min(j, N-4)` (1-based).

| Option | Windows | Combination |
|---|---|---|
| 0 (SASSI2000) | all candidate windows | weighted mean, `w_m = max(1e-3, 1 - \|f - c_m\|/h_m)` (c_m, h_m centre and half-width of window m) |
| 1 (SASSI 1982) | starts 1, 5, 9, ... sharing end points; the last window is the last five points | single window |
| 2 | all candidate windows | arithmetic mean |
| 3 | the three candidates with the smallest `\|f - c_m\|` | arithmetic mean |
| 4 | starts 2, 6, 10, ...; the first interval uses m = 1 | single window |
| 5 | starts 3, 7, 11, ...; the first two intervals use m = 1 | single window |
| 6 | independent not-a-knot cubic splines of Re H and Im H through the computed points and the f = 0 anchor | — |

With fewer than five SSI frequencies options 0-5 become option 6, with fewer than three complex linear
interpolation. Every scheme returns the computed values exactly at the SSI frequencies (TFI = TFU,
enforced bit-exactly). The weighting of option 0 and the shifts of options 4 and 5 are
reconstructions (R1 §10 item 4).

**Outside the computed range** (D-MOT-03): above the last SSI frequency f_N the transfer function is
zero (the cut-off is a low-pass filter); below the first one a seismic transfer function is linear
between the rigid-body anchor `H(0)` and `H_1` (H(0) = projection of the control direction on the
output DOF: 1 parallel to the input, cos a / sin a for rotated horizontals, 0 otherwise and for
rotations), a vibration transfer function is held at `H_1`.

### 10.3 Smoothing and phase adjustment

**Smoothing** (MOTION `<smo>` = S, options 0-5): in [f_j, f_j+1], with the complex linear reference
L(f) and the band `[A_min, A_max]` of the two neighbouring amplitudes,

```
r = max(0, |H|/A_max - 1) + max(0, 1 - |H|/max(A_min, 1e-12 A_max))
H_s = L + (H - L)/(1 + S r)                                                            (10.3)
```

Values inside the band are unchanged; S = 0 is the identity. Genuine resonance peaks between
computed points are clipped too, hence S = 0 for coherent analyses.

**Phase adjustment** (MOTION `<pzadj>` = 1; decision D-MOT-05 as amended by D-W1-11):

```
φ = unwrap(arg H) along the grid from f = 0,   φ0 = phase at f = 0 (or the first non-zero value)
H' = |H| exp(i(φ0 + ρ(φ - φ0))),   ρ = 1/(1 + S) (options 0-5),  ρ = 0 (option 6)            (10.4)
```

Referring the phase to φ0 keeps negative rigid-body anchors (a reversed input direction reverses the
response). With ρ = 0 the transfer function has zero differential phase, the "upper bound" approach of
incoherent analyses; VP-33 checks both cases exactly.

**SRSS transfer function** (MOTIONX `<srss>` with `SRSSTF.txt`, D-MOT-06): the modal FILE8 transfer
functions are interpolated and combined as `|H| = √(Σ|H_k|²)` with zero phase, or with the phase of
the coherent solution (VP-32).

---

## 11. Convolution, response spectra, relative displacements and stresses

**Code:** `sassi/core/signal.py`, `sassi/core/spectra.py`, `sassi/modules/motion.py`,
`sassi/modules/reldisp.py`, `sassi/core/stress_lib.py`, `sassi/modules/stress.py`.
**Verified by:** VP-29 (identity convolution), VP-30 (response spectra), VP-31 (baseline correction),
VP-34 (relative displacement), VP-35, VP-S1, VP-S2 (stress recovery), VP-40.

### 11.1 Convolution

The control acceleration a(t) (in g, records rec1 ... rec2, scaled by `<mult>` or to `<max>`) is
zero-padded to NFFT points and transformed, `A = rfft(a)`. The frequency step must equal the FILE8
step (relative 1e-6). Then

```
seismic acceleration:      a_k(t) = irfft( H_k(f) A(f) )                (in g)
quantities linear in displacement (relative displacements, strains, stresses):
                           U_g(f) = -g A(f)/ω²   (f = 0 term set to 0),   r(t) = irfft( H_r U_g )     (11.1)
vibration:                 u = irfft(H F),  v = irfft(iωH F),  a = irfft(-ω²H F)
```

with F the FFT of the reference load history. Because FILE8 holds the dimensionless ratio
`U/U_cp`, the same transfer function serves for accelerations and displacements (D-CNV-06). The
signal is periodic with period `NFFT·Δt`; the zero-padded tail (quiet zone) must be long enough for
the free vibration to decay. The `.ACC` history is written for `min(NFFT, 1.2·dur/Δt)` samples.
VP-29 checks that H ≡ 1 returns the input to 1e-12.

### 11.2 Baseline correction (MOTION `<bl>` = 1)

Hudson-Housner time-domain correction (D-MOT-08): integrate the acceleration (trapezoidal rule) to
the velocity v; fit `v ≈ c1 t + c2 t²/2` by least squares (the velocity produced by an acceleration
baseline error `c1 + c2 t`); subtract `c1 + c2 t` from the acceleration and integrate again. VP-31
removes a constant acceleration offset to 1e-6.

### 11.3 Response spectra (Nigam-Jennings)

For an oscillator of frequency ω and damping ζ, `ẍ + 2ζωẋ + ω²x = -a(t)`, with the ground
acceleration linear between samples, the exact recurrence of Nigam and Jennings (1969) advances the
state `s = (x, ẋ)`:

```
s_(k+1) = A s_k + B0 a_k + B1 a_(k+1)                                                  (11.2)
```

with A, B0, B1 functions of ω, ζ, Δt (R2 G.1; `spectra._nj_coefficients`). The spectral acceleration
is the maximum absolute *absolute* acceleration `|2ζωẋ + ω²x|` (written to `.RS`); SV and SD are the
maxima of the relative velocity and displacement, `PSA = ω²SD`. Two refinements keep the peak
sampling error small: for oscillator frequencies above 0.1/Δt the record is FFT-upsampled four times
(D-MOT-07), and the piecewise-linear input is sub-divided so that every oscillator period has at
least 64 integration steps (D-W1-07; peak error ≤ 1 - cos(π/64) = 0.12 %). The spectra are computed
over the full NFFT record at `<fstep>` log-spaced frequencies from `<freq1>` to `<freq2>`. VP-30 checks
the closed forms: step input `PSA/a0 = 1 + exp(-ζπ/√(1-ζ²)) = 1.854468` at ζ = 5 %, harmonic input at
resonance `A/(2ζ)`, impulse `PSV/ΔV = 0.92669`.

### 11.4 Relative displacements (RELDISP)

From the complex `.TFI` of a node and of a reference (another node's `.TFI`, or the free field: unit
amplitude and zero phase in the input direction, D-RDP-04):

```
D(f) = (H_node(f) - H_ref(f)) · U_g(f),      d(t) = irfft(D)                              (11.3)
```

This avoids the drift of double integration. Above the cut-off the node TFs are zero while the
free-field reference is still 1, so a free-field-referenced relative displacement contains the small
ground displacement above f_N: exactly `irfft((H - 1)U_g)` (VP-34). The `.TFD` file is the complex
relative-displacement transfer function per unit control acceleration (length per g).

### 11.5 Stresses and forces (STRESS)

HOUSE stores per output element a complex **recovery operator** S, so that a component is `S_c · u_e`
(SOLID/PLANE centroid stresses `D*B u` in global axes, SHELL membrane stresses and moments per unit
length `D*_m B_m u`, `D*_b κ` in local axes, BEAMS local end forces `k_L T u` exerted on the element,
SPRING `k*(u_J - u_I)`). STRESS then

1. forms the **stress transfer functions** at the SSI frequencies, `STF_c(f_j) = S_c (U_e(f_j) - r_e)`,
   where `r_e` is the rigid-body part (the control motion projected on the element DOFs; S annihilates
   it in exact arithmetic, and subtracting it keeps the low-frequency values accurate);
2. **interpolates the STF** (not the nodal TFs) with the MOTION schemes; below the first SSI frequency
   a seismic STF goes linearly to 0 at f = 0 (a rigid-body motion has no stress);
3. **convolves** with the control displacement spectrum (Eq. 11.1);
4. computes the non-linear derived quantities in the time domain, e.g. the octahedral shear stress
   `τ_oct = ⅓√[(σxx-σyy)² + (σyy-σzz)² + (σzz-σxx)² + 6(τxy² + τyz² + τxz²)]`.

VP-35 checks that a rigid-body motion gives zero stress, `τ_oct = (√6/3)τ` for pure shear and the
equilibrium of beam end forces; VP-S1 the base moment of a cantilever stick against an independent
calculation; VP-S2 the shear stress of a soil column under vertical SV against `G* × free-field strain`.

---

## 12. The SHAKE equivalent-linear method (SOIL)

**Code:** `sassi/core/shake.py`, `sassi/modules/soil.py`. **Verified by:** VP-02a (closed-form layer
transfer functions, 1e-8), VP-04 (SHAKE91 sample problem), example 4.

### 12.1 Wave solution and recursion

In sublayer m (local depth z_m measured down from its top) a vertically propagating shear wave is the
sum of an up-going wave E_m and a down-going wave F_m (Schnabel, Lysmer and Seed 1972; SHAKE91):

```
u_m = E_m exp(i(ωt + k*_m z_m)) + F_m exp(i(ωt - k*_m z_m)),    k*_m = ω/V*_m,   V*_m = √(G*_m/ρ_m)   (12.1)
```

Continuity of displacement and shear stress `τ = G* ∂u/∂z` at each interface gives

```
E_(m+1) = ½ E_m (1 + α_m) exp(ik*_m h_m) + ½ F_m (1 - α_m) exp(-ik*_m h_m)
F_(m+1) = ½ E_m (1 - α_m) exp(ik*_m h_m) + ½ F_m (1 + α_m) exp(-ik*_m h_m)                (12.2)
α_m = ρ_m V*_m / (ρ_(m+1) V*_(m+1))
```

starting from the free-surface condition `E_1 = F_1 (= 1)`. The **within** motion at the top of
sublayer m is `E_m + F_m`, the **outcrop** motion `2E_m`. Transfer functions are ratios of these
quantities; the input is applied at the top of the control sublayer as outcrop or within motion, and
the last SPRO entry is the half-space (D-SOL-02). VP-02a checks the closed form
`1/[cos k*H + iα* sin k*H]` of a uniform layer to 1e-8 for both damping forms.

### 12.2 Strains and the equivalent-linear iteration

The shear strain at mid-height of each sublayer, from the acceleration-wave amplitudes, is

```
γ_m(ω) = ( E_m exp(ik*h/2) - F_m exp(-ik*h/2) ) / (iωV*_m)                               (12.3)
```

transformed to the time domain to obtain `γ_max = max|γ(t)|`. The **effective strain** is
`γ_eff = R_γ γ_max` (SOIL `<ratio>`, typically 0.5-0.7; SHAKE91 suggests (M - 1)/10 for magnitude M).
New properties are read from the strain-dependent curves:

```
G = G_max · (G/G_max)(γ_eff),      β = D(γ_eff)                                           (12.4)
```

interpolated **linearly in log10(γ)** with constant extrapolation (SHAKE91's rule; D-SOL-03). SOIL
runs exactly `<iter>` iterations (SHAKE91 has no convergence test; 8 are recommended) and lists the
changes of G and β per iteration; the half-space is not iterated. For vertical input (SOILX
`<indir>` = 1) the P-wave velocity and damping are used without iterations. Fourier components above
`EDUOPT,SOILCUTOFF` are removed (D-SOL-08; 25 Hz reproduces the SHAKE91 sample problem).

### 12.3 Hand-off to SSI

FILE88 holds per sublayer the thickness, the effective strain, G, Vs, βs, Vp and βp. With `SITEX,1`
SITE replaces the TOPL layer properties by these, position by position (the layer counts must match,
EDU-07). For Vp and βp the default policy keeps Poisson's ratio (`EDUOPT,VPPOLICY,NU`; `VP` keeps Vp
for saturated soils) and sets βp = βs (D-SOL-06). The same complex-modulus form must be used in SOIL
and SITE, otherwise the SSI free field would not reproduce the SOIL motion at the control point.

**VP-04 and the SHAKE91 sample problem.** SOIL reproduces the published SHAKE91 results (strains,
moduli, damping, transfer functions, peak accelerations) within the planned tolerances. The printed
SHAKE91 *spectra* turned out to be spectra of a record corrupted by a defect of SHAKE91 subroutine
DRCTSP (samples that set a new running maximum are not converted from g); emulating that defect
reproduces all 151 printed values within 0.6 %, which is the VP-04 criterion, and the correct spectra
are reported in the notes (lead decision D-W1-06).

---

## 13. Element formulations (HOUSE)

**Code:** `sassi/elements/` (`solid.py`, `plane.py`, `beam.py`, `shell.py`, `tshell.py`, `spring.py`,
`general.py`, `base.py`, `assemble.py`), `sassi/core/renumber.py` (node optimizer),
`sassi/modules/house.py`. **Verified by:** VP-03, VP-37, VP-38, VP-38T, VP-39, VP-39b, VP-54T, VP-TS1,
VP-H1, VP-H2, VP-A1, VP-O1.

HOUSE assembles four frequency-independent sparse matrices on one global DOF map: `K*_s` (complex)
and `M_s` of the structure (including near-field soil) and `K*_e`, `M_e` of the excavated soil.
Because damping is hysteretic, none depends on the frequency; ANALYS combines them at each frequency.

### 13.1 The element library

| Element | Stiffness | Mass |
|---|---|---|
| SOLID | 8-node isoparametric hexahedron, `D* = λ* D_λ + G* D_G`; Gauss rule 2×2×2 / 3×3×3 / 4×4×4 for EINT 0/1/2; for structural solids 9 Wilson incompatible modes (1-ξ², 1-η², 1-ζ² per direction) with the Taylor centroid-Jacobian correction, condensed statically `K = K_uu - K_ua K_aa⁻¹ K_au` (passes the patch test); prisms and pyramids by repeated nodes; positive Jacobian required (EDU-05) | ½ lumped (row sum) + ½ consistent |
| PLANE | 4-node plane strain in the X-Z plane, unit thickness, 2×2 Gauss, 4 incompatible modes for structural elements; either node orientation accepted | ½ + ½ |
| BEAMS | 2-node 3D Timoshenko frame, `φ = 12EI/(G A_s L²)` (A_s = 0: Euler-Bernoulli); bending in the 1-2 plane with I3 and As2, in the 1-3 plane with I2 and As3; end releases by static condensation of the local matrices; exact for end loads | consistent (translations from the Timoshenko interpolation, axial `ρAL/6[2 1; 1 2]`, torsion `ρ(I2+I3)L/6[2 1; 1 2]`, no rotary inertia) |
| SHELL | flat facet: membrane = plane-stress Q4 with 4 incompatible modes (CST for triangles); plate = DKQ / DKT discrete-Kirchhoff elements (Batoz et al. 1980, 1982); **zero drilling stiffness**; warped quads projected on the mean plane | lumped `ρtA/n` on translations, no rotary inertia |
| TSHELL | flat facet, Mindlin-Reissner plate (section 13.2): membrane as SHELL; transverse shear by the MITC4 assumed natural strains (Bathe and Dvorkin 1985) for both EINT values; bending EINT 1 = 2×2 Gauss, EINT 0 = one point with hourglass stabilisation (coefficient 0.1); triangles: CST membrane, MITC3 shear with a condensed rotation bubble and a stabilised curl part; shear factor 5/6; drilling penalty 1e-4 × the smallest membrane diagonal × area | lumped `ρtA/n` on translations and the Mindlin rotary inertia `ρt³A/(12n)` on the two bending rotations; none on the drilling rotation |
| SPRING | six uncoupled global springs `k*_d = k_d c(damp)` between I and J; zero length allowed | none |
| GENERAL | user 12×12 `K_R + iK_I` and mass (÷g in weight units); with a K node, local axes as BEAMS and `T = blockdiag(Λ, Λ, Λ, Λ)` | user |
| nodal masses | MT/MR on the diagonal (÷g with MUNITS = 1), ignored on fixed DOFs | |

**Beam stiffness in one bending plane** (deflection v, rotation θ = dv/dx; Przemieniecki 1968):

```
           E I        |  12      6L         -12     6L        |
k = ---------------   |  6L   (4+φ)L²       -6L   (2-φ)L²     |                             (13.1)
     (1 + φ) L³       | -12     -6L          12    -6L        |
                      |  6L   (2-φ)L²       -6L   (4+φ)L²     |
```

with the sign of the coupling terms flipped in the 1-3 plane (θ2 = -du3/dx). Local axes: e1 = I → J,
e2 toward K in the plane IJK, e3 = e1 × e2; BEAMS end forces are the forces exerted on the element.

**Excavated soil.** Excavated SOLID/PLANE elements take the L-table layer of their MSET index
(D-ELM-12): `G = ρVs²`, `M = ρVp²`, ρ = weight/g, with the complex moduli of Eq. (2.3); never
incompatible modes. ETYPE 0 elements are excavated when all nodes are at or below the ground
elevation and the centroid is strictly below it (D-HOU-01). VP-H2 checks the excavated mass `Σ ρV` and
the stiffness identity with the layer properties.

**DOF management** (D-ELM-09). The active DOFs of a node are the union of the DOFs of its elements
minus the fixed ones; DOFs that no element defines (rotations of solid-only nodes) are removed
automatically, DOFs that are defined but not stiffened (shell drilling rotations) are kept and need
FIXROT/FIXSHLROT. Interaction DOFs are the translations of the interaction nodes.

**What the element VPs show.** VP-37: Euler-Bernoulli and Timoshenko cantilevers, end releases,
simply supported plates and bars, patch tests (exact constant stress). VP-38: NAFEMS free-vibration
benchmarks FV12, FV16 (SHELL bending; the lumped translational mass converges from below at O(h²),
hence the 32 × 32 mesh of lead decision D-W1-05) and FV32 (membrane). VP-03: the dispersion of
discrete soil columns with consistent, lumped and mixed mass (the mixed mass gives ω_h/ω = 0.99451 at
five elements per wavelength). The TSHELL problems are described in section 13.2.

### 13.2 The TSHELL thick shell (Mindlin-Reissner)

The thin SHELL element (Kirchhoff) neglects transverse shear deformation. For thick basemats, the
shear walls of nuclear buildings or the higher modes of any plate the Mindlin-Reissner theory is the
better model: the normal to the mid-surface stays straight but not normal, so the rotations β of the
normal are independent of the slope of w and the transverse shear strain carries the shear forces.
In the local axes x'y'z' of the facet (defined as for SHELL), with right-hand rotations θ':

```
u(z) = u' + z β_x,   v(z) = v' + z β_y,   β_x = θ_y',   β_y = -θ_x'

membrane   ε = [u,x   v,y   u,y + v,x]             N = D_m ε,   D_m = E t/(1-ν²) [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
bending    κ = [β_x,x   β_y,y   β_x,y + β_y,x]     M = D_b κ,   D_b = E t³/(12(1-ν²)) [same pattern]
shear      γ = [w,x + β_x   w,y + β_y]             Q = D_s γ,   D_s = κ_s G t I,   κ_s = 5/6       (13.2)
```

so M_xx is positive with tension on the +z' face and the plate equilibrium reads
`Q_x = M_xx,x + M_xy,y`. In the Kirchhoff limit γ → 0 gives β = -∇w, the SHELL convention.

**Shear locking and the MITC4 field.** With bilinear w and β a thin plate cannot make γ = 0
everywhere and the element becomes far too stiff (with plain 2×2 integration of the shear, a cantilever
plate 1000 times longer than thick, 8 × 2 elements, deflects about 6500 times too little; VP-TS1). The
quadrilateral therefore uses the assumed natural strains of Bathe and Dvorkin (1985): the covariant
shear strains `γ_ξ = w,ξ + β·x,ξ` and `γ_η = w,η + β·x,η` are sampled at the edge mid-points
A (0, 1), C (0, -1), B (1, 0), D (-1, 0) and interpolated,

```
γ̃_ξ = ½(1 + η) γ_ξ(A) + ½(1 - η) γ_ξ(C),     γ̃_η = ½(1 + ξ) γ_η(B) + ½(1 - ξ) γ_η(D),     γ = J⁻¹ [γ̃_ξ, γ̃_η]    (13.3)
```

integrated with 2×2 Gauss points for both EINT values. At the centre γ̃ is the one-point (tying
average) value of the manual's "1-point transverse shear", and the linear terms act as an
assumed-strain stabilisation of it: a literal one-point shear would leave the w = ±1 hourglass mode
without stiffness. The MITC4 field neither locks nor has a spurious mode.

**Bending.** EINT 1 (*selective*): 2×2 Gauss integration. EINT 0 (*reduced*, the manual's default):
one point at the centre plus a physical hourglass stabilisation,

```
K_b = A B₀ᵀ D_b B₀ + ε_hg Σ_g w_g |J_g| (B_g - B₀)ᵀ D_b (B_g - B₀),     ε_hg = 0.1          (13.4)
```

`B_g - B₀` vanishes for every linear rotation field, so constant curvature is exact on any
quadrilateral (patch test) and only the two rotation hourglass modes get ε_hg times their 2×2
stiffness. **Membrane** (both EINT values): the SHELL membrane, a Q4 with the four Wilson-Taylor
incompatible modes condensed.

**Triangles.** CST membrane. Plate: w linear, rotations linear plus a cubic bubble
`f4 = 27 r s (1 - r - s)` whose two DOFs are condensed inside the element; transverse shear from the
MITC3 field of the corner DOFs (Lee and Bathe 2004), `γ̃ = a + c(-(y - y_c), x - x_c)`, whose
tangential component on each side equals the side average of the displacement-based strain; the bubble
enters through its element mean, `γ₀ = a + (27/60) β_b`. The shear energy is
`A γ₀ᵀ D_s γ₀ + w_lin ∫ (c rot)ᵀ D_s (c rot) dA` with the Lyly-Stenberg-Vihinen (1993) weight
`w_lin = t²/(t² + α h²)`, α = 0.2 and h the longest side, on the linear ("curl") part only. The plain
MITC3 triangle locks: a simply supported plate with t/a = 1/1000 on a 16 × 16 mesh deflects 46 % of
Navier's value.

**Drilling rotation.** "A small rotational stiffness is automatically added in HOUSE" (manual): the
nodal drilling rotations are tied to the in-plane rotation of the membrane at the centre,
`ω = (v,x - u,y)/2`, by the penalty `k_d Σ_i (θ_z',i - ω)²` with `k_d = 10⁻⁴ min(diag K_m) A`. It
vanishes for a rigid in-plane rotation, so the free element keeps its six rigid-body modes, and it is
too small to alter the membrane response. TSHELL models therefore need no FIXROT.

**Mass, damping, recovery.** Lumped `ρtA/n` on each translation and `ρt³A/(12n)` on the two bending
rotations (required to reproduce the Mindlin frequencies of NAFEMS Test 21, which include rotary
inertia); no drilling inertia (an inertia on a penalty DOF would create spurious low-frequency modes).
`E* = E c(β_s)` with ν real, so `K* = c(β_s) K₀`. STRESS recovers at the centre the membrane forces
NXX NYY NXY, the transverse shear forces QXZ QYZ (MITC field) and the moments MXX MYY MXY per unit
length. With `THSHLSTR,1` it adds the face stresses of D-TSH-01 from the maxima of the eight
components, which occur at different times with unknown signs:

```
SXX = s_N NXX/t + s_M 6 MXX/t²,   SYY = s_N NYY/t + s_M 6 MYY/t²,   (s_N, s_M) = ++, --, +-, -+
TXY = NXY/t + 6 MXY/t² (++ only),   TXZ = 1.5 QXZ/t,   TYZ = 1.5 QYZ/t (mid-surface maxima)      (13.5)
```

with the principal stresses and the plane-stress strains of every permutation.

**Verification.** VP-TS1: membrane, bending and transverse-shear patch tests exact (1e-10) on
distorted quadrilaterals and triangles; no locking (cantilever plate t/L = 1/1000 within 1 % of
Kirchhoff on 8 × 2 elements and 0.02 % on 32 × 4; simply supported plates t/a = 1/1000 and 1/10 with
both triangulation patterns within 2 % (16 × 16) and 0.5 % (32 × 32) of Navier, observed order 2);
exactly six zero-energy modes; `K* = c(β) K₀`. VP-38T: NAFEMS Test 21 (simply supported thick plate
10 × 10 × 1 m, 45.897, 109.44 (×2), 167.89 Hz ...) within 2 % on 16 × 16 and 32 × 32 meshes with an
observed order of 2; the NAFEMS 8 × 8 rows are informative (lead decision D-W3-01: the lumped mass
converges from below at O(h²), -5.4 % on mode 4 with 8 × 8); the thin limit FV12/FV16 within 0.5 % of
the Kirchhoff SHELL on the same mesh. VP-54T: the THSHLSTR face stresses N/t and ±6M/t² of pure
membrane and pure bending states.

### 13.3 The node-numbering optimizer (HOUSEX `<optimize>` = 1)

A profile (skyline) solver, the solver of the original SASSI, stores every row of the matrix from its
first non-zero to the diagonal, so its cost grows with the profile `Σ_i (i - first_i)` and the bandwidth.
The optimizer (D-HOU-03) renumbers the nodes by reverse Cuthill-McKee (Cuthill and McKee 1969; George
and Liu 1981): the breadth-first search starts from **all interaction nodes at once**, so they get the
last numbers, ascending in their bottom-up table order (EDU-21 by construction), and the dense `X_ff`
block that ANALYS adds stays in the last rows of the profile. Components without interaction nodes are
numbered first by classical RCM from a pseudo-peripheral node; nodes without element DOFs come last.
The new order is accepted only if neither the bandwidth nor the profile grows (spec 05b test 8).
HOUSE writes `<model>.hounew` and the `old new` pairs `<model>.map`; FILE4 and FILE8 are in the new
numbering. SASSI-EDU's sparse solver chooses its own fill-reducing ordering, so the results do not
change: VP-O1 maps the FILE8 of an optimised, deliberately badly numbered TSHELL basement back through
the `.map` and finds the reference transfer functions to 2e-13 (tolerance 1e-10), with the equation
profile reduced to 59 % of the original.

---

## 14. Incoherency, wave passage and multiple excitation

**Code:** `sassi/core/coherency.py` (coherency models), `sassi/core/incoherency.py` (decomposition,
synthesis, wave passage, multiple excitation, FILE77), `sassi/modules/house.py` (writes FILE77),
`sassi/core/ssi_solver.py` (`incoherent_seismic_load`), `sassi/modules/analys.py`,
`sassi/modules/motion.py` (SRSS TF), `sassi/prep/commands/incoherency_cmds.py` (INCOH, WPASS, ME, AMP,
HOUSEX, ANALYSX, BUILDFILE77). **Data:** `sassi/data/coherency/*.json`. **Verified by:** VP-26, VP-27,
VP-I1, VP-I2 (requirements §4.4 item 6, §4.6 item 5; decisions D-INC-01 ... D-INC-12, D-W3-05,
D-W3-06).

### 14.1 What is modelled

Two points of the free-field surface a horizontal distance D apart do not move alike. Part of the
difference is a **delay** (the wave sweeps across the foundation at an apparent velocity: *wave
passage*), the rest is random (*incoherence*, from scattering). The **unlagged (plane-wave) coherency**
γ(f, D) is the correlation coefficient of the Fourier amplitudes of the two motions at frequency f
after the delay has been removed: 1 at D = 0 or f = 0, decaying with frequency and distance. SASSI
applies the two effects separately: HOUSE builds complex **incoherency factors** `s_i(ω)` per SSI
frequency, motion component (X, Y, Z) and interaction node, from the coherency matrix of the
interaction nodes and from the wave-passage delays, and writes them to FILE77; ANALYS multiplies the
coherent free field of SITE by them:

```
u_i,inc(ω) = s_i(ω) u_i,coh(ω)                                                             (14.1)
```

### 14.2 Coherency models

`WPASS <cohf>` selects the model (D-INC-01). The coefficients of the Abrahamson models are **data**,
transcribed from the primary reports into `sassi/data/coherency/*.json` with the citation (report,
equation, table, page) and points read from the published figures, which the unit tests reproduce;
they are never typed from memory (D-INC-04). A model without a verified file refuses to run with a
"coefficients not available" error that names the alternative.

| `<cohf>` | Model | γ(f, D) | Status and source |
|---|---|---|---|
| 1 | Luco and Wong (1986) | `exp[-(γ_c ω D / V_s)²]`, ω = 2πf in rad/s (D-INC-05); γ_c = INCOH `<gammax>`, `<gammay>`, `<gammaz>` per component (≥ 0.1, Error 58), V_s = INCOH `<alpha>`, the mean shear-wave velocity (Error 59); D = plain horizontal distance | available; form confirmed from the EPRI/NRC CLASSI incoherency presentation (NRC ADAMS ML072620217) |
| 2 | Abrahamson (1993), all soil types | - | refused: the seminar paper could not be retrieved |
| 3 | Abrahamson (2005), all sites, surface foundations | Eq. (14.2) with `a2 fc` in the second factor | available: EPRI 1012968, Eq. 3-1, Tables 3-4 and 3-5 (also EPRI 1015110 Eq. 5-2) |
| 4 | Abrahamson (2006), all sites, embedded foundations | - | refused: EPRI 1014101 could not be retrieved (EPRI 1015110 recommends the hard-rock model for embedded foundations) |
| 5 | Abrahamson (2007), hard rock | Eq. (14.2) | available: EPRI 1015110 chapter 6, Eq. 6-1, Tables 6-1 and 6-2 |
| 6 | Abrahamson (2007), soil sites, surface foundations | Eq. (14.2) | available: EPRI 1015110 chapter 7, Eq. 7-1, Tables 7-1 and 7-2 |
| 7 | user tables | `COHXUSER`, `COHYUSER`, `COHZUSER` on the grids `FREQCOH` × `DISTCOH`: bilinear in (f, D), clamped at the table ends, γ(f, 0) = 1 (D-INC-11) | available |

The Abrahamson models share the plane-wave form (f in Hz, separation ξ in metres; feet are converted
when the gravity is > 20):

```
γ_pw(f, ξ) = [1 + (f tanh(a3 ξ) / (a1 fc(ξ)))^n1]^(-1/2) · [1 + (f tanh(a3 ξ) / (a2 [fc(ξ)]))^n2]^(-1/2)   (14.2)
```

with separate horizontal (X, Y) and vertical (Z) coefficient sets; a1, a2, a3, n1, n2 and fc are
constants or functions of ξ given in the files (for example `fc = 27.9 - 4.82 ln(ξ+1) + 1.24 [ln(ξ+1) -
3.6]²` for the horizontal hard-rock model), and the factor `[fc]` is present in the 2005 model only.
Beyond the end of the formula range (150 m for the 2007 models, whose coefficient expressions become
non-physical further out) the coherency is held at its end value, the more coherent and conservative
choice, and HOUSE warns when separations leave the validated range 0-150 m.

Models 2-7 use the **directional distance** of the manual, with ΔX', ΔY' the separations along and
across Line D (the horizontal line at `WPASS <ang>` from the X axis) and α = INCOH `<alpha>` the
directionality factor:

```
D_ij = √(2 (α ΔX'² + (1 - α) ΔY'²))                                                      (14.3)
```

α = 0.5 gives the Euclidean distance; α = 0.1 weights the Y' separations three times more than the X'
separations (VP-27 checks √(2α), √(2(1 - α)) and the ratio 1/3). These models are applied only with
the wave-passage option on (HOUSE `<wpass>` = 1, EDU-26 otherwise); `V_app` = 10⁹ suppresses the delays.

### 14.3 Incoherent spatial modes

For every SSI frequency and component HOUSE forms the coherency matrix of the N interaction nodes,
`Σ_ij = γ(f, D_ij)` (horizontal projections, unit diagonal, real symmetric), and decomposes it
spectrally (LAPACK syevd):

```
Σ = Φ Λ Φᵀ,   λ1 ≥ λ2 ≥ ... ≥ 0,   Σ_k λ_k = trace Σ = N                               (14.4)
```

Column φ_k is the k-th **incoherent spatial mode** and λ_k its variance; negative round-off eigenvalues
are clipped to 0. The listing reports the trace check `|Σλ - N|/N < 1e-8` (manual Eq. 6.1; a warning
when it fails), the contribution `100 λ_1 / N` of the first mode and the number of modes that carry 90 %
of the variance (Eq. 6.4). When the
X and Y matrices are identical one decomposition serves both.

Eigenvectors are defined up to their sign. D-INC-03 orients every mode so that `Σ_i φ_ik ≥ 0`
(`EDUOPT,INCOHSIGN,ADJUST`, the default; `RAW` keeps the eigensolver's signs); a mode whose sum
vanishes (an antisymmetric mode of a symmetric layout) is oriented by its first significant component.
**Repeated eigenvalues** (mirror-image modes of a symmetric layout) leave the eigenvectors undetermined
within their eigenspace, and the deterministic sums below depend on that choice. With ADJUST the
eigenspace gets a canonical basis that depends on the subspace only (D-W3-06): first the projection of
the uniform field 1 (the coherent motion), then pivoted orthogonalisation in node order, using only the
projector `P = QQᵀ`. The factors are then reproducible and independent of round-off and of the
eigensolver. Coincident plan positions of nodes on different levels give equal rows of Σ; EDU-17 warns
when projections from different levels are closer than 0.1 × the mesh size, and `EDUOPT,INCOHMERGE,1`
merges coincident positions (D-INC-09).

### 14.4 Synthesis of the factors

INCOH `<nmodes>` selects the modes (0 = all, k > 0 = the k largest, -k = mode k only) and HOUSEX
`<supmode>` the superposition (0 Linear, 1 Quadratic):

```
stochastic simulation (SS):   s⁽ʳ⁾ = Σ_k √λ_k φ_k exp(iθ_k⁽ʳ⁾),   θ_k⁽ʳ⁾ uniform in [-RandPhz, +RandPhz]
algebraic sum (AS):           s = Σ_k √λ_k φ_k                    (zero modal phases)          (14.5)
single mode k (SRSS):         s⁽ᵏ⁾ = √λ_k φ_k                      (INCOH <nmodes> = -k)
```

The simulation is stochastic when a seed is non-zero and RandPhz > 0 (D-INC-02). With RandPhz = 180°
the phases cover the circle, `E[exp(i(θ_k - θ_l))] = δ_kl`, and every sample has the target coherency
on average: `E[s sᴴ] = Φ Λ Φᵀ = Σ`. VP-I2 runs 50 samples (FILE77001 ... FILE77050) and checks that the
mean of `s sᴴ` over Ns samples converges to Σ with the exact variance of the estimator, i.e. an RMS
error proportional to 1/√Ns (observed log-log slope -0.499). The random streams are reproducible
(D-INC-07): X and Y use the two children of `SeedSequence(HSeed)`, Z `SeedSequence(VSeed)`, one PCG64
generator per sample drawing one phase per mode at every frequency in ascending order, so sample r does
not depend on how many samples are computed together.

As f → 0, Σ → 11ᵀ, λ1 = N, φ1 = 1/√N and the AS factor is s = 1: the coherent motion and a
zero-frequency ATF of 1.00 are recovered (the manual's check, VP-27). The single modes k > 1 vanish as
f → 0, so MOTION interpolates their transfer functions to H(0) = 0.

| Approach | HOUSE input | Runs | Results |
|---|---|---|---|
| SS, "Simulation Mean" (the reference approach) | HSeed, VSeed ≠ 0, RandPhz > 0; HOUSEX `<nsim>` = Ns ≤ 50; Linear | one HOUSE run; ANALYS `<simul>` = Ns: FILE1X/Y/Z × FILE77sss → FILE8{3(s-1)+d} | MOTION per sample; mean of the ISRS (`AVERAGE`) |
| AS (deterministic) | zero seeds or RandPhz = 0; Linear | one run | FILE8 used directly; EDU-13: valid for rigid foundations only |
| SRSS TF (Quadratic) | `<nmodes>` = -k, `<supmode>` = 1 | one HOUSE + ANALYS run per mode | MOTIONX `<srss>` with `SRSSTF.txt`: `\|H\| = √(Σ_k \|H_k\|²)`, zero phase or the coherent phase (D-MOT-06) |
| SRSS FRS (Linear, single modes) | `<nmodes>` = -k, `<supmode>` = 0 | one run per mode | SRSS of the modal ISRS (line mathematics `SRSS`) |

For deeply embedded foundations the per-level approach decomposes the coherency level by level and
`BUILDFILE77` assembles the per-level FILE77 files into one; its accuracy rests on the consistent mode
signs of D-INC-03.

### 14.5 Wave passage and multiple excitation

Wave passage (HOUSE `<wpass>` = 1) delays the motion along Line D, measured from the ANALYS control
point (D-INC-06):

```
τ_i = ((x_i - x_c) cos a + (y_i - y_c) sin a) / V_app,      s_i ← s_i exp(-iωτ_i)              (14.6)
```

with `a` = WPASS `<ang>` and `V_app` = WPASS `<appv>` (> 0, Error 113). With exp(+iωt) a later arrival
is a phase lag; for coherent motion the factors are pure delays, |s| = 1 (VP-27 checks the phase
difference -ωL/V_app of two nodes L apart). **Multiple excitation** (HOUSE `<me>` = 1, which needs
`<wpass>` = 1) multiplies, after incoherency and wave passage, the factors of the interaction nodes of
zone k (ME: nodes `nfirst` ... `nlast`) by its spectral amplification ratio `SAR_k(ω)`, one value per
SSI frequency (AMP; complex with HOUSE `<cmplxspec>` = 1; modulus in [0, 10]; SAR = 1 outside the zones;
Errors 115-119; D-INC-10).

### 14.6 Application in ANALYS: free-field load (FFL) or free-field motion (FFM)

With S = diag(s) (the factor of the control-motion direction multiplies the three translations of a
node):

```
FFL (ANALYSX <ffm> = 0, default):   b_f = A_fᵀ (s ⊙ X_ff U'_f)  =  A_fᵀ S X_ff U'_f
FFM (ANALYSX <ffm> = 1):            b_f = A_fᵀ X_ff (s ⊙ U'_f)  =  A_fᵀ X_ff S U'_f                   (14.7)
```

FFL multiplies the coherent free-field load of every interaction node by its factor (the form
validated by EPRI); FFM multiplies the motion before the impedance acts on it and is meant for surface
foundations (a warning with embedded interaction nodes, D-INC-12). With s = 1 both reproduce the
coherent load exactly (VP-26: difference 0).

**FFL versus FFM (lead decision D-W3-05).** The two loads differ by the commutator,
`b_FFL - b_FFM = A_fᵀ [S, X_ff] U'_f`. For a rigid massless surface foundation under vertically
propagating waves the free field is a rigid-body translation, `U'_f = T t_d` (T of Eq. 8.3, t_d the unit
vector of the input direction), and the foundation response is `q = K_G⁻¹ Tᵀ b` with `K_G = Tᵀ X_ff T`,
so

```
K_G (q_FFL - q_FFM) = Tᵀ [S, X_ff] T t_d,        (Tᵀ [S, X_ff] T)ᵀ = -Tᵀ [S, X_ff] T             (14.8)
```

because S and X_ff are both (complex) symmetric. An antisymmetric matrix has a zero diagonal, hence:
the generalised force in the input direction, `t_dᵀ K_G q`, is the same for FFL and FFM; so is the
response of a foundation DOF that K_G does not couple to the others (the vertical translation of a
doubly symmetric mat); and coherent factors (S = I) make the commutator vanish. The horizontal
translations, coupled to rocking through K_G, and the rotations do differ: the manual's statement that
FFL and FFM give identical results for a surface rigid foundation on rock holds exactly only in these
cases. VP-26 (20 m mat on rock, hard-rock coherency, AS) checks the identities (≤ 3e-10) and Eq. (14.8)
evaluated independently from FILE11 and FILE77 (residual 4e-10); the foundation ATF differences, about
1e-3 (X) and 3e-2 (Y), are reported as informative. VP-I1 compares the ANALYS ATF of a rigid mat under
Luco-Wong coherency with the rigid-foundation average `K_G⁻¹ Tᵀ b` computed independently from FILE3,
FILE1 and FILE77, for FFL and FFM (3e-10, tolerance 1e-8): the ATF is 1 at f → 0 and decreases with
frequency (0.86 at 7.3 Hz).

**Validity limits** (manual §2.7 and §6.5.4): 3D models only, no SYMM planes (EDU-26), interaction nodes
numbered bottom-up (EDU-21 is an error for incoherent analyses), at least 200 SSI frequencies from the
start, interpolation option 6 recommended; the deterministic approaches (AS, SRSS) are validated only for
stick models on rigid mats (EDU-13).

---

## 15. Two-dimensional SSI and symmetry planes

**Code:** `sassi/core/ssi2d.py`, `sassi/core/strip2d.py` (POINT2), `sassi/core/symmetry.py`,
`sassi/modules/analys.py`, `sassi/modules/house.py`. **Verified by:** VP-42, VP-T1, VP-T2, VP-T3
(decisions D-PNT-01, D-W2-07, D-ANL-12 as corrected by D-W3-13, D-W3-12).

### 15.1 Plane-strain SSI (HOUSE `<dim>` = 1)

A 2D model lies in the global X-Z plane and is a slice of unit thickness of an infinitely long structure
and site: PLANE elements are plane-strain elements of unit thickness, and masses, stiffnesses and loads
are per unit length. The flexible-volume equation (3.1), the Schur solver, the restarts and the FILE8
layout are those of 3D; what changes is the soil:

* SITE is unchanged (one FILE1/FILE2 serves both dimensions); the in-plane input is SV (x'), P (z') or
  Rayleigh waves. SH and Love waves are anti-plane and cannot excite the in-plane model (PLANE elements
  have no UY DOF): anti-plane SH analysis is not implemented.
* POINT2 (`<dim>` = 1 is passed from HOUSE to the POINT deck) computes unit **line** loads on a
  plane-strain strip with Waas-Lysmer transmitting boundaries (section 7.6). FILE3 must come from POINT2
  for a 2D model and from POINT3 for a 3D one (D-PNT-01).
* HOUSE: the interaction DOFs are UX and UZ (FILE4 `int_dofs` = [1, 3], D-W2-07); BEAMS, SPRING and
  GENERAL nodes keep their six DOFs, whose out-of-plane components must be fixed.
* ANALYS: the flexibility is the 2 × 2 P-SV block per node pair, with the odd coupling terms multiplied
  by `sgn(x_i - x_j)` (section 8.1); the impedance is the inverse of the in-plane block, stored in the 3D
  layout with zero UY rows and columns. The free field is evaluated at `(x, y_c, z)` (the y coordinate is
  ignored) and the coordinate transformation angle must be 0 or 180°. ANALYS `<simul>` = 1 analyses the
  X and Z cases (FILE8X, FILE8Z).
* MOTION, RELDISP and STRESS (PLANE: SXX SZZ TXZ) work on FILE8 and FILE4 unchanged.

Not available in 2D: the global impedance (3D only) and incoherency, wave passage and multiple
excitation (manual §6.5.4). The manual does not recommend 2D SSI for design: a plane-strain model
radiates energy along the whole length of the structure and overestimates radiation damping, which is
unconservative.

**Verification.** VP-T1: the POINT2 far field against Kausel's exact line-load Green functions of the
same layered column (R1 V7: 0.33 % at |x| = 2R0, 0.3 % beyond; lead decision D-W3-12). VP-T3: the
2D zero-SSI identity, an embedded PLANE excavation whose structure equals the soil, gives U = U'_f at
every node for vertical SV and P and an inclined SV wave (5e-15, tolerance 1e-8). VP-42: a rigid strip
of half-width B on a layer of thickness H over a rigid base, loaded through FORCE (three load cases),
against Jakub and Roesset (1977), ν = 0.30:

```
K_x / G = 1.175 (1 + 2.15 B/H),        K_φ / (G B²) = 2.394 (1 + 0.17 B/H)        (1/8 ≤ B/H ≤ 1/2)   (15.1)
```

The FE stiffness converges as O(h) (singular edge stresses of a rigid punch), so the criterion uses the
Richardson value `2K(32) - K(16)` (observed mesh-difference ratio about 2): within 5 % (observed 0.3 % to
4.5 %); the chain FORCE + ANALYS equals the direct rigid-strip impedance `Tᵀ X_ff T` of FILE3, and K_x
decreases as H grows (it vanishes for a half-space in plane strain).

### 15.2 Symmetry and antisymmetry planes (SYMM)

A structure symmetric about a vertical plane, on a horizontally layered site, under a loading that is
symmetric or antisymmetric about that plane, has a response with the same property, so a half model
(one plane) or a quarter model (two orthogonal planes, parallel to XZ and YZ; in 2D a line parallel to Z)
can be analysed. With n the plane normal, the reflection `P = I - 2nnᵀ` and the sign `s = +1`
(symmetric, SYMM type 0) or `-1` (antisymmetric, type 1), the field obeys `u(x') = s P u(x)` for the
mirror image x' of x. HOUSE fixes on the plane nodes the DOFs that this makes zero: for symmetric
loading the normal translation and the two in-plane rotations, for antisymmetric loading the two
in-plane translations and the normal rotation (rotations are pseudo-vectors, D-ANL-12). A horizontal X
input is antisymmetric about a plane normal to X and symmetric about a plane normal to Y; a vertical
input is symmetric about both.

**The soil: image superposition (D-W3-13).** ANALYS needs the flexibility of the *reduced* interaction
set. The interaction forces of the full model have the symmetry of the loading, `f(j_S) = s_S P_S f(j)`,
so the displacement at a reduced-model node i is the sum over the 2^k combinations S of the k planes
(the empty combination is the node itself):

```
F_red(i, j) = Σ_S s_S F(i, j_S) P_S                                                          (15.2)
            = F(i, j) + s F(i, j') P                    (one plane)
```

with j_S the image of j in the planes of S and s_S, P_S the products of their signs and reflections. (The
earlier text `F(i,j) ± P F(i,j') P` of D-ANL-12 was wrong; for a layered site `F(i, j') P = P F(i', j)`,
so the reflection belongs on the load side.) A node on a plane carries half (on two planes a quarter) of
the full-model force, exactly as the structure, the excavated soil, the masses and the loads of a half
model carry half of that node's full-model values. For the DOFs the symmetry constrains, `(I + sP)` is
zero: their rows and columns of F_red vanish and they are removed **before** the inversion
`X_red = F_red⁻¹`; F_red is symmetrised like F_ff. ANALYS checks at every frequency that the free field
has the declared symmetry, `U'(image) = s P U'`, and stops otherwise.

**Modelling rules.** The model lies on one side of every plane; elements, masses and loads lying *in* a
plane carry 1/2 (1/4 on two planes) of their full-model values. Incoherency, wave passage, multiple
excitation and the global impedance are not allowed with SYMM (manual §2.7, G-17). **VP-T2** analyses a
flexible surface mat with a stick and an embedded FV basement with a stick as full, half and quarter
models under X and Z input: the reduced models reproduce the full model at every DOF (≤ 4e-11, tolerance
1e-8), and the DOFs the symmetry fixes vanish in the full model.

---

## 16. Nonlinear soil SSI iterations

**Code:** `sassi/core/shake.py` (free field), `sassi/core/nlsoil.py`,
`sassi/prep/commands/nlsoil_cmds.py` (near field), `sassi/modules/house.py`, `sassi/modules/stress.py`.
**Verified by:** VP-04 (SOIL), VP-N1 (the near-field iterations reproduce SOIL), VP-N2 (restart
consistency), example 6.

The frequency-domain solution is linear. Soil nonlinearity is approximated by the **equivalent-linear
method**: each soil element gets a shear modulus and a damping ratio compatible with its own effective
strain, found by iteration (manual §4.1.2 item 18).

* **Primary (free-field) nonlinearity**: SOIL iterates the layer properties of the site (section 12);
  SITE and the excavated soil use them.
* **Secondary (near-field) nonlinearity**: the extra straining caused by the structure (soft backfill
  next to an embedded wall, soil under heavy neighbouring buildings). The near-field soil is modelled
  with SOLID (3D) or PLANE (2D) elements of the **structure** (ETYPE 1), next to the excavated free-field
  soil of the flexible-volume method, and iterated element by element.

One near-field iteration (requirements §2.4, D-NLS-01 ... D-NLS-04):

```
HOUSE    iteration 0 (.liq ≠ 1):  G = GFAC · G_layer,   β_s = DFAC · β_s,layer,   β_p = DFAC · β_p,layer
         later (.liq = 1):        G = G_max · (G/G_max)(γ_eff),   β_s = β_p = D(γ_eff)   (curve ICURVE of FILE73)
         one internal material per element (density and Poisson's ratio of the element material);
         writes FILE78 (the properties used)
ANALYS   New Structure restart (<mode> 1): the impedance X_ff of the initiation run (COOXqqq) is reused
STRESS   (<iter> = 1, once per input direction) strain histories of the nonlinear elements,
         γ_eff = ESF · max_t |γ(t)|  → FILE74 (one per direction: FILE74X, FILE74Y, FILE74Z)
COMBXYZSTRAIN   γ_eff = √(γ_X² + γ_Y² + γ_Z²) per element  → FILE74
converged when  max |ΔG/G| < 2 %  and  max |Δβ| < 0.5 % (absolute),  at most 8 iterations            (16.1)
```

**What the material and the .pin factors mean** (requirements §4.4 item 7). The M-table material of a
nonlinear element is the **low-strain** near-field soil: density, Poisson's ratio and `G_max = ρVs²`,
the reference of the G/G_max curve (never strain-compatible values, which the curve would soften a
second time). `G_layer`, `β_layer` are the free field at the element: the TOPL layer that contains the
element's mid-depth, with its L properties. GFAC = DFAC = 1 ("same shear modulus as in free-field")
starts the near field from the free-field state; a looser backfill starts lower (GFAC < 1). GFAC and
DFAC choose the starting point only: as in SHAKE, the converged strain-compatible state is set by G_max,
the curves and the motion. When G changes, Poisson's ratio is kept (the constrained modulus follows G)
and the P-wave damping equals the shear damping, as in SOIL (D-SOL-06). Curve interpolation is that of
SOIL: linear in log10(γ), constant outside the tabulated range.

**Strain measures** (ISTR; D-NLS-03): 0 = the largest engineering shear-strain component
`max(|γ_xy|, |γ_xz|, |γ_yz|)` (|γ_xz| in 2D), the SHAKE measure; 1 = the octahedral shear strain
`(2/3)√[(ε_x-ε_y)² + (ε_y-ε_z)² + (ε_z-ε_x)² + 1.5(γ_xy² + γ_yz² + γ_xz²)]` (the maximum shear strain
`√[(ε_x-ε_z)² + γ_xz²]` in 2D). For simple shear γ, ISTR 0 gives γ and ISTR 1 gives `√(2/3)γ`.

**Convergence measure.** Row k of `NLSOIL_CONVERGENCE.TXT` compares the properties *used* in iteration
k (FILE78) with the strain-compatible properties of its response (FILE74): `ΔG/G = (G_new - G_used)/G_new`
(relative to the new value, as in SHAKE) and `Δβ = β_new - β_used` in percentage points.
`NLSSIITER,<variables>,<maxit>` runs the command lists until (16.1) holds; `NLSSIRESET` starts a new
analysis (`<model>.liq` = 0). The files carry the FILE4 hash of the properties they belong to, so a
stale strain file of another iteration is refused.

Each iteration is a New Structure restart (2-4 times faster than an initiation run). **VP-N1**: a
laterally uniform column of near-field soil inside an FV excavation of a uniform deposit, with the free
field and the excavated soil at the strain-compatible properties of a SOIL run of the same column,
motion and curves. Without a structure there is no secondary nonlinearity, so (a) the iterations started
from the low-strain properties converge monotonically (max |ΔG/G| = 318, 42, 10, 2.6, 0.65 %) to the
SHAKE profile, within 5 % on G/G_max and 10 % on the damping of every layer (observed 1.3 % and 0.5 %),
and (b) started from the free field (GFAC = DFAC = 1) they are converged at once. **VP-N2**: an iteration
whose properties do not change reproduces the initiation FILE8 through the New Structure restart (0,
tolerance 1e-10), and repeating HOUSE with the same strains reproduces FILE78, FILE8 and FILE74.

---

## 17. Option NON: nonlinear structures by equivalent linearisation

**Code:** `sassi/core/hysteresis.py` (hysteresis models, equivalent-linear properties),
`sassi/core/panels.py` (panel kinematics, SHEAR, BBCGEN), `sassi/modules/nonlinear.py` (NONLINEAR
module, COMB_XYZ_THD), `sassi/prep/commands/nonlinear_cmds.py`. **User guide:**
[OPTION_NON.md](../user/OPTION_NON.md). **Verified by:** VP-45, VP-46, VP-NON1, example 7
(requirements §4.15; decisions D-NON-01 ... D-NON-12, D-W3-08, D-W3-09).

The SHAKE idea applied to the structure: a nonlinear wall panel or spring is replaced by an
equivalent-linear element whose stiffness and damping depend on its deformation amplitude, and the SSI
analysis is repeated until they no longer change (manual §1.5.4, Fig. 1.2).

**Deformations.** A panel is a group of coplanar SHELLs in a vertical plane; only its four corner nodes
(BL, BR, TR, TL in the panel axes e_h horizontal and e_v = Z) are used, with `u = d·e_h`, `w = d·e_v`,
the width L and the height H:

```
γ   = ½[(u_TL - u_BL) + (u_TR - u_BR)]/H + ½[(w_BR - w_BL) + (w_TR - w_TL)]/L        shear strain (Disp. Opt 1)
κ   = [(w_TR - w_TL) - (w_BR - w_BL)] / (L H)                                       curvature (Disp. Opt 2)
ε_v = ½[(w_TL - w_BL) + (w_TR - w_BR)]/H                                            axial strain (informative)   (17.1)
```

A rigid translation or in-plane rotation gives zero for all three. A nonlinear spring uses its
elongation `u_J - u_I` along one DOF. The histories come from RELDISP (`.THD`, free-field reference,
which cancels in every difference) for the X, Y and Z inputs, added time step by time step by
COMB_XYZ_THD (D-NON-08), because the nonlinear behaviour must be driven by the simultaneous
three-component input (manual §1.5.4).

**Equivalent-linear properties.** For every element (requirements §4.15):

```
x_eq   = EDF · max_t |x(t)|                                   EDF ≈ 0.8 (0.7-0.9)
E_new  = E_el · K_sec(x_eq) / K_el,   K_sec = F_bb(x_eq)/x_eq,   K_el = Y1/X1  (first BBC slope, D-NON-11)
ξ_h    = E_D / (4π E_S),   E_S = x_eq F(x_eq)/2                (stabilised loop at x_eq)
ξ      = min(cutoff, scale · ξ_h + [ξ_el])                     (scale 0 → 1, cutoff 0 → none; D-NON-04)    (17.2)
```

F_bb is the backbone curve (BBC): odd, piecewise linear through the origin and the user points (point 1
= the cracking point, the end of the elastic range), held constant beyond the last point. For a panel the
modulus of its own material is scaled with Poisson's ratio kept, so shear, axial and bending stiffness
degrade together (each panel needs a material of its own, D-NON-12); for a spring the spring constant of
its DOF is scaled. E_D is the area of the third of three symmetric cycles at ±x_eq through the hysteresis
model. The ductility `μ = max|x| / x_cr` and the force-reduction factor `F_μ = K_el max|x| / |F(t*)|` at
the time t* of max|x| (the inelastic absorption factor of ASCE 43-05) are reported.

**Hysteresis models** (BBC type / Force Opt, D-NON-03):

* **GMR, General Masing Rule** (springs): virgin loading on the BBC; a branch from the reversal point
  (x_r, F_r) is `F = F_r + 2 F_bb((x - x_r)/2)` (Masing 1926), with the extended memory rules (Pyke 1979;
  Kramer 1996): a branch that reaches the previous reversal continues on the branch it left, one that
  reaches the largest past excursion joins the backbone. Symmetric loops are closed with
  `E_D = 8 ∫₀ˣ F dx - 4 x F(x)`; for an elastic-perfectly-plastic BBC `ξ_h = 2(x - x_y)/(πx)` and
  `K_sec = F_y/x` (VP-45, to round-off).
* **CMS, Cheng-Mertz Shear** (panels; low-rise RC walls): the S1 model transcribed rule by rule from
  subroutine HYST04 of INRESB-3D-SUP (Cheng and Mertz 1989): elastic until cracking, loading on the
  multi-linear backbone; unloading in three force bands with the degrading stiffnesses
  `S1 = SI min(1.4675 (DC/DMAX)^0.345, 1)`, `S2 = SI min(0.7761 (DC/DMAX)^0.5195, 1)`,
  `S3 = SI min(0.0707 (DC/DMAX)^1.369, 1)` (SI = PC/DC the initial stiffness); pinched reloading below
  0.75 PC with the slip stiffness `SR = SI min((DC/DMAX)^1.02, 1)`; reloading toward 0.95 PMAX on the
  first unloading branch and then toward (1.04 DMAX, PMAX); up to ten stored small loops. The original is
  load-incremental; here every rule is a straight branch followed exactly in displacement control.
* **TAK, Takeda** (experimental, `EDUOPT,NONEXT,1`): trilinear backbone through the cracking and yield
  points, flat beyond; unloading after yield with `K_u = K_cy (D_y/D_max)^0.4`,
  `K_cy = (P_y + P_c)/(D_y + D_c)` (Takeda, Sozen and Nielsen 1970; Otani 1974). The manual allows TAK
  only with the elastic damping, so `ξ = ξ_el`.
* CMB (Cheng-Mertz Bending) is not included in the manual's version and is refused.

**Iteration and convergence.** NONLINEAR writes the new properties into `<model>_new.hou`; the next SSI
analysis is a New Structure restart. The state machine of the manual: without `PANEL.NON` /
`SPRING.NON` the run is the *elastic run* and creates the file; with `.NON` = 1 it is an *iteration*
that compares the properties of the analysis just made (`*_EQL_Matl_Prop.txt`) with the new ones.
Convergence (D-NON-06): `max |E_new - E_old| / E_old < 2 %` (relative to the previous value) and
`max |ξ_new - ξ_old| < 0.5 %` (absolute), at most 10 iterations. When the wall shear is set by the
inertia forces (force control) the fixed point is where the backbone force equals the demand F_d; one
iteration `x ← F_d x / F_bb(x)` reduces the distance to it only by the factor `1 - x F_bb'/F_bb`, close to
1 just after cracking, where the backbone is nearly flat: hence the manual's advice to use smooth BBCs.

**Shear capacities and backbones** (SHEAR, BBCGEN; spec 10 §3.11). ACI 318-08, Wood (1990), Barda et
al. (1977) and Gulec and Whittaker (2009) capacities with √f'c in psi (the formulas are listed in
[OPTION_NON.md §5](../user/OPTION_NON.md#5-shear-capacities-and-backbone-generation-shear-bbcgen)). BBCGEN
(D-NON-09) writes a 22-point CMS curve: point 1 = cracking `(V_cr/(G A_W), V_cr)` with
`V_cr = 3√f'c A_W` (ASCE 4-17) or `CFL·V_u`; points 2-21 equally spaced in strain up to the yield point
`(0.004, V_u)` on `V = V_cr + (V_u - V_cr)[1 - (1 - ξ)²]`, `ξ = (γ - γ_cr)/(γ_y - γ_cr)` (zero slope at
yield); point 22 = failure `(0.02, 1.02 V_u)`. Its first slope is `G A_W` by construction.

**Verification.** VP-45: GMR loops of an elastic-perfectly-plastic BBC at 1.25 ... 50 x_y (closed, ξ_h and
K_sec to 1e-6, observed 2e-16) and the NONLINEAR state machine. VP-46: the SHEAR numbers of the spec 10
worked example in British units (1e-6) and SI (within half the printed last digit, 0.05 kN; lead
decision D-W3-08) and the BBCGEN curve. VP-NON1: an SDOF with a GMR spring on a rigid base under harmonic
input, iterated through HOUSE → ANALYS restart → MOTION → RELDISP → NONLINEAR: every response equals the
closed-form steady state and the iterations follow an independent closed-form fixed-point iteration
(1e-8, tolerance 1e-6), converging in 5 iterations to the secant stiffness and Masing damping at the
converged amplitude.

---

## 18. SOIL-NON: nonlinear time-domain site response

**Code:** `sassi/core/soilnon.py` (decisions SN-1 ... SN-12 in its docstring), `sassi/modules/soil.py`,
`sassi/prep/commands/soilnon_cmds.py` (NLSOIL, NLSLAYER, DELNLS; side file `<model>.nls` written by
AFWRITE). **Verified by:** VP-SN1, VP-SN2, VP-SN3 (informative) (lead decision D-W3-10).

`NLSOIL,1` makes SOIL integrate the soil column in the time domain instead of iterating
equivalent-linear properties (SOIL-EQL, section 12). The manual cites DEEPSOIL (Hashash and Park 2001)
as the theory; the details it leaves open are taken from that literature.

**Model.** The column above the half-space is a chain of shear elements with lumped masses. Each element
follows the modified hyperbolic (MKZ) backbone of the manual (Matasovic and Vucetic 1993), with the
extended Masing rules for unloading and reloading (after a reversal at (γ_R, τ_R):
`τ = τ_R + 2F((γ - γ_R)/2)`; a branch that reaches the previous reversal continues on the curve of the
previous cycle; beyond the largest past strain the backbone is followed):

```
τ = G0 γ / (1 + β (|γ|/γ_r)^s)   [+ η dγ/dt]                                                (18.1)
```

NLSLAYER gives β, s, γ_r (percent) and the viscosity η per SOIL sublayer, or `curvefit` = 1 fits the
backbone to the sublayer's G/G_max curve. β and γ_r enter only through β/γ_r^s, so the fit sets β = 1 and
returns γ_r and s, with s restricted to [0.2, 1] so that the stress never decreases (SN-3). The
hysteretic damping of a Masing loop of amplitude γ_a is

```
D = (2/π) [2 ∫₀^γa F dγ / (γ_a F(γ_a)) - 1],    ∫₀^γa F dγ = G0 γ_a²/2 · ₂F₁(1, 2/s; 1 + 2/s; -z),   z = β(γ_a/γ_r)^s   (18.2)
```

(for s = 1 the closed form of the hyperbolic model, `D = (2/π)[2(1+z)(z - ln(1+z))/z² - 1]`). It is zero
at small strains, so a small-strain viscous damping ξ_min (the DYNP damping at the smallest strain of the
curve) is added (SN-4) as NLDampType 1, a classical modal damping matrix
`C = M Φ diag(2ξ_n ω_n) Φᵀ M` of the small-strain column on a fixed base with strain-energy weighted
modal ratios (frequency independent at the modal frequencies, SN-5); 2, element dashpots η/h; or 3,
Rayleigh `C = α M + β K0` (matched at f1 and 5 f1 when the multipliers are 0).

**Base and equation of motion.** BedInt 0: rigid base moving with the input (a *within* motion). BedInt 1:
elastic half-space by the Joyner and Chen (1975) dashpot `c = ρ_r V_r` per unit area with the *outcrop*
motion as input ("currently disabled" in ACS SASSI V3, implemented here, SN-7). Relative to the input
frame:

```
M ü + C u̇ + R(u) + c_b u̇_base = -M 1 a_in(t)                                                (18.3)
```

integrated by the implicit Newmark average-acceleration method (γ = 1/2, β = 1/4) with Newton-Raphson
iterations and the consistent tangent of the Masing branches (force and displacement norms ≤ 1e-6, at
most 25 iterations, sub-steps bisected up to 2¹² times; SN-10). Sub-increments: NSTimeSunInc fixed
steps, or (0) flexible steps of at most 1/400 s, repeated when a strain increment exceeds 0.005 % (SN-9).
Each sublayer is divided into an odd number of elements of thickness `h ≤ Vs/(10 f_max)`,
f_max = min(25 Hz, Nyquist) (SN-11). Outputs: `ACCxxx.TH`, `SNxxx.TH`, `SSxxx.TH` and FILE88 with the
secant G and the damping ξ_min + Masing damping at the effective strain, so SITE can use the properties
of the model that was integrated (SN-12). Every soil sublayer needs an NLSLAYER set (Error 125), the
half-space does not.

**Verification.** VP-SN1: in the linear limit (γ_r = 10⁶ %) the time-domain column reproduces the
frequency-domain linear SOIL solution of the same profile and motion (elastic base, outcrop input:
surface PGA and 5 % spectra within 2 %, observed ≤ 0.4 %; the discrete column's amplification peak within
1 % of the continuum). VP-SN2: the Masing loop of one element has the backbone secant modulus and the
loop-area damping of Eq. (18.2) (1e-3; the FILE88 damping 1e-6). VP-SN3 (informative): the SHAKE91
sample problem with SOIL-NON and SOIL-EQL; at this moderate shaking (shear-strain index 0.04 %) the
surface PGA differs by about 12 %, the order of difference that Kim et al. (2016) report.

---

## 19. Water: hydrodynamic mass by soft solids

**Code:** `sassi/prep/water_lib.py`, `sassi/prep/commands/water.py` (FILLPOOL, REFINEMODEL,
LISTPOOLINTER, MERGEPOOL, POOLDATA). **User guide:** [WATER.md](../user/WATER.md). **Verified by:**
VP-54W, VP-W1 (decisions D-WAT-01 as amended by D-W3-07).

ACS SASSI models the water of pools and tanks with ordinary SOLID elements that have the properties of
water (FILLPOOL fills a pool sub-model with them and connects them to the walls with springs). An
isotropic elastic solid obeys `ρü = (K + 4G/3) ∇(∇·u) - G ∇×(∇×u)`; for G → 0 the only stress is the
pressure `p = -K ∇·u`, and

```
ρ ü = -∇p,     p̈ = c² ∇²p,     c = √(K/ρ) ≈ 1480 m/s                                          (19.1)
```

the linear acoustic equation of an inviscid compressible fluid, written in displacements. Far below the
first compression frequency of the water column, `c/(4H)` (74 Hz for H = 5 m), the water is
incompressible, `∇²p = 0`, with `∂p/∂n = -ρ a_n` on the walls and the floor (imposed by the interface
springs along the wall normals, frictionless along the walls) and `p = 0` at the free surface: the
potential-flow problem of Westergaard (1933) and Housner (1963). The water pushed by the walls adds the
**impulsive mass**; for a rigid rectangular tank of length 2L in the direction of motion and depth H the
exact series is

```
m_i / m = 2/(L H²) Σ_n tanh(λ_n L)/λ_n³,     λ_n = (2n + 1)π/(2H)                               (19.2)
```

(exactly ½ for L = H; Housner's closed form `tanh(√3 L/H)/(√3 L/H)` is 6-11 % conservative).

**Why G = 10⁻⁸ K (D-W3-07).** D-WAT-01 proposed ν = 0.49, i.e. G = 0.0201 K. Such a "water" is an elastic
solid whose shear modes lie in the seismic band (about 10 Hz for 5 m of water); below them it moves
rigidly with the tank and gives the **total** mass instead of the impulsive mass (F/(m a) = 1.00 at
0.2-2 Hz for L = H, where the exact value is 0.50). The smaller G, the lower the spurious shear modes:
with G = 10⁻⁸ K (Vs = 10⁻⁴ Vp) they lie below about 0.05 Hz for metre-sized pools and the water behaves
as an inviscid fluid from 0.2 Hz up; K/G = 10⁸ is far from the double-precision limit, and at every
non-zero frequency the inertia controls the near-zero-stiffness shear deformations. The water SOLIDs
must use **incompatible modes** (FILLPOOL and MERGEPOOL set `MOPT,0`): a trilinear hexahedron with 2×2×2
Gauss points locks volumetrically for a nearly incompressible material (with `MOPT,1` F/(m a) is again
about 1.0 below 1 Hz and wrong at higher frequencies). Oblique walls use the diagonal of the coupled interface spring matrix
`Stiff P_N + stiff2 (I - P_N)`, because zero-length GENERAL elements are rejected (Error 8, D-CHK-07).
Sloshing (the convective mass, which needs the restoring force ρg at the free surface) is not modelled;
FILLPOOL prints Housner's estimate of it.

**Verification.** VP-54W: REFINEMODEL counts, areas and volumes; FILLPOOL water volume `a·b·(H - m·h)`,
springs, offsets, the water material (K = 2.2 GPa, G/K = 10⁻⁸), LISTPOOLINTER and MERGEPOOL. VP-W1: a rigid
tank 10 m × 2 m filled to 5 m on a practically rigid site, through SITE → POINT → HOUSE → ANALYS, with SHELL
walls and with SOLID walls on a slab: the water mass equals ρV (1e-10), the interface forces balance the
inertia of the water (1e-6), the tank follows the control motion (1e-3), and the base shear ratio
F/(m a) = 0.507-0.509 at 2-8 Hz equals the exact impulsive mass ratio 0.5 within 5 % (the +1.5 to +1.9 %
are the mesh error of 0.5 m elements and a small compressibility effect); Housner and Westergaard are
informative.

---

## 20. Option A: SSI loads for a second-step ANSYS analysis

**Code:** `sassi/modules/loadgen.py` (module LOADGEN), `sassi/prep/commands/loadgen_cmds.py` (LOADGEN,
LOADGENDYN, LGFILE, LGNODE, LGTIME, LGMAP, LGOPT, LGLIST, RUNLOADGEN). **User guide:**
[OPTION_A.md](../user/OPTION_A.md). **Verified by:** VP-LA1, VP-LA2 (lead decision D-W3-11;
interpretations D-LGN-01 ... D-LGN-08).

Option A is the two-step approach: step 1 is the SSI analysis; step 2 a refined ANSYS model of the
structure loaded by the motion that step 1 computed. Split the total displacement into the motion of a
reference frame and the motion relative to it, `u = u_r + ι u_g`, with `ι u_g` a rigid-body translation
(the free-field ground motion, or the motion of a reference node). A rigid-body translation strains
nothing (`K ι = 0`), so with the interface DOFs b prescribed the structural DOFs s obey

```
M_ss ü_r,s + C u̇_r,s + K_ss u_r,s = -M_ss ι_s a_g(t) - K_sb u_r,b(t)                             (20.1)
```

* **Dynamic second step** (LOADGENDYN, method REL): an ANSYS transient analysis with `ACEL = a_g(t)` (the
  acceleration of the reference frame; ANSYS applies -M·ACEL to every mass) and the interface
  displacements relative to the ground prescribed with D (RELDISP `.THD` with the free-field reference,
  started at rest: the zero-mean offset of the periodic FFT solution is removed). LOADGEN writes them as
  APDL TABLE arrays. Method ACC fixes the interface and drives it with the absolute acceleration of a
  reference node (valid for a rigid foundation).
* **Equivalent static second step** (LOADGEN): freeze Eq. (20.1) at a critical time t* and drop the
  damping and velocity terms, `K_ss u_s = -M_ss a_s(t*) - K_sb u_b(t*)`: the nodal inertia forces of the
  **absolute** accelerations, `F = -m a(t*)` (MOTION `.ACC`), plus the interface displacements at t*.
  Relative and total interface displacements differ by a rigid translation, which produces no stress.
  Without the displacements ("Acceleration") the interface is fixed: the classical fixed-base equivalent
  static analysis. The critical times are the largest local maxima of the base shear (or of a moment or a
  node response).

The masses are the row sums of the HOUSE structure mass matrix per direction,
`m_(i,d) = Σ_j M[(i,d), (j,d)]` (exact for lumped masses, the total mass per direction preserved for
consistent ones), or load nodes of the ANSYS model slaved to master nodes of the SSI model,
`a_k = a_M + α_M × (x_k - x_M)` (a coarse stick → a refined ANSYS floor). Rayleigh damping can be given
by a damping ratio ζ at two frequencies: `α = 2ζω1ω2/(ω1 + ω2)`, `β = 2ζ/(ω1 + ω2)`.

**What is approximate.** For a linear structure with the same mass, stiffness and damping as the SSI
model the decomposition is exact. The dynamic step replaces the frequency-independent hysteretic damping
of SASSI by Rayleigh damping, `ζ(ω) = α/(2ω) + βω/2`, exact at two frequencies only ("a significant
limitation", manual); the static step represents one instant without damping forces; and anything added
in the ANSYS model that changes the interface motion violates the assumption of Option A.

**Verification.** VP-LA1: a three-mass stick on a rigid mat; the inertia forces of the APDL sum to minus
the STRESS base shear and their moment to minus the base moment at the critical time (1 %, observed
0.04 %: MOTION interpolates nodal transfer functions, STRESS the element force transfer function), the
critical time is that of the largest STRESS base shear, the generated masses equal the MT masses (1e-12)
and the interface D equal RELDISP (1e-10). VP-LA2: the APDL tables reproduce the control motion, MOTION
`.ACC` and RELDISP `.THD` (≤ 5e-12, tolerance 1e-10), and the exported second step, solved with the HOUSE
matrices at every Fourier frequency, reproduces the SSI absolute accelerations to 5e-9 (tolerance 1e-6),
the exactness of Eq. (20.1). The APDL is checked by a parser; it has not been run in ANSYS.

---

## 21. Spectrum-compatible motions (EQUAKE)

**Code:** `sassi/core/equake_lib.py`, `sassi/modules/equake.py`. **Verified by:** VP-36 (RG 1.60,
SRP 3.7.1 criteria), VP-30 (response spectra).

1. **Target**: a design spectrum (Hz, SA in g) at damping ζ, interpolated log-log; constant spectral
   displacement below its first frequency, constant SA above its last.
2. **Initial motion**: a stationary Gaussian process with random phases (PCG64, seed `<rand>`) and
   Fourier amplitudes from the SIMQKE power-spectral-density estimate (Gasparini and Vanmarcke 1976),
   multiplied by an envelope with a stationary strong-motion part; or the Fourier phases of a seed
   record (`<accopt>` = 1).
3. **Frequency-domain matching** (Levy-Wilkinson): Fourier amplitudes multiplied by
   `SA_target/SA_computed` (Nigam-Jennings spectra at 100 points per decade), phases kept; with random
   phases each iteration also clips the stationary process at the target zero-period acceleration
   (iterative clipping and filtering), which lowers the crest factor so that the PGA can equal the
   target ZPA.
4. **Time-domain refinement** (P1): Al Atik and Abrahamson (2010) tapered-cosine wavelets added at the
   response-peak times to move the peaks to the target.
5. **Drift control and baseline correction** in every step (D-EQK-03, D-W1-09): a zero-phase high-pass
   filter on the frequency-domain updates, zero-mean acceleration, and a degree-5 displacement
   polynomial fitted by constrained least squares (zero final velocity and displacement).
6. **Checks** (D-EQK-04): the manual's criteria and SRP 3.7.1 Rev. 4 Option 1 Approach 2, reported
   separately: duration ≥ 20 s; Δt ≤ 0.005 s (manual) and Nyquist ≥ 50 Hz (SRP); spectrum at ≥ 100
   points per decade never more than 10 % below or 30 % above the target, at most 9 adjacent points
   below; strong-motion duration (Arias 5-75 %) ≥ 6 s; cross-correlation |ρ| ≤ 0.16 between
   components; PSD ≥ 80 % of the target PSD over 0.3-24 Hz when a target PSD is given.

---

## 22. Code and verification map

| Topic | Theory section | Code | Verification problems |
|---|---|---|---|
| complex moduli | 2 | `conventions.py`, `elements/base.py` | VP-01, VP-02a, VP-39 |
| flexible-volume equation, method variants | 3 | `core/ssi_solver.py`, `modules/analys.py`, `modules/house.py` | VP-15, VP-16, VP-E1, VP-S2 |
| thin-layer method | 4 | `core/tlm.py`, `modules/site.py` | VP-05, VP-06, VP-43 |
| half-space | 5 | `core/tlm.py` | VP-02b, VP-08, VP-11 |
| free field, normalisation | 6 | `core/freefield.py`, `modules/site.py` | VP-02b, VP-15, VP-16, VP-24 |
| POINT3 / POINT2 | 7 | `core/axisym.py`, `core/strip2d.py`, `core/greens_tlm.py`, `modules/point.py` | VP-07, VP-08, VP-09, VP-T1 |
| flexibility, impedance, global impedance | 8 | `core/flexibility.py`, `core/ssi_solver.py` | VP-07, VP-10, VP-11, VP-13, VP-14, VP-17 |
| ANALYS solution, restarts, cases | 9 | `core/ssi_solver.py`, `modules/analys.py`, `modules/force.py`, `modules/combin.py` | VP-22, VP-23, VP-24, VP-25, VP-40, VP-41 |
| interpolation | 10 | `core/interp.py` | VP-28, VP-32, VP-33, VP-53 |
| convolution, spectra, RELDISP, STRESS | 11 | `core/signal.py`, `core/spectra.py`, `modules/motion.py`, `modules/reldisp.py`, `core/stress_lib.py`, `modules/stress.py` | VP-29, VP-30, VP-31, VP-34, VP-35, VP-S1, VP-S2 |
| SHAKE | 12 | `core/shake.py`, `modules/soil.py` | VP-02a, VP-04 |
| elements | 13.1 | `elements/*.py`, `modules/house.py` | VP-03, VP-37, VP-38, VP-39, VP-39b, VP-H1, VP-H2, VP-A1 |
| TSHELL thick shell, THSHLSTR | 13.2 | `elements/tshell.py`, `modules/stress.py` | VP-38T, VP-54T, VP-TS1 |
| node-numbering optimizer | 13.3 | `core/renumber.py`, `modules/house.py` | VP-O1 |
| incoherency, wave passage, multiple excitation, FFL/FFM | 14 | `core/coherency.py`, `core/incoherency.py`, `modules/house.py`, `core/ssi_solver.py`, `modules/analys.py`, `data/coherency/*.json` | VP-26, VP-27, VP-I1, VP-I2 |
| 2D plane-strain SSI | 15.1 | `core/ssi2d.py`, `core/strip2d.py`, `modules/analys.py` | VP-42, VP-T1, VP-T3 |
| symmetry planes (SYMM) | 15.2 | `core/symmetry.py`, `modules/house.py`, `modules/analys.py` | VP-T2 |
| nonlinear soil SSI iterations | 16 | `core/nlsoil.py`, `prep/commands/nlsoil_cmds.py` | VP-N1, VP-N2 |
| Option NON | 17 | `core/hysteresis.py`, `core/panels.py`, `modules/nonlinear.py`, `prep/commands/nonlinear_cmds.py` | VP-45, VP-46, VP-NON1 |
| SOIL-NON | 18 | `core/soilnon.py`, `modules/soil.py`, `prep/commands/soilnon_cmds.py` | VP-SN1, VP-SN2, VP-SN3 |
| water modelling | 19 | `prep/water_lib.py`, `prep/commands/water.py` | VP-54W, VP-W1 |
| Option A (LOADGEN) | 20 | `modules/loadgen.py`, `prep/commands/loadgen_cmds.py` | VP-LA1, VP-LA2 |
| EQUAKE | 21 | `core/equake_lib.py`, `modules/equake.py` | VP-36 |

The independent references built for the impedance problems (a static boundary-element method with
Boussinesq/Cerruti kernels and a dynamic one with exact Lamb kernels) are described in
[docs/verification/impedance_study.md](../verification/impedance_study.md).

---

## 23. References

The full lists, with links where open copies exist, are in R1 §11 and R2 §L. The main sources of this
manual:

* Lysmer, J., Tabatabaie-Raissi, M., Tajirian, F., Vahdani, S., Ostadan, F. (1981). *SASSI - A System
  for Analysis of Soil-Structure Interaction*. Report UCB/GT/81-02, University of California, Berkeley.
* Lysmer, J., Ostadan, F., Chin, C.C. (1999). *SASSI2000 Theoretical and User's Manual*. UC Berkeley.
* Kausel, E. (1981). *An Explicit Solution for the Green Functions for Dynamic Loads in Layered Media*.
  MIT Research Report R81-13 (Kausel and Peek 1982, BSSA 72(5)).
* Kausel, E., Roesset, J.M. (1981). Stiffness matrices for layered soils. *BSSA* 71(6).
* Lysmer, J., Waas, G. (1972). Shear waves in plane infinite structures. *J. Eng. Mech. Div. ASCE* 98(EM1).
* Lysmer, J., Kuhlemeyer, R.L. (1969). Finite dynamic model for infinite media. *J. Eng. Mech. Div.
  ASCE* 95(EM4).
* Tajirian, F.F., Tabatabaie, M. (1985). Vibration analysis of foundations on layered media. ASCE
  Convention, Detroit.
* Kagawa, T., Mejia, L.H., Seed, H.B., Lysmer, J. (1981). *TLUSH*. Report UCB/EERC-81/14.
* Schnabel, P.B., Lysmer, J., Seed, H.B. (1972). *SHAKE*. Report EERC 72-12 (with the 1973 Udaka-Lysmer
  revision); Idriss, I.M., Sun, J.I. (1992). *User's Manual for SHAKE91*.
* Ostadan, F., Deng, N., Roesset, J.M. (2004). Estimating total system damping for soil-structure
  interaction systems. 3rd UJNR Workshop on SSI.
* Nigam, N.C., Jennings, P.C. (1969). Calculation of response spectra from strong-motion earthquake
  records. *BSSA* 59(2).
* Veletsos, A.S., Verbic, B. (1973). Vibration of viscoelastic foundations. *EESD* 2; Pais, A.,
  Kausel, E. (1988). Approximate formulas for dynamic stiffnesses of rigid foundations. *Soil Dyn.
  Earthq. Eng.* 7(4).
* Ghiocel, D.M. (2019). *ACS SASSI Modeling*, NRC training slides; Mertz, G.E. et al. (2011). SASSI
  subtraction method effects at various DOE projects; US DOE (2011). OE-3 2011-02.
* Bathe, K.J., Dvorkin, E.N. (1985). A four-node plate bending element based on Mindlin/Reissner plate
  theory and a mixed interpolation. *IJNME* 21; Batoz, J.L., Bathe, K.J., Ho, L.W. (1980). A study of
  three-node triangular plate bending elements. *IJNME* 15.
* Gasparini, D.A., Vanmarcke, E.H. (1976). *SIMQKE*. MIT; Al Atik, L., Abrahamson, N. (2010). An
  improved method for nonstationary spectral matching. *Earthquake Spectra* 26(3).
* Lee, P.S., Bathe, K.J. (2004). Development of MITC isotropic triangular shell finite elements.
  *Computers & Structures* 82; Lyly, M., Stenberg, R., Vihinen, T. (1993). A stable bilinear element for
  the Reissner-Mindlin plate model. *Comput. Methods Appl. Mech. Eng.* 110 (TSHELL, section 13.2).
* Cuthill, E., McKee, J. (1969). Reducing the bandwidth of sparse symmetric matrices. *Proc. 24th ACM
  National Conference*; George, A., Liu, J.W.H. (1981). *Computer Solution of Large Sparse Positive
  Definite Systems*. Prentice-Hall (node optimizer, section 13.3).
* Luco, J.E., Wong, H.L. (1986). Response of a rigid foundation to a spatially random ground motion.
  *EESD* 14, 891-908; Abrahamson, N.A. (2006). *Spatial Coherency Models for Soil-Structure
  Interaction*. EPRI 1012968 (NRC ADAMS ML060400268); Abrahamson, N.A. (2007). *Effects of Spatial
  Incoherence on Seismic Ground Motions*. EPRI 1015110 (NRC ADAMS ML090960409) (coherency models,
  section 14; the transcriptions and page references are in `sassi/data/coherency/*.json`).
* Jakub, M., Roesset, J.M. (1977). MIT Research Report R77-36 (rigid strip on a layer, R2 B.4; section
  15).
* Cheng, F.Y., Mertz, G.E. (1989). *A computer program for inelastic analysis of 3-dimensional
  reinforced-concrete and steel seismic buildings*. Civil Engineering Study 89-31, University of
  Missouri-Rolla; Takeda, T., Sozen, M.A., Nielsen, N.N. (1970). Reinforced concrete response to
  simulated earthquakes. *J. Struct. Div. ASCE* 96(ST12); Otani, S. (1974). Report UILU-ENG-74-2029,
  University of Illinois; the shear-capacity equations of ACI 318-08, Wood (1990), Barda et al. (1977)
  and Gulec and Whittaker (2009) (Option NON, section 17).
* Masing, G. (1926); Pyke, R. (1979); Kramer, S.L. (1996). *Geotechnical Earthquake Engineering*,
  §6.4 (Masing and extended Masing rules, sections 17 and 18).
* Matasovic, N., Vucetic, M. (1993). Cyclic characterization of liquefiable sands. *J. Geotech. Eng.*
  119(11); Hashash, Y.M.A., Park, D. (2001). Non-linear one-dimensional seismic ground motion
  propagation in the Mississippi embayment. *Engineering Geology* 62; Phillips, C., Hashash, Y.M.A.
  (2009). Damping formulation for nonlinear 1D site response analyses. *Soil Dyn. Earthq. Eng.* 29(7);
  Joyner, W.B., Chen, A.T.F. (1975). Calculation of nonlinear ground response in earthquakes. *BSSA*
  65(5) (SOIL-NON, section 18).
* Westergaard, H.M. (1933). Water pressures on dams during earthquakes. *Trans. ASCE* 98; Housner, G.W.
  (1963). The dynamic behavior of water tanks. *BSSA* 53(2) (water, section 19).
* USNRC (2014). SRP 3.7.1 Rev. 4; RG 1.60 Rev. 2.
* Ghiocel Predictive Technologies (2017). *ACS SASSI Version 3 User Manual* (nomenclature only).
