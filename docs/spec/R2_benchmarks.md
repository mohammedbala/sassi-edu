# R2 — Verification Benchmarks for the Educational SASSI-type SSI Program

**Purpose.** This is a catalogue of reference solutions for checking each module of the educational
SASSI re-implementation: SITE/SOIL (site response), POINT (layered Green functions), the impedance
and SSI solution (ANALYS/HOUSE), MOTION/STRESS (transfer functions and response spectra), EQUAKE
(spectrum-compatible motions) and the HOUSE element library. Every item states:

* **Type**
  * **E** = exact closed form.
  * **P** = published numerical values (tabulated).
  * **A** = published approximate formula (the stated accuracy is the authors').
  * **D** = derived here. This means values evaluated from an E or A formula, or a discrete-model
    exact result. They are reproducible with the script logic given and are *not* independent
    published data.
* **Assumptions** and the **recommended tolerance** for an automated test.
* **Source**, with a URL wherever an open copy exists.

Where a value could not be confirmed from an accessible primary source, it is marked **[VERIFY]**.

---

## 0. Global conventions (apply to every benchmark below)

| Item | Convention used in this document |
|---|---|
| Time factor | `e^{+iωt}` (SASSI/SHAKE convention). A displacement lagging the force has a **negative** imaginary part. Outgoing waves are `e^{-ikr}`. |
| Complex modulus, SASSI/SHAKE91 form ("Model 2") | `G* = G(1 − 2β² + 2iβ√(1−β²))`, so `|G*| = G`. The complex wave speed is `V* = V(√(1−β²) + iβ)`. This is the form SHAKE91 uses: see the EERA manual §2.2.2, Eq. 18 ([Bardet et al. 2000](http://www.ce.memphis.edu/7137/PDFs/EERA2/EERAManual.pdf)). |
| Complex modulus, Kelvin–Voigt-like form ("Model 1") | `G* = G(1 + 2iβ)`, so `V* = V√(1+2iβ)`. Keep it selectable for benchmarking only. |
| Dimensionless frequency | `a0 = ωR/Vs` (circle of radius R), or `a0 = ωB/Vs` (rectangle of half-width B). |
| Impedance normalisation | `K(ω) = K0 [k(a0) + i a0 c(a0)]`. Here K0 is the static stiffness, k is the dynamic stiffness coefficient and c the damping coefficient. The dashpot is `C = K0 R c / Vs`. The radiation damping ratio is `β_rad = a0 c /(2k)`. |
| Rectangular foundation | Plan 2L × 2B with L ≥ B (NIST GCR 12-917-21 convention). |
| Units in published data | SHAKE91 example: ft, kcf, ft/s, g; spectra in cm and cm/s (g = 981 cm/s²). |

**Why the damping form matters (E).** For an SDOF with hysteretic damping β, the peak of the
total-acceleration transfer function is:

* SASSI form: `|H|max = 1/(2β√(1−β²))`, at `ω = ω0√(1−2β²)`.
* `(1+2iβ)` form: `|H|max = √(1+4β²)/(2β)`, at `ω = ω0`.

Ostadan, Deng & Roesset (2004, Table 1) analysed a fixed-base 4 Hz SDOF with SASSI2000. They
recovered damping `1/(2Umax)` = **5.0, 9.9, 14.7, 19.5 %** for input β = 5, 10, 15, 20 %. The SASSI
form predicts **4.99, 9.95, 14.83, 19.60 %**. The (1+2iβ) form predicts 4.98, 9.81, 14.37, 18.57 %.
So the SASSI form is confirmed as the SASSI2000 convention. ([Ostadan et al. 2004, 3rd UJNR SSI workshop](https://www.pwri.go.jp/eng/ujnr/tc/a/ssi_w3/Contributions/Ostadan.pdf).)

| β | SASSI-form peak (D) | ω_peak/ω0 | (1+2iβ)-form peak (D) |
|---|---|---|---|
| 0.02 | 25.0050 | 0.99960 | 25.0200 |
| 0.05 | 10.0125 | 0.99750 | 10.0499 |
| 0.10 | 5.0252 | 0.98995 | 5.0990 |
| 0.15 | 3.3715 | 0.97724 | 3.4801 |
| 0.20 | 2.5516 | 0.95917 | 2.6926 |

Tolerance: 1e-6 relative (exact algebra) for a fine frequency grid. Otherwise the error is set by
the frequency step.

---

## A. One-dimensional site response (vertically propagating SH) — SITE/SOIL

### A.1 Exact transfer functions (E)

Layer m has thickness h_m, density ρ_m, complex speed V*_m and wavenumber `k*_m = ω/V*_m`. The
local coordinate z_m runs downward from the top of the layer. The displacement is

`u_m = [A_m e^{i k*_m z_m} + B_m e^{−i k*_m z_m}] e^{iωt}`

with A_m the up-going (incident) wave and B_m the down-going (reflected) wave.

* Free surface: `A_1 = B_1`.
* Recursion (Kramer 1996, Eqs. 7.33–7.34; SHAKE), reproduced in [BNL N6112-051208 (Costantino
  2009) Eqs. 8–11](https://www.nrc.gov/docs/ML0919/ML091980384.pdf):

```
A_{m+1} = ½ A_m (1+α*_m) e^{+i k*_m h_m} + ½ B_m (1−α*_m) e^{−i k*_m h_m}
B_{m+1} = ½ A_m (1−α*_m) e^{+i k*_m h_m} + ½ B_m (1+α*_m) e^{−i k*_m h_m}
α*_m = ρ_m V*_m / (ρ_{m+1} V*_{m+1})            (complex impedance ratio)
```

Motion definitions (the source of most "outcrop" errors):

* **within** motion at the top of layer m = `A_m + B_m`.
* **outcrop** motion of layer m = `2A_m`.
* surface motion = `2A_1`.

Closed forms for a uniform soil layer (thickness H, speed Vs, damping ξ):

| Case | Transfer function |
|---|---|
| Rigid base: surface / base (within) | `F1(ω) = 1 / cos(k* H)` |
| Elastic half-space (ρr, V*r): surface / rock **outcrop** | `F2(ω) = 1 / [cos(k*_s H) + i α* sin(k*_s H)]`, with `α* = ρs V*s /(ρr V*r)` |
| Elastic half-space: surface / rock **within** (top of rock) | `1/cos(k*_s H)` (same as F1) |
| Natural frequencies (undamped) | `f_n = (2n−1) Vs /(4H)` |
| Peak near f_n, small ξ (A) | rigid base: `|F1| ≈ 1/[ξ(2n−1)π/2]`; half-space: `|F2| ≈ 1/[α + ξ(2n−1)π/2]` (Kramer 1996) |

Important behaviour checks:

* Below f1 the rigid-base ratio surface/within-base is exactly `1/|cos|`.
* The within motion at depth z is the surface motion times `cos(k* z)`. Under SASSI's
  "within/outcrop" input options this must be reproduced exactly.
* An embedded control point gives zero motion (infinite amplification) at the column frequencies
  `f = (2n−1)V/(4z)` when the soil is undamped. Elsabee & Morray (1977) state the same for the
  embedded transfer function. ([MIT R77-33](https://nehrpsearch.nist.gov/static/files/NSF/PB286493.pdf))
* With an outcrop definition that keeps the overlying layers, the amplitude at an interface is
  `2A_m`. This is **not** a true free-surface outcrop: see BNL N6112-051208 §2. A test should check
  both definitions.

### A.2 Numerical check values (D, from A.1)

Layer: H = 30 m, Vs = 200 m/s, ρ = 2.0 t/m³, ξ = 5 %. Rock: Vr = 1000 m/s, ρr = 2.2 t/m³, ξr = 1 %.
This gives f1 = 1.6667 Hz and α0 = ρVs/(ρrVr) = 0.18182.

| f (Hz) | SASSI form \|surf/base-within\| | SASSI form \|surf/outcrop\| | KV form \|surf/base-within\| | KV form \|surf/outcrop\| |
|---|---|---|---|---|
| 0.5 | 1.1216 | 1.1150 | 1.1209 | 1.1143 |
| 1.0 | 1.6931 | 1.6207 | 1.6878 | 1.6161 |
| 1.5 | 5.7710 | 3.4116 | 5.6735 | 3.3916 |
| 1.6667 | 12.7152 | 3.8327 | 12.7633 | 3.8344 |
| 2.0 | 3.1156 | 2.4115 | 3.1590 | 2.4313 |
| 3.0 | 1.0411 | 1.0087 | 1.0436 | 1.0110 |
| 5.0 (= 3 f1) | 4.2038 | 2.3554 | 4.2202 | 2.3610 |
| 8.3333 (= 5 f1) | 2.4815 | 1.6702 | 2.4918 | 1.6756 |

The Kramer approximations at f1 give 12.73 (rigid) and 3.841 (half-space), in agreement.
Tolerance: 1e-8 relative for the closed form. A multi-sublayer model of the same uniform layer must
reproduce these values to round-off, because sub-dividing a uniform layer is exact in the recursion.

### A.3 Discrete soil columns (D, exact for the discrete model)

1-D linear bar/shear elements with element size h and wavenumber k give a discrete-to-exact
frequency ratio `ω_h/ω`:

* consistent mass: `√[6(1−cos kh)/(2+cos kh)]/(kh)`.
* lumped mass: `√[2(1−cos kh)]/(kh)`.
* **50 % lumped + 50 % consistent** (the ACS SASSI solid-element default; see 01_capabilities §1.5):
  `√[12(1−cos kh)/(5+cos kh)]/(kh)`. This equals `1 − (kh)⁴/480 + …`, i.e. 4th-order accurate.

| Elements per wavelength N | kh | consistent | lumped | 50/50 |
|---|---|---|---|---|
| 4 | 1.5708 | 1.10266 | 0.90032 | 0.98625 |
| 5 (SASSI λ/5 rule) | 1.2566 | 1.06632 | 0.93549 | 0.99451 |
| 6 | 1.0472 | 1.04607 | 0.95493 | 0.99739 |
| 8 | 0.7854 | 1.02586 | 0.97450 | 0.99919 |
| 10 | 0.6283 | 1.01652 | 0.98363 | 0.99967 |
| 20 | 0.3142 | 1.00412 | 0.99589 | 0.99998 |

For a fixed-free column of Ne equal elements (ω1 exact = πVs/2H), the eigen-analysis gives:

| Ne | consistent ω1/ω1,exact | lumped | 50/50 | 50/50 ω2/ω2,exact |
|---|---|---|---|---|
| 2 | 1.02586 | 0.97450 | 0.99919 | 0.92712 |
| 4 | 1.00644 | 0.99359 | 0.99995 | 0.99578 |
| 8 | 1.00161 | 0.99839 | 1.00000 | 0.99975 |

This is the most useful single regression test for HOUSE solid/plane elements against SITE, and for
an ANSYS cross-check model (SOLID185 with rollers on a rigid base). Tolerance: 1e-6 for the
eigenvalues of the discrete model.

---

## B. Static stiffness of rigid foundations — impedance low-frequency limit

### B.1 Rigid circular disk on a homogeneous elastic half-space (E, relaxed contact)

| DOF | K0 | ν = 0.25 | ν = 1/3 | ν = 0.40 | ν = 0.45 |
|---|---|---|---|---|---|
| Vertical | `4GR/(1−ν)` | 5.3333 GR | 6.0000 GR | 6.6667 GR | 7.2727 GR |
| Horizontal | `8GR/(2−ν)` | 4.5714 GR | 4.8000 GR | 5.0000 GR | 5.1613 GR |
| Rocking | `8GR³/[3(1−ν)]` | 3.5556 GR³ | 4.0000 GR³ | 4.4444 GR³ | 4.8485 GR³ |
| Torsion | `16GR³/3` | 5.3333 GR³ | 5.3333 GR³ | 5.3333 GR³ | 5.3333 GR³ |

Accuracy notes:

* Vertical and torsion are exact. Vertical is the Boussinesq rigid punch. Torsion is
  Reissner–Sagoci.
* Rocking is exact for **relaxed** (frictionless) contact.
* Horizontal `8GR/(2−ν)` is the relaxed approximation (Bycroft/Veletsos). The welded solution
  couples horizontal translation and rocking weakly, and its stiffness differs by a few percent.
  Richart, Hall & Woods (1970) use `32(1−ν)GR/(7−8ν)` (4.92 GR at ν = 1/3).
* A FE-based (welded) SASSI result should therefore be within about 3–5 % of the table, after mesh
  extrapolation.

Luco (1976) and Luco et al. (1978, Table 3.2) give torsion for a disk (h = 0) as **5.33 GR³**
(integral equations) and **5.53** (FE). This shows the typical FE overestimate.
([PB296617](https://nehrpsearch.nist.gov/static/files/NSF/PB296617.pdf))

### B.2 Rectangular rigid footings on a half-space (A)

Source: NIST GCR 12-917-21, Table 2-2a ([pdf](https://www.nehrp.gov/pdf/nistgcr12-917-21.pdf)).
It quotes Pais & Kausel (1988, *Soil Dyn. Earthq. Eng.* 7(4):213–227) and Gazetas (1991). Plan is
2L × 2B, L ≥ B; x is along L, y is along B, z is vertical.

**Pais & Kausel (1988):**
```
K_z  = GB/(1−ν) [3.1 (L/B)^0.75 + 1.6]
K_y  = GB/(2−ν) [6.8 (L/B)^0.65 + 0.8 (L/B) + 1.6]
K_x  = GB/(2−ν) [6.8 (L/B)^0.65 + 2.4]
K_zz = GB³ [4.25 (L/B)^2.45 + 4.06]
K_yy = GB³/(1−ν) [3.73 (L/B)^2.4 + 0.27]
K_xx = GB³/(1−ν) [3.2 (L/B) + 0.8]
```
**Gazetas (1991)** (with I_x, I_y the area moments, J_t = I_x + I_y):
```
K_z  = 2GL/(1−ν) [0.73 + 1.54 (B/L)^0.75]
K_y  = 2GL/(2−ν) [2 + 2.5 (B/L)^0.85]
K_x  = K_y − 0.2/(0.75−ν) · GL (1 − B/L)
K_zz = G J_t^0.75 [4 + 11 (1 − B/L)^10]
K_yy = G/(1−ν) I_y^0.75 [3 (L/B)^0.15]
K_xx = G/(1−ν) I_x^0.75 (L/B)^0.25 [2.4 + 0.5 (B/L)]
```

**Check values for a square, L = B, ν = 1/3 (D).** Units are G·B for translation and G·B³ for
rotation.

| | K_z | K_y | K_x | K_zz | K_yy | K_xx |
|---|---|---|---|---|---|---|
| Pais & Kausel | 7.0500 | 5.5200 | 5.5200 | 8.3100 | 6.0000 | 6.0000 |
| Gazetas | 6.8100 | 5.4000 | 5.4000 | 8.3471 | 5.5836 | 5.3975 |
| Circle of equal area/inertia | 6.7703 (R = 1.1284B) | 5.4162 | 5.4162 | 7.9326 (R = 1.1415B) | 5.9490 (R = 1.1415B) | 5.9490 |

Built-in checks:

* Pais & Kausel is exactly symmetric at L = B (K_x = K_y, K_xx = K_yy). This is a good unit test of
  the implementation.
* Gazetas is not exactly symmetric at L = B (the rocking terms differ by about 3 %), because the
  formulas are fitted.
* A square is 4–5 % stiffer than the equal-area circle.

Accuracy: the fits are claimed by their authors to be within a few percent of rigorous solutions
over the fitted L/B range. This was not re-verified here, so use a ±5 % tolerance.

### B.3 Embedment correction factors (A)

`K_emb = η K_sur` (NIST GCR 12-917-21, Table 2-2b). D is the embedment depth, d_w the effective
side-wall contact height, z_w the depth to the centroid of the wall contact, and A_w the wall
contact area.

**Pais & Kausel (1988):**
```
η_z  = 1 + (0.25 + 0.25/(L/B)) (D/B)^0.8
η_y  = 1 + (0.33 + 1.34/(1 + L/B)) (D/B)^0.8 ;   η_x ≈ η_y
η_zz = 1 + (1.3 + 1.32/(L/B)) (D/B)^0.9
η_yy = 1 + D/B + (1.6/(0.35 + (L/B)^4)) (D/B)²
η_xx = 1 + D/B + (1.6/(0.35 + L/B)) (D/B)²
coupling (as printed): K_emb,rx = (D/3) K_emb,x ;  K_emb,ry = (D/3) K_emb,y   [VERIFY axis pairing]
```
**Gazetas (1991):**
```
η_z  = [1 + D/(21B) (1 + 1.3 B/L)] [1 + 0.2 (A_w/(4BL))^(2/3)]
η_y  = [1 + 0.15 √(D/B)] [1 + 0.52 (z_w A_w/(B L²))^0.4]       (η_x: same form, A_w term changes)
η_zz = 1 + 1.4 (1 + B/L) (d_w/B)^0.9
η_yy = 1 + 0.92 (d_w/B)^0.6 [1.5 + (d_w/D)^1.9 (B/L)^−0.6]
η_xx = 1 + 1.26 (d_w/B) [1 + (d_w/B)(d_w/D)^−0.2 √(B/L)]
```

**Published static torsional stiffness of embedded rigid cylinders (P).** Values are K_t/(G a³),
from Luco et al. (1978, Table 3.2; integral equations of Luco 1976 vs FE):

| h = E/a | 0 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|
| Integral equation (Luco 1976) | 5.33 | 13.23 | 19.89 | 32.75 |
| Finite element (UCSD) | 5.53 | 13.17 | 19.65 | 32.03 |
| Formula `16/3 (1 + 8E/(3a))` (A, Kausel-type) | 5.33 | 12.44 | 19.56 | 33.78 |

### B.4 Foundations on a soil layer over a rigid base (A)

Circular footing of radius R on a layer of thickness H. Kausel (1974) is quoted by
[Elsabee & Morray 1977, MIT R77-33, p. 18–19](https://nehrpsearch.nist.gov/static/files/NSF/PB286493.pdf).
Valid for R/H ≲ 1/2.
```
K_xx = 8GR/(2−ν) (1 + R/(2H))
K_φφ = 8GR³/[3(1−ν)] (1 + R/(6H))
K_xφ = −0.03 R K_xx          (small coupling)
```
Embedded to depth E in the layer (Elsabee & Morray 1977, valid R/H ≤ 1/2, E/R ≤ 1, E/H ≤ 1/2):
```
K_xx = 8GR/(2−ν) (1 + R/(2H)) (1 + 2E/(3R)) (1 + 5E/(4H))
K_φφ = 8GR³/[3(1−ν)] (1 + R/(6H)) (1 + 2E/R) (1 + 0.7E/H)
K_xφ = (0.4 E/R − 0.03) R K_xx
```
Vertical and torsion on a layer (Kausel 1974; Kausel & Ushijima 1979). These are taken from
secondary sources and not seen in a primary source here, so they are marked **[VERIFY]**:
```
K_v = 4GR/(1−ν) (1 + 1.28 R/H)          (some reproductions write 1.3)
K_t ≈ 16GR³/3                           (layer effect negligible for H/R ≳ 2)
embedded: K_v = 4GR/(1−ν)(1+1.28R/H)(1+E/(2R))[1+(0.85−0.28E/R)(E/H)/(1−E/H)] ;  K_t = 16GR³/3 (1+8E/(3R))
```
Square footing (half-width B) on a layer (Gonzalez 1977) and a plane-strain strip on a layer (per
unit length). Both are from [Jakub & Roesset 1977, MIT R77-36, pp. 12–14](https://nehrpsearch.nist.gov/static/files/NSF/PB286504.pdf),
valid 1/8 ≤ B/H ≤ 1/2:
```
square: K_x = 4.2 GB/(2−ν) (1 + 2B/H)        K_φ = 3.24 GB³/(1−ν) (1 + 0.2 B/H)
strip : K_x = 2.1 G/(2−ν)  (1 + 2B/H)        K_φ = 1.62 GB²/(1−ν) (1 + 0.2 B/H)
strip, Poisson-specific (Table 1/2): K_x/G = 0.904(1+2.5B/H) [ν=0], 1.017(1+2.37B/H) [0.15],
       1.175(1+2.15B/H) [0.30], 1.419(1+1.95B/H) [0.45];
       K_φ/GB² = 1.891(1+0.17B/H) [0], 2.094(1+0.17B/H) [0.15], 2.394(1+0.17B/H) [0.30], 2.907(1+0.24B/H) [0.45]
```
A strip on a half-space has zero static horizontal stiffness. This is a useful 2-D sanity check:
the computed K_x must go to 0 as H → ∞.

---

## C. Dynamic impedance of rigid massless foundations

### C.1 Exact high-frequency asymptotes (E; derived from plane-wave radiation, confirmed by published values)

As a0 → ∞ the dashpots tend to the plane-wave values: shear traction ρVs·velocity and normal
traction ρVp·velocity.

**Surface disk, normalised `c(∞)`:**

| | c(∞) formula | ν = 0 | ν = 1/3 | ν = 0.5 |
|---|---|---|---|---|
| horizontal | `π(2−ν)/8` | 0.785 | 0.654 | 0.589 |
| vertical | `π(1−ν)/4 · Vp/Vs` | 1.111 | 1.047 | ∞ (use trapped mass) |
| rocking | `(3π/32)(1−ν) Vp/Vs` | 0.417 | 0.393 | ∞ |
| torsion | `3π/32` | 0.295 | 0.295 | 0.295 |

These agree to within 2 % with the Veletsos–Verbic limits in C.2: α1 = 0.775/0.65/0.60,
γ4 + γ1γ2 = 1.10/1.03, and β1β2 = 0.42/0.40.

**Embedded rigid cylinder (radius a, depth h·a) in a half-space, ν = 1/4 (α/β = √3).** Dashpots are
normalised by ρVs a² (translation) or ρVs a⁴ (rotation about the base centre). The formulas are:

* `C_HH = π[1 + h(1+√3)]`
* `C_VV = √3π + 2πh`
* `C_TT = π/2 + 2πh`
* `C_MM = √3π/4 + π(1+√3)h³/3 + πh`

| h | C_HH (D / P) | C_VV (D / P) | C_TT (D / P) | C_MM (D / P) |
|---|---|---|---|---|
| 0 | 3.142 / 3.14 | 5.441 / 5.44 | 1.571 / 1.57 | 1.360 / 1.36 |
| 0.5 | 7.433 / 7.43 | 8.583 / 8.58 | 4.712 / 4.71 | 3.289 / 3.29 |
| 1.0 | 11.725 / 11.72 | 11.725 / 11.72 | 7.854 / 7.85 | 7.363 / 7.36 |
| 2.0 | 20.308 / 20.29 | 18.008 / 18.00 | 14.137 / 14.14 | 30.532 / 30.51 |

"P" is Luco et al. (1978) Table 3.3, "analytic a0 → ∞". The same table gives FE values at a0 = 6
(e.g. h = 1: 7.01, 7.40, 11.54, 11.49). [PB296617](https://nehrpsearch.nist.gov/static/files/NSF/PB296617.pdf).
Tolerance: computed `Im K /(ω)` at a0 ≳ 6–8 should be within about 5 % of these values.

### C.2 Veletsos & Verbic (1973) analytical approximations — rigid disk, relaxed contact (A)

Veletsos, A.S. & Verbic, B. (1973), *Earthquake Eng. Struct. Dyn.* 2(1):87–102. The formulas and
coefficients below are taken from the reproduction in Shye & Robinson (1980), *Dynamic
Soil-Structure Interaction*, UIUC SRS 484, Eqs. 4.46–4.57 and Table 4.1
([PB81144503](https://nehrpsearch.nist.gov/static/files/NSF/PB81144503.pdf)).

```
K* = K0 (k + i a0 c),  a0 = ωR/Vs
Horizontal: k = 1,                                     c = α1
Rocking   : k = 1 − β1 (β2 a0)²/[1+(β2 a0)²] − β3 a0², c = β1 β2 (β2 a0)²/[1+(β2 a0)²]
Vertical  : k = 1 − γ1 (γ2 a0)²/[1+(γ2 a0)²] − γ3 a0², c = γ4 + γ1 γ2 (γ2 a0)²/[1+(γ2 a0)²]
Mechanical analog (vertical): c1 = K0 γ1 γ2 R/Vs, m2 = K0 γ1 (γ2 R/Vs)², m3 = K0 γ3 (R/Vs)², c4 = K0 γ4 R/Vs
(rocking analog: replace γ by β)
```

| Coefficient | ν = 0 | ν = 1/3 | ν = 0.45 | ν = 0.5 |
|---|---|---|---|---|
| α1 | 0.775 | 0.65 | 0.60 | 0.60 |
| β1 | 0.8 | 0.8 | 0.8 | 0.8 |
| β2 | 0.525 | 0.5 | 0.45 | 0.4 |
| β3 | 0 | 0 | 0.023 | 0.027 |
| γ1 | 0.25 | 0.35 | — | 0 |
| γ2 | 1.0 | 0.8 | — | 0 |
| γ3 | 0 | 0 | — | 0.17 |
| γ4 | 0.85 | 0.75 | — | 0.85 |

The vertical coefficients for ν = 0.45 are not given in the reproduction. Get them from the
original paper or interpolate, flagged as an approximation.

**Evaluated values (D).** Horizontal: k_h = 1 and c_h = α1 at every a0.

| a0 | ν=0 k_r | ν=0 c_r | ν=0 k_v | ν=0 c_v | ν=1/3 k_r | ν=1/3 c_r | ν=1/3 k_v | ν=1/3 c_v | ν=0.5 k_r | ν=0.5 c_r | ν=0.5 k_v | ν=0.5 c_v |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 1.0000 | 0.0000 | 1.0000 | 0.8500 | 1.0000 | 0.0000 | 1.0000 | 0.7500 | 1.0000 | 0.0000 | 1.0000 | 0.8500 |
| 0.25 | 0.9865 | 0.0071 | 0.9853 | 0.8647 | 0.9877 | 0.0062 | 0.9865 | 0.7608 | 0.9904 | 0.0032 | 0.9894 | 0.8500 |
| 0.50 | 0.9484 | 0.0271 | 0.9500 | 0.9000 | 0.9529 | 0.0235 | 0.9517 | 0.7886 | 0.9625 | 0.0123 | 0.9575 | 0.8500 |
| 0.75 | 0.8926 | 0.0564 | 0.9100 | 0.9400 | 0.9014 | 0.0493 | 0.9074 | 0.8241 | 0.9188 | 0.0264 | 0.9044 | 0.8500 |
| 1.00 | 0.8271 | 0.0907 | 0.8750 | 0.9750 | 0.8400 | 0.0800 | 0.8634 | 0.8593 | 0.8627 | 0.0441 | 0.8300 | 0.8500 |
| 1.25 | 0.7592 | 0.1264 | 0.8476 | 1.0024 | 0.7753 | 0.1124 | 0.8250 | 0.8900 | 0.7978 | 0.0640 | 0.7344 | 0.8500 |
| 1.50 | 0.6938 | 0.1608 | 0.8269 | 1.0231 | 0.7120 | 0.1440 | 0.7934 | 0.9152 | 0.7275 | 0.0847 | 0.6175 | 0.8500 |
| 2.00 | 0.5805 | 0.2202 | 0.8000 | 1.0500 | 0.6000 | 0.2000 | 0.7483 | 0.9513 | 0.5798 | 0.1249 | 0.3200 | 0.8500 |

For ν = 0.45 (rocking only): k_r = 1, 0.9886, 0.9557, 0.9053, 0.8423, 0.7718, 0.6978, 0.5500 and
c_r = 0, 0.0045, 0.0173, 0.0368, 0.0606, 0.0865, 0.1127, 0.1611 at the same a0.

Accuracy:

* Veletsos & Verbic fitted these to the Veletsos–Wei (1971) and Luco–Westmann (1971) rigorous
  solutions (relaxed boundary, elastic half-space). The approximations are reliable to about 5–10 %
  for 0 ≤ a0 ≤ 1.5. Agreement degrades for a0 > 2 and for k_h, which is truly mildly
  frequency-dependent.
* **Recommended use:** test that an educational SASSI model of a rigid massless disk (fine
  interaction-node mesh, many thin layers plus a damped half-space) gives K/K0 within ±10 % of this
  table for a0 ≤ 1.5, and approaches the C.1 asymptotes.
* Viscoelastic soil: for a half-space with damping β, multiply the elastic impedance approximately
  by `(1+2iβ)` (correspondence principle at fixed a0). Veletsos & Verbic (1973) give the exact
  correspondence-principle form.

The rigorous tabulated solutions were not accessible online. Obtain them from the original papers
and add them as [P] benchmarks:

* Luco & Westmann (1971), *J. Eng. Mech. Div. ASCE* 97(EM5):1381–1395 (ν = 0, 1/4, 1/2).
* Veletsos & Wei (1971), *J. SMFD ASCE* 97(SM9):1227–1248 (ν = 0, 1/3, 1/2).
* Wong & Luco (1978, USC CE 78-15; and 1985, *Soil Dyn. Earthq. Eng.* 4(2):64–81, square
  foundations on layered media).

### C.3 Rectangular foundations: Pais & Kausel (1988) dynamic modifiers and radiation damping (A)

From NIST GCR 12-917-21 Table 2-3a. `k_j = K_j·α_j`, with damping ratio β_j added to the soil
hysteretic damping. a0 = ωB/Vs and ψ = √(2(1−ν)/(1−2ν)) ≤ 2.5.
```
α_z  = 1 − (0.4 + 0.2/(L/B)) a0² / [10/(1+3(L/B−1)) + a0²]
α_y  = α_x = 1
α_zz = 1 − (0.33 − 0.03√(L/B−1)) a0² / [0.8/(1+0.33(L/B−1)) + a0²]
α_yy = 1 − 0.55 a0² / [(0.6 + 1.4/(L/B)³) + a0²]
α_xx = 1 − (0.55 + 0.01√(L/B−1)) a0² / [(2.4 − 0.4/(L/B)³) + a0²]
β_z  = [4ψ(L/B)/(K_z,sur/GB)] [a0/(2α_z)]
β_y  = [4(L/B)/(K_y,sur/GB)] [a0/(2α_y)] ;  β_x = [4(L/B)/(K_x,sur/GB)] [a0/(2α_x)]
β_zz = [(4/3)((L/B)³+(L/B)) a0² / ((K_zz,sur/GB³)(1.4/(1+3(L/B−1)^0.7) + a0²))] [a0/(2α_zz)]
β_yy = [(4ψ/3)(L/B)³ a0² / ((K_yy,sur/GB³)(1.8/(1+1.75(L/B−1)) + a0²))] [a0/(2α_yy)]
β_xx = [(4ψ/3)(L/B) a0² / ((K_xx,sur/GB³)(2.2 − 0.4/(L/B)³ + a0²))] [a0/(2α_xx)]
```
Table 2-3b gives the embedded versions.

**Square, ν = 1/3 (ψ = 2), evaluated (D):**

| a0 | α_z | β_z | α_yy=α_xx | β_yy=β_xx | α_zz | β_zz | β_x=β_y |
|---|---|---|---|---|---|---|---|
| 0.25 | 0.9963 | 0.1424 | 0.9833 | 0.0019 | 0.9761 | 0.0018 | 0.0906 |
| 0.50 | 0.9854 | 0.2879 | 0.9389 | 0.0144 | 0.9214 | 0.0132 | 0.1812 |
| 1.00 | 0.9455 | 0.6001 | 0.8167 | 0.0972 | 0.8167 | 0.0819 | 0.3623 |
| 1.50 | 0.8898 | 0.9565 | 0.7088 | 0.2613 | 0.7566 | 0.1961 | 0.5435 |
| 2.00 | 0.8286 | 1.3695 | 0.6333 | 0.4840 | 0.7250 | 0.3279 | 0.7246 |

Consistency check: the equal-area circle (R = 1.128B, a0,circle = 1.128 a0) with Veletsos–Verbic at
a0 = 1 gives β_z ≈ 0.59. Pais & Kausel give 0.60.

### C.4 Gazetas (1991) dynamic coefficients — arbitrary basemat on a half-space (A)

Gazetas, *J. Geotech. Eng.* 117(9):1363–1381, Table 1
([open pdf](http://ssi.civil.ntua.gr/downloads/journals/1991-ASCE_FORMULAS%20AND%20CHARTS%20FOR%20IMPEDANCES%20OF%20SURFACE%20AND%20EMBEDDED%20FOUNDATIONS.pdf)).
Valid for 0 ≤ a0 = ωB/Vs ≤ 2.

* Dynamic stiffness coefficients:
  * k_x ≈ 1.
  * k_rx ≈ 1 − 0.20 a0.
  * k_ry ≈ 1 − 0.26 a0 (ν < 0.45); k_ry ≈ 1 − 0.26 a0 (L/B)^0.30 (ν ≈ 0.5).
  * k_t ≈ 1 − 0.14 a0.
  * k_z and k_y: charts (Fig. 2).
* Dashpots, with Lysmer's analog velocity `V_La = 3.4 Vs/[π(1−ν)]`:
  * `C_z = ρ V_La A_b c̃_z`
  * `C_y = ρ Vs A_b c̃_y`
  * `C_x = ρ Vs A_b`
  * `C_rx = ρ V_La I_bx c̃_rx`
  * `C_ry = ρ V_La I_by c̃_ry`
  * `C_t = ρ Vs J_t c̃_t`
* **Worked example (P):** basemat of Fig. 4a, G = 120 MPa, ρ = 1.85 Mg/m³, ν = 0.40, f = 18 Hz,
  a0 = 1.23, χ = A_b/4L² = 0.26, L/B = 3.2:
  * K_z(static) = 4.13×10⁶ kN/m; k_z ≈ 0.92; K_z = 3.8×10⁶ kN/m.
  * Radiation C_z = 56.9×10³ kN·s/m; total C_z = 60×10³ kN·s/m (β = 5 %).
  * K_y ≈ 3.9×10⁶ and K_x = 3.0×10⁶ kN/m.
  * K_rx(static) ≈ 25×10⁶ kN·m; k_rx ≈ 0.754.
  * K_t(static) ≈ 109×10⁶ kN·m; k_t ≈ 0.83.

  This is a regression test for a "Gazetas formula" helper utility.

### C.5 Lumped-parameter analogs (A; classical)

From Richart, Hall & Woods (1970), *Vibrations of Soils and Foundations*, ch. 10; Lysmer's vertical
analog and Hall's analog. Rigid circular footing of mass m, radius R:

| Mode | Spring | Dashpot | Mass/inertia ratio | Damping ratio |
|---|---|---|---|---|
| Vertical | 4GR/(1−ν) | 3.4R²√(ρG)/(1−ν) | B_z = (1−ν)m/(4ρR³) | D_z = 0.425/√B_z |
| Horizontal | 32(1−ν)GR/(7−8ν) | 18.4(1−ν)R²√(ρG)/(7−8ν) | B_x = (7−8ν)m/[32(1−ν)ρR³] | D_x = 0.288/√B_x |
| Rocking | 8GR³/[3(1−ν)] | 0.8R⁴√(ρG)/[(1−ν)(1+B_ψ)] | B_ψ = 3(1−ν)I_ψ/(8ρR⁵) | D_ψ = 0.15/[(1+B_ψ)√B_ψ] |
| Torsion | 16GR³/3 | (frequency-dependent) | B_θ = I_θ/(ρR⁵) | D_θ = 0.50/(1+2B_θ) |

Check: the Lysmer vertical dashpot gives c_v = 3.4/4 = 0.85, which equals Veletsos–Verbic γ4 at
ν = 0 and ν = 0.5.

### C.6 Wolf cone models (A; static stiffness exact by construction)

Wolf (1994), *Foundation Vibration Analysis Using Simple Physical Models*, Prentice Hall; Wolf &
Deeks (2004), *Foundation Vibration Analysis: a Strength-of-Materials Approach*, Elsevier; Meek &
Wolf (1992), *J. Geotech. Eng.* 118(5). Surface disk of radius r0, area A0 = πr0², moments
I0 = πr0⁴/4 (rocking) or πr0⁴/2 (torsion).

| Mode | Aspect ratio z0/r0 | Wave speed c | Static stiffness | Trapped mass (ν > 1/3) |
|---|---|---|---|---|
| Horizontal | (π/8)(2−ν) | Vs | ρc²A0/z0 = 8Gr0/(2−ν) | — |
| Vertical | (π/4)(1−ν)(c/Vs)² | Vp (ν ≤ 1/3); 2Vs (1/3 < ν ≤ 1/2) | ρc²A0/z0 = 4Gr0/(1−ν) | ΔM = 2.4(ν−1/3)ρA0r0 |
| Rocking | (9π/32)(1−ν)(c/Vs)² | same as vertical | 3ρc²I0/z0 = 8Gr0³/[3(1−ν)] | ΔM_θ = 1.2(ν−1/3)ρI0r0 |
| Torsion | 9π/32 | Vs | 3ρc²I0/z0 = 16Gr0³/3 | — |

Dynamic stiffness (with `b0 = ωz0/c`):

* Translational: `S(ω) = K + iωρcA0 − ω²ΔM`.
* Rotational: `S(ω) = K_θ[1 − (1/3) b0²/(1+b0²) + i (1/3) b0³/(1+b0²)] − ω²ΔM_θ`. The discrete form
  is a spring K_θ, a dashpot C_θ = ρcI0, plus an internal rotational mass M_θ = ρI0z0 coupled
  through −C_θ.

The static values are exact by construction; this was checked here symbolically. The dynamic values
agree with the rigorous solution to within about 20 % (Wolf 1994). The dynamic-stiffness algebra
above should be confirmed against Wolf 1994, Table 2-3A **[VERIFY]**.

### C.7 SASSI-specific published impedance benchmark (P, graphical)

Tajirian & Tabatabaie (1985), "Vibration Analysis of Foundations on Layered Media", ASCE Convention,
Detroit (Bechtel). This paper verifies the **flexible volume method in SASSI** against LUCON (Luco's
Green-function program).

* Problem: rigid massless circular disk, **r = 0.5**, on a uniform layer **H = 3** (H/r = 6) over a
  **rigid base**.
* Soil: Vs = 1, Vp = 2 (ν = 1/3), ρ = 1, β = **5 % and 15 %**.
* SASSI model: 24 sublayers of 0.125.
* Output: vertical, horizontal, rocking and torsional compliances vs A0 = ωr/Vs from 0 to 3. Results
  are in Figs. 4–11, graphical only.
* Agreement: SASSI and LUCON agree closely, except the real part of the 5 % vertical compliance for
  A0 < 0.2.

Exact features to test:

* Vertical layer resonances `f = (2n−1)Vp/(4H)` = 0.167, 0.500, 0.833, i.e. **A0 = 0.524, 1.571,
  2.618**.
* Horizontal resonances `(2n−1)Vs/(4H)` = 0.083, 0.250, 0.417, i.e. **A0 = 0.262, 0.785, 1.309**.
* Below the first layer frequency the radiation damping is nil, so Im K comes only from material
  damping. Elsabee & Morray (1977) give approximate transition formulas.
* Rocking and torsion are almost equal to the half-space values, because H/r > 3.

The half-space reference curves in the figures are consistent with Veletsos–Verbic ν = 1/3. For
example, at A0 = 1.2 the normalised vertical compliance is 0.46 − 0.58i.

This is the recommended **first end-to-end benchmark** for the educational code. The input file is
trivial to build. Digitise Figs. 4–11 of the open copy (Bechtel/MTR scan) to get numeric targets
[action item].

Related (P): Ostadan, Deng & Roesset (2004) used SASSI2000 for a surface rigid circular foundation
on (Case 2) a 5 %-damped half-space and (Case 3) a layer with H/R = 3 on a rigid base, with masses
tuned to 2–10 Hz. Effective system damping by the decay method:

* Case 2: **10.3, 15.2, 21.3, 27.4, 30.7 %**.
* Case 3: **5.1, 5.4, 18.1, 28.8, 39.6 %**, at 2, 4, 6, 8, 10 Hz.

Geometry details are in their Fig. 1 [VERIFY radius/Vs before use].

---

## D. Green functions (POINT / layered-soil flexibility)

### D.1 Static point loads (E)

Notation: μ = G, R = √(r²+z²), z positive into the solid.

**Kelvin (full space)**, with γ_i = x_i/R:

`u_i = P_j/(16πμ(1−ν)R) [(3−4ν)δ_ij + γ_iγ_j]`.

**Boussinesq** (vertical P at the surface, interior point):

* `u_z = P/(4πμR) [2(1−ν) + z²/R²]`
* `u_r = P/(4πμR) [rz/R² − (1−2ν) r/(R+z)]`

At the surface: `u_z = P(1−ν)/(2πμr)` and `u_r = −P(1−2ν)/(4πμr)`. Surface points move toward the
load.

**Cerruti** (horizontal Q along x at the surface, surface point):

* `u_x = Q/(2πμr)[(1−ν) + νx²/r²]`
* `u_y = Qνxy/(2πμr³)`
* `u_z = +Q(1−2ν)x/(4πμr²)` (into the solid for x > 0)

The sign of u_z is fixed by Betti reciprocity with the Boussinesq u_r. Test that the code's
flexibility matrix is symmetric, i.e. reciprocity G_zx = G_xz.

In cylindrical components, with θ measured from the load direction: `u_r = Q cosθ/(2πμr)` and
`u_θ = −Q(1−ν) sinθ/(2πμr)`.

Flexible uniform circular load of radius a and pressure p (classical, E):

* centre settlement `w0 = 2pa(1−ν²)/E = pa(1−ν)/μ`;
* edge settlement `(2/π)w0`;
* average settlement `16pa(1−ν²)/(3πE) = 0.849 w0`.

The disk-load option of a Green-function module can be tested against these. Rigid disk: `w = P(1−ν)/(4μa)`, which is exact
and consistent with B.1.

### D.2 Stokes dynamic full-space Green function (E; static limit verified here)

```
G_ij(r,ω) = 1/(4πμr) [ ψ δ_ij − χ γ_i γ_j ]
ψ = (1 − i/x_s − 1/x_s²) e^{−i x_s} − (Vs/Vp)² (−i/x_p − 1/x_p²) e^{−i x_p}
χ = (1 − 3i/x_s − 3/x_s²) e^{−i x_s} − (Vs/Vp)² (1 − 3i/x_p − 3/x_p²) e^{−i x_p}
x_s = ωr/Vs,  x_p = ωr/Vp   (e^{+iωt})
```
As ω → 0, `ψ → (3−4ν)/[4(1−ν)]` and `χ → −1/[4(1−ν)]`, i.e. Kelvin. This was checked numerically
here: ν = 0.25 and x_s = 10⁻³ give ψ = 0.666666 − 0.00073i and χ = −0.333333. Follows Kausel (2006),
*Fundamental Solutions in Elastodynamics*, CUP.

### D.3 Dynamic surface Green functions of a uniform half-space (P — tabulated)

Source: Apsel (1979), PhD dissertation, UCSD, Table 5.1
([PB81178204](https://nehrpsearch.nist.gov/static/files/NSF/PB81178204.pdf)). It reproduces
**Wong (1975)**, Caltech EERL 75-01, exact contour integration, perfectly elastic, **ν = 0.33**.

Setup: harmonic point force at the free surface, receiver at the surface, `r0 = ωr/Vs`.

Normalisation, deduced here: each entry is `μ r u / P` (Re, Im). This is confirmed exactly by the
static limits at r0 = 0:

* Boussinesq: −(1−2ν)/(4π) = −0.0271 and (1−ν)/(2π) = 0.1066.
* Cerruti: 1/(2π) = 0.1592 and −(1−ν)/(2π) = −0.1066.

Columns:

* **U_r0**: radial displacement, vertical force.
* **U_z0**: vertical displacement, vertical force.
* **U_r1**: radial displacement, horizontal force (× cosθ).
* **U_θ1**: tangential displacement, horizontal force (× sinθ).

| r0 | U_r0 (V-force) | U_z0 (V-force) | U_r1 (H-force) | U_θ1 (H-force) |
|---|---|---|---|---|
| 0.0 | −.027, .000 | .106, .000 | .159, .000 | −.106, .000 |
| 0.5 | −.032, .007 | .087, −.062 | .146, −.058 | −.089, .058 |
| 1.0 | −.033, .025 | .037, −.102 | .112, −.105 | −.045, .099 |
| 1.5 | −.020, .047 | −.029, −.108 | .063, −.133 | .015, .112 |
| 2.0 | .006, .060 | −.087, −.077 | .011, −.137 | .075, .093 |
| 2.5 | .041, .058 | −.120, −.017 | −.034, −.120 | .118, .045 |
| 3.0 | .074, .035 | −.114, .053 | −.064, −.090 | .132, −.020 |
| 3.5 | .092, −.005 | −.070, .110 | −.076, −.054 | .111, −.084 |
| 4.0 | .087, −.054 | −.001, .134 | −.075, −.024 | .059, −.132 |
| 4.5 | .056, −.096 | .072, .118 | −.065, −.004 | −.014, −.144 |
| 5.0 | .005, −.120 | .127, .064 | −.054, .004 | −.088, −.129 |
| 5.5 | −.054, −.116 | .144, −.013 | −.048, .003 | −.145, −.077 |

Apsel's own wavenumber-integration results with Q_α = Q_β = 5000 (0.01 % damping) match to ±0.005.
For example, at r0 = 3.0: U_r1 = −.068, −.090 and U_θ1 = .128, −.020.

The sign pattern (negative imaginary parts at small r0) is consistent with `e^{+iωt}`. Confirm this
when comparing; flip the sign of Im if the code uses `e^{−iωt}`.

**Use:** the thin-layer/POINT module, run with many thin layers over a paraxial or viscous
half-space boundary (or a very deep damped layer), should reproduce this table to ±0.005 absolute
for r0 ≤ 5.5. This is the single most direct check of SASSI-type layered Green functions.

### D.4 Lamb's problem, time domain (E, ν = 1/4)

Vertical step load P·H(t) at the surface. The vertical surface displacement at distance r is
positive in the direction of the load. Pekeris (1955a), *PNAS* 41:469–480. The structure is
confirmed by Richards (1979), *BSSA* 69(4):947–956, Eqs. 9–13
([pdf](https://www.ldeo.columbia.edu/~richards/my_papers/Richards_BSSA1979_LambsProblem.pdf)).

Define τ = Vs t/r and γ_R = √((3+√3)/4) = 1.087664, so that c_R/Vs = 0.919402.
```
w·πμr/P = 0                                                               τ < 1/√3
        = (1/32)[6 − √3/√(τ²−1/4) − √(3√3+5)/√(3/4+√3/4−τ²) + √(3√3−5)/√(τ²−3/4+√3/4)]   1/√3 < τ < 1
        = (1/16)[6 − √(3√3+5)/√(3/4+√3/4−τ²)]                              1 < τ < γ_R
        = 3/8   (static Boussinesq value (1−ν)/2)                          τ > γ_R
```
Check values (D): τ = 0.6: −0.018957; τ = 0.8: −0.010234; τ = 0.9: −0.028501; τ = 1⁻ and 1⁺:
−0.09151 (continuous at S); τ = 1.05: −0.32834; τ = 1.08: −1.17338; τ ≥ 1.0877: +0.375.

Self-checks:

* Continuity at P (to 1e-6), continuity at S, and the static limit are verified here.
* The constants match Richards' partial-fraction roots R1 = 1/4, R2 = (3−√3)/4 and R3 = (3+√3)/4.
* The remaining risk is a global sign of the dynamic part; confirm it against Pekeris' Fig. 1
  **[VERIFY]**.

Use: inverse-FFT check of the frequency-domain Green function, and a check of the Rayleigh arrival.

### D.5 Surface waves (E for the secular equations; D for the numbers)

**Rayleigh (half-space):** `(2 − x²)² = 4√(1−x²)√(1−κ²x²)`, with x = c_R/Vs and κ² = (1−2ν)/(2(1−ν)).

| ν | 0 | 0.2 | 0.25 | 0.3 | 1/3 | 0.35 | 0.4 | 0.45 | 0.49 | 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| c_R/Vs | 0.874032 | 0.910996 | 0.919402 | 0.927413 | 0.932526 | 0.935013 | 0.942195 | 0.948960 | 0.954074 | 0.955313 |

For ν = 1/4 the exact value is `√(2 − 2/√3)` = 0.919402. Viktorov's approximation
`(0.862+1.14ν)/(1+ν)` is good to about 1.5 %.

**Love waves (layer H, β1, μ1 over half-space β2, μ2):** `tan[ωH√(1/β1²−1/c²)] = μ2√(1/c²−1/β2²) /
[μ1√(1/β1²−1/c²)]`, with β1 < c < β2.

* Cut-off of mode n: `f_n = nβ1/[2H√(1−β1²/β2²)]`.
* Example: β1 = 200 and β2 = 400 m/s, equal ρ, H = 10 m. Cut-offs are 0, 11.547 and 23.094 Hz.
  Fundamental c = 387.536 (2 Hz), 300.026 (5 Hz), 224.175 (10 Hz), 205.949 m/s (20 Hz). Mode-1
  c = 279.222 m/s at 20 Hz.
* Use: the thin-layer eigenvalue problem (SH part) must converge to these values. Typical error is
  O((k h)²), approached from above, as the number of sublayers increases.

**Thin-layer method.** Kausel & Peek (1981, MIT R81-13;
[PB82147893](https://nehrpsearch.nist.gov/static/files/NSF/PB82147893.pdf)) compare the discrete
(thin-layer) and "exact" (wavenumber-integration) Green functions. They use a homogeneous stratum
with surface disk loads, discretised with N = 4, 6 and 12 sublayers (graphs only). The discretisation
always **stiffens** the system: computed compliance is below exact. Use this monotonic-convergence
property as a test: compliance should increase towards the exact value as N doubles.

---

## E. Published SSI benchmark problems relevant to SASSI-type codes

| ID | Problem | What is available | Source |
|---|---|---|---|
| E1 | Rigid massless disk on layer over rigid base, flexible volume method vs LUCON | Full definition (C.7). Graphical results. | Tajirian & Tabatabaie 1985, ASCE Convention preprint, Detroit (Bechtel); no stable URL found |
| E2 | **NUREG/CR-6896** (BNL 2006): deeply embedded structures (DEB), SASSI vs CARES, plus LS-DYNA for wall pressures | Full definition below. Results graphical. | [NUREG/CR-6896, ML060820521](https://www.nrc.gov/docs/ML0608/ML060820521.pdf) |
| E3 | DOE subtraction-method anomaly test cases (embedded box, uniform and layered sites); direct vs subtraction vs modified subtraction | Problem description, qualitative and graphical | Mertz, Costantino, Houston, Maham (2011), DOE NPH Workshop; DOE OE-3 2011-02 "SASSI Software Problem"; DOE SSI Report July 2011 |
| E4 | SASSI2010 SM/ESM/direct vs **ANSYS full harmonic** and SAP2000 time history (coarse 3-D quarter model; complex 3-D model) | Problem description, graphical | Anderson & Ostadan (2014), DOE NPH Workshop |
| E5 | SDOF damping estimation, fixed base and on half-space/layer (SASSI2000) | Tabulated damping (A.0 and C.7) | Ostadan, Deng & Roesset (2004), 3rd UJNR SSI Workshop |
| E6 | Kinematic interaction formulas (embedment, base-slab averaging) | Closed-form approximations (E7 below) | NIST GCR 12-917-21 ch. 3 |
| E7 | Seismic Safety Margins Research Program: CLASSI vs FLUSH comparisons (SMACS chain) | Reports | [NUREG/CR-2015 Vol. 1](https://www.nrc.gov/docs/ml1309/ml13093a113.pdf); Vol. 4 (SSI, Project III) |
| E8 | IAEA KARISMA (Kashiwazaki-Kariwa Unit 7, 2007 NCOE earthquake) — validation against records | Records, multi-team results | [IAEA-TECDOC-1722 (2013)](https://www-pub.iaea.org/MTCD/Publications/PDF/TE-1722_web.pdf) |
| E9 | Lotung and Hualien large-scale SSI experiments (EPRI) — validation | Records (EPRI-licensed) | EPRI NP-6154 (Lotung), EPRI TR-103407 Hualien [VERIFY numbers] |

**E2 definition (NUREG/CR-6896 §3.1.1, §5).**

Structure (sample problems):

* Cylinder 46 m tall, outer diameter 27 m, uniform 2 m wall: A = 157.1 m², I = 12,272 m⁴.
* Weights: wall plus equipment 92,202 kN; basemat 40,474 kN; roof 4,448 kN.
* Concrete: f'c = 27,579 kPa, ν = 0.2, damping 2 %, 23.2 kN/m³.
* Depth of burial E = 11.5, 23, 34.5 and 46 m (E/R = 0.85, 1.7, 2.55, 3.4).

Soils, all to bedrock at **80 m**:

* Column A: uniform, Vs = 250 m/s.
* Column B: uniform, Vs = 1000 m/s.
* Column C: `Vs = 250 + 250.78 z^0.25` (z in m; gives 1000 m/s at 80 m).
* Common properties: γ = 17.28 kN/m³, ν = 0.3, ξ = 4 %.
* Bedrock: Vs = 2000 m/s, 17.28 kN/m³, 2 %.
* No strain degradation. Input is a rock-outcrop motion matched to an RG 1.165 / NUREG/CR-6728
  spectrum, M7 at 25 km, PGA 0.3 g.

Benchmark set (§5): the same cylinder with internal walls (1 m shell, 3 m basemat, 1 m roof, reactor
and PCU vessels as lumped beams, total 506,159 kN, 5 % damping), and an ABWR-like 60 × 60 × 83 m box
embedded 33.75 m.

Closed-form checks in the report:

| DOB (m) | Column A | Column B | Column C |
|---|---|---|---|
| 11.5 | 5.4 | 21.7 | 13.5 |
| 23 | 2.7 | 10.9 | 7.5 |
| 34.5 | 1.8 | 7.3 | 5.3 |
| 46 | 1.4 | 5.4 | 4.2 |

These are soil-column frequencies V/(4·DOB) in Hz. Use them as the expected locations of embedment
"notches" in kinematic transfer functions.

**E3/E4 lessons, which define a test.**

* The subtraction method deviates from the direct (flexible-volume) method at and above the first
  natural frequency of the excavated soil volume. This happens when the structure is embedded, the
  foundation is wide and shallow, and the response frequencies are at or above that excavation
  frequency (DOE OE-3 2011-02). It is not a mesh (λ/5) effect.
* Generic study (Mertz et al. 2011): Western US site, 120 × 120 × 30 ft excavation, uniform 6 ft
  bricks (λ/5 = 29.2 Hz). The anomaly appears near 10 Hz, at the excavation frequency.
* Anderson & Ostadan (2014), coarse model:
  * Uniform soil Vs = 2500 fps, Vp = 4677 fps, 5 % damping.
  * Foundation at −140 ft, roof +50 ft, bedrock at −240 ft (fixed).
  * Soil box 600 × 600 ft with roller sides, cut-off about 15 Hz.
  * SASSI (rigid base) matches the SAP2000/ANSYS box model.

  This is an **ANSYS-reproducible** verification: a full harmonic analysis of a large soil box on a
  rigid base gives the same transfer functions as SASSI with a rigid-base layered site. The
  condition is that the box is large and damped enough for lateral reflections to be negligible.
  Recommended test: "flexible volume (all excavated nodes interacting) ≡ ANSYS harmonic box";
  "subtraction ≈ flexible volume below f_excavation".

**E6 kinematic interaction formulas (A)** (NIST GCR 12-917-21 Eqs. 3-1 to 3-6; Elsabee & Morray
1977; Kausel et al. 1978; Day 1978; Veletsos & Prasad 1989):

* Embedment, translation: `H_u = cos(Dω/Vs)` for Dω/Vs < 1.1, else 0.45.
* Embedment, rocking: `H_θ = θL/u_g = 0.26[1 − cos(Dω/Vs)]` for Dω/Vs < π/2, else 0.26.
* Base-slab averaging, wave passage (rectangle): `H_u = sin(a0^k Vs/V_app)/(a0^k Vs/V_app)` with
  `a0^k = ωB_e/Vs`, `B_e = √(BL)`; and H_u = 2/π beyond a0^k = πV_app/(2Vs).
* Incoherent field (Veletsos & Prasad): `H_u = {b0⁻²[1 − e^{−2b0²}(I0(2b0²) + I1(2b0²))]}^{1/2}`
  with `b0 = √(4/π) κ_a a0^k`. Map κ_a to the code's coherency parameter before use **[VERIFY]**.
* Exact limit: a **surface** rigid massless foundation under vertically incident coherent SH waves
  has foundation input motion equal to the free-field motion (H_u = 1, H_θ = 0) at every frequency.
  Any deviation in SASSI is discretisation error. This is a mandatory zero-tolerance-ish test,
  |H_u − 1| < 1 %.

---

## F. Structure on a flexible base (inertial SSI)

### F.1 Period lengthening and effective damping (E for massless foundation with real springs; A otherwise)

From Veletsos & Meek (1974) and NIST GCR 12-917-21 Eqs. 2-5 to 2-11:
```
(T̃/T)² = 1 + k/K_x + k h²/K_θ             (exact for massless rigid foundation, constant real springs)
β0 = β_f + β_i /(T̃/T)^n,   n = 3 (viscous structural damping), n = 2 (hysteretic)  (Givens 2013)
β_f ≈ [((T̃/T)^2 − 1)/(T̃/T)^2] β_s + β_x/(T̃/T_x)^2 + β_yy/(T̃/T_yy)^2 ,
       T_x = 2π√(m/K_x),  T_yy = 2π√(m h²/K_yy)          (Wolf 1985 form, Eq. 2-11)
```
**Worked example (D).**

* Inputs: circular foundation R = 10 m on a half-space, Vs = 100 m/s, ρ = 2 t/m³, ν = 1/3
  (G = 20,000 kPa). SDOF with h = 10 m, mass ratio m/(ρπR²h) = 0.15 (m = 942.48 t), fixed-base
  T = 0.5 s (k = 148,830 kN/m).
* Springs: K_x = 8GR/(2−ν) = 9.60×10⁵ kN/m; K_θ = 8GR³/(3(1−ν)) = 8.00×10⁷ kN·m.
* Veletsos–Meek: **T̃/T = 1.15805** (T̃ = 0.5790 s); structural damping term 0.05/(T̃/T)³ = 0.0322.
* Exact frequency-domain 3-DOF with Veletsos–Verbic impedances (massless foundation, β_struct = 5 %
  hysteretic): resonance at **1.717 Hz (T̃ = 0.5824 s)**, peak |ü_mass/ü_g| = 6.71.

The two T̃ values differ because k_r < 1 at a0 = 1.085. This is a good test of the coupling of
frequency-dependent impedances.

### F.2 Exact 3-DOF frequency-domain model (E, for an impedance-driven code path)

Unknowns relative to the free field: u0 (foundation translation), θ (rocking) and u (structural
deformation). Inputs: structure m, k* = k(1+2iβ) at height h; foundation mass m0 and inertia I0;
impedances S_x(ω) and S_θ(ω); input ü_g, with ü_total = ü_g + ü0 + hθ̈ + ü.
```
row 1 (structure)      : (k* − ω²m) u  −  ω²m u0            −  ω²mh θ                    = −m ü_g
row 2 (total shear)    :  −ω²m u      + (S_x − ω²(m+m0)) u0 −  ω²mh θ                    = −(m+m0) ü_g
row 3 (moment at base) :  −ω²mh u     −  ω²mh u0            + (S_θ − ω²(mh²+I0)) θ        = −mh ü_g
total acceleration of the mass: ü_t = ü_g − ω²(u0 + hθ + u)
```
With m0 = I0 = 0 and real constant S, the resonance reproduces F.1 exactly. This verifies the
assembly of (structure + interaction impedance + seismic load vector) in ANALYS. Feed it with
impedances computed by the code itself; tolerance 1e-8.

### F.3 Cone and lumped models

Use C.5 and C.6 to build closed-form SDOF-on-soil checks: a rigid block of mass m on a half-space.
For example, the vertical natural frequency `f ≈ (1/2π)√(K_z/(m + ΔM))` and damping
`D_z ≈ 0.425/√B_z`.

---

## G. Response spectra and spectrum-compatible motions (MOTION/EQUAKE)

### G.1 Exact piecewise-linear integration (E) — Nigam & Jennings (1969)

Nigam, N.C. & Jennings, P.C. (1969), *BSSA* 59(2):909–922. Solve `ẍ + 2ζωẋ + ω²x = −a(t)` with a(t)
linear inside each Δt: `{x,ẋ}_{i+1} = A{x,ẋ}_i + B{a_i, a_{i+1}}`.

Let `ω_d = ω√(1−ζ²)`, `E = e^{−ζωΔt}`, `S = sin ω_dΔt`, `C = cos ω_dΔt`, `r = ζ/√(1−ζ²)`,
`t1 = (2ζ²−1)/(ω²Δt)` and `t2 = 2ζ/(ω³Δt)`.
```
a11 = E(rS + C)                 a12 = E S/ω_d
a21 = −(ω/√(1−ζ²)) E S          a22 = E(C − rS)
b11 = E[(t1 + ζ/ω) S/ω_d + (t2 + 1/ω²) C] − t2
b12 = −E[t1 S/ω_d + t2 C] − 1/ω² + t2
b21 = E[(t1 + ζ/ω)(C − rS) − (t2 + 1/ω²)(ω_d S + ζωC)] + 1/(ω²Δt)
b22 = −E[t1 (C − rS) − t2 (ω_d S + ζωC)] − 1/(ω²Δt)
absolute acceleration = −(2ζωẋ + ω²x)
```
These coefficients were verified here: a step input gives `max|x − x_exact|·ω² < 3e-12`.

**Exact spectral values for test inputs (E):**

| Input | Quantity | Exact value | NJ result here |
|---|---|---|---|
| Step a(t) = a0 (t ≥ 0), from rest | SD | `(a0/ω²)(1 + e^{−ζπ/√(1−ζ²)})`; ζ = 5 %: PSA/a0 = **1.854468** at every frequency | 1.854466 (1 Hz), 1.854461 (5, 10 Hz; Δt = 1 ms) |
| Velocity impulse ΔV (a = ΔVδ(t)) | SD | `(ΔV/ω) exp[−ζ arccos ζ/√(1−ζ²)]`; ζ = 5 %: PSV/ΔV = 0.92669 | — |
| Harmonic a = A sin ωt at resonance, long duration | PSA → steady state | `A/(2ζ)` (ζ = 5 %: 10.000; 2 %: 25.000) | 9.9992; 24.998 (200 s) |
| Same | max \|abs. acc.\| → | `A√(1+4ζ²)/(2ζ)` (5 %: 10.0499) | 10.0489 |
| Harmonic, ratio r = Ω/ω, steady state | SD | `(A/ω²)/√((1−r²)² + (2ζr)²)` | — |
| Any record | f → ∞ | SA → PGA (ZPA) | — |
| Any record | f → 0 | SD → PGD (peak ground displacement) | — |

Tolerances: 1e-5 relative for the step at Δt ≤ T/100. Report both PSA = ω²SD and true absolute SA
(they differ by O(ζ) at low frequency).

### G.2 SRP 3.7.1 (Rev. 4, Dec 2014) acceptance criteria (P — regulatory text)

[NUREG-0800 SRP 3.7.1 Rev. 4, ML14198A460](https://www.nrc.gov/docs/ML1419/ML14198A460.pdf)

**General requirements:**

* Three orthogonal components, statistically independent: |correlation coefficient| ≤ **0.16**.
  Shifting the start time does not make a new history.
* Strong-motion duration (Arias intensity 5 % to 75 %) ≥ **6 s**.
* Acceleration, velocity and displacement must be compatible, with no baseline drift.
* V/A and AD/V² should be consistent with the controlling events.

**Option 1, Approach 1 (enveloping):**

* The spectrum envelops the target for all damping values. "Envelops" means no more than **5
  points** below the target and no point more than **10 %** below.
* Table 3.7.1-1 frequency increments (Hz):

  | Range | Increment |
  |---|---|
  | 0.2–3.0 | 0.10 |
  | 3.0–3.6 | 0.15 |
  | 3.6–5.0 | 0.20 |
  | 5.0–8.0 | 0.25 |
  | 8.0–15.0 | 0.50 |
  | 15–18 | 1.0 |
  | 18–22 | 2.0 |
  | 22 up | 3.0 |

* A single time history must also satisfy a PSD check.

**Option 1, Approach 2 (mean-based fit):**

* (a) Nyquist ≥ 50 Hz (Δt ≤ 0.010 s); total duration ≥ 20 s (zero-padding allowed).
* (b) 5 % spectrum at ≥ **100 points per frequency decade**, log-uniform, from 0.1 Hz to 50 Hz (or
  Nyquist).
* (c) Never more than **10 % below** the target at any frequency. Below-target points are allowed
  only within a window ≤ ±10 % of frequency, i.e. **≤ 9 adjacent points** below.
* (d) Never more than **30 % above** (factor 1.3). The PSD must show no significant gaps.

**Option 2 (multiple sets):** ≥ 4 histories for linear analysis. (a) and (b) apply to each history;
(c) and (d) apply to the suite mean.

**Appendix A — minimum PSD for RG 1.60 (1 g; scale by PGA²).** The one-sided PSD is
`S0(ω) = 2|F(ω)|²/(2πT_D)`. The ±20 %-band average must exceed **80 %** of the target between
0.3 and 24 Hz:
```
f < 2.5 Hz    : S0 = 0.419 m²/s³ · (f/2.5)^0.2
2.5–9.0 Hz    : S0 = 0.419 m²/s³ · (2.5/f)^1.8
9.0–16.0 Hz   : S0 = 418 cm²/s³ · (9.0/f)^3
f > 16 Hz     : S0 = 74.2 cm²/s³ · (16.0/f)^8
```

### G.3 RG 1.60 design spectra (P) — a standard target for EQUAKE tests

[RG 1.60 Rev. 2, ML13210A432](https://www.nrc.gov/docs/ML1321/ML13210A432.pdf). Spectra are anchored
to 1.0 g with PGD = 36 in. Amplification factors are linear between control points on log–log
(tripartite) axes.

| Damping % | Horiz. A (33 Hz) | Horiz. B (9 Hz) | Horiz. C (2.5 Hz) | Horiz. D (0.25 Hz, displ.) | Vert. A (33 Hz) | Vert. B (9 Hz) | Vert. C (3.5 Hz) | Vert. D (0.25 Hz, displ.) |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 1.0 | 4.96 | 5.95 | 3.20 | 1.0 | 4.96 | 5.67 | 2.13 |
| 2 | 1.0 | 3.54 | 4.25 | 2.50 | 1.0 | 3.54 | 4.05 | 1.67 |
| 5 | 1.0 | 2.61 | 3.13 | 2.05 | 1.0 | 2.61 | 2.98 | 1.37 |
| 7 | 1.0 | 2.27 | 2.72 | 1.88 | 1.0 | 2.27 | 2.59 | 1.25 |
| 10 | 1.0 | 1.90 | 2.28 | 1.70 | 1.0 | 1.90 | 2.17 | 1.13 |

Derived (D), 5 %, at 1 g:

* Horizontal SA at D (0.25 Hz) = ω²·(2.05 × 36 in) = **0.4716 g**. Below 0.25 Hz the spectrum is
  constant-displacement at 73.8 in.
* Vertical SA at D = **0.3152 g** (49.32 in).

Test: (i) a target generator reproduces these points exactly; (ii) a generated motion passes the G.2
checks.

---

## H. Site-response reference data (SHAKE91 example; modulus and damping curves)

### H.1 SHAKE91 sample problem (P) — the primary SOIL/SITE benchmark

Idriss & Sun (1992), *User's Manual for SHAKE91*, UC Davis, Tables B-1 and B-2. Sources:

* Manual and files are public domain via [github.com/ocrickard/SHAKE16](https://github.com/ocrickard/SHAKE16).
  Its `test-data` output is identical to the original executable.
* EERA reproduces the problem to within 0.2 % in strain and 0.4 % in Sa
  ([EERA manual App. A–B](http://www.ce.memphis.edu/7137/PDFs/EERA2/EERAManual.pdf)).
* Local copies: `docs/spec/R2_benchmark_data/shake91_example/` (INP.DAT, DIAM.ACC, output1).

Input:

* Profile: 150 ft over a half-space, 17 sublayers, English units. Unit weight is in kcf; mass uses
  g = 32.2 ft/s². For example G1 = 0.125/32.2 × 1000² = 3882 ksf.
* Initial damping 5 % in soil, 1 % in rock.
* Motion: DIAM.ACC, 1989 Loma Prieta, Diamond Heights H1_90. dt = 0.02 s; 1900 points read;
  NFFT = 4096. Scaled from 0.11289 g to **0.10 g**; frequencies above **25 Hz** removed.
* Applied as an **outcrop** motion at sublayer 17 (the half-space).
* Iterations: **8**. Effective strain ratio **0.50**.
* Spectra at 5 % damping for layer 1 (outcrop surface), with g = 981 cm/s².

| No. | Type | Thick. (ft) | Depth mid (ft) | Unit wt (kcf) | Vs (ft/s) | Gmax (ksf) |
|---|---|---|---|---|---|---|
| 1 | 2 (sand) | 5 | 2.5 | 0.125 | 1000 | 3882 |
| 2 | 2 | 5 | 7.5 | 0.125 | 900 | 3144 |
| 3 | 2 | 10 | 15 | 0.125 | 900 | 3144 |
| 4 | 2 | 10 | 25 | 0.125 | 950 | 3503 |
| 5 | 1 (clay) | 10 | 35 | 0.125 | 1000 | 3882 |
| 6 | 1 | 10 | 45 | 0.125 | 1000 | 3882 |
| 7 | 1 | 10 | 55 | 0.125 | 1100 | 4697 |
| 8 | 1 | 10 | 65 | 0.125 | 1100 | 4697 |
| 9 | 2 | 10 | 75 | 0.130 | 1300 | 6823 |
| 10 | 2 | 10 | 85 | 0.130 | 1300 | 6823 |
| 11 | 2 | 10 | 95 | 0.130 | 1400 | 7913 |
| 12 | 2 | 10 | 105 | 0.130 | 1400 | 7913 |
| 13 | 2 | 10 | 115 | 0.130 | 1500 | 9084 |
| 14 | 2 | 10 | 125 | 0.130 | 1500 | 9084 |
| 15 | 2 | 10 | 135 | 0.130 | 1600 | 10335 |
| 16 | 2 | 10 | 145 | 0.130 | 1800 | 13081 |
| 17 | 3 (rock, half-space) | — | — | 0.140 | 4000 | 69565 |

Curves (strain in %):

| γ (%) | Clay G/Gmax (Seed & Sun 1989, upper) | Clay D % (Idriss 1990) | Sand G/Gmax (Seed & Idriss 1970, upper) | Sand D % (Idriss 1990) |
|---|---|---|---|---|
| 0.0001 | 1.000 | 0.24 | 1.000 | 0.24 |
| 0.0003 | 1.000 | 0.42 | 1.000 | 0.42 |
| 0.001 | 1.000 | 0.80 | 0.990 | 0.80 |
| 0.003 | 0.981 | 1.40 | 0.960 | 1.40 |
| 0.01 | 0.941 | 2.80 | 0.850 | 2.80 |
| 0.03 | 0.847 | 5.10 | 0.640 | 5.10 |
| 0.1 | 0.656 | 9.80 | 0.370 | 9.80 |
| 0.3 | 0.438 | 15.50 | 0.180 | 15.50 |
| 1 | 0.238 | 21.0 | 0.080 | 21.0 |
| 3 (damping at 3.16) | 0.144 | 25.0 | 0.050 | 25.0 |
| 10 | 0.110 | 28.0 | 0.035 | 28.0 |

Rock (Schnabel 1973, "average"):

* G/Gmax at γ = 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 1.0 % is 1.000, 1.000, 0.9875,
  0.9525, 0.900, 0.810, 0.725, 0.550.
* Damping at γ = 0.0001, 0.001, 0.01, 0.1, 1 % is 0.4, 0.8, 1.5, 3.0, 4.6 %.
* Interpolation is SHAKE91's linear interpolation in log(strain).

**Published results (SHAKE91 manual Table B-2, final iteration 8):**

| Sublayer | Depth (ft) | Effective strain (%) | Damping | G (ksf) | G/Gmax | Max strain (%) | Max stress (psf) |
|---|---|---|---|---|---|---|---|
| 1 | 2.5 | 0.00077 | 0.007 | 3851.5 | 0.992 | 0.00154 | 59.43 |
| 2 | 7.5 | 0.00295 | 0.014 | 3020.0 | 0.960 | 0.00591 | 178.41 |
| 3 | 15 | 0.00634 | 0.023 | 2803.8 | 0.892 | 0.01267 | 355.31 |
| 4 | 25 | 0.00976 | 0.028 | 2985.8 | 0.852 | 0.01952 | 582.73 |
| 5 | 35 | 0.01099 | 0.030 | 3621.7 | 0.933 | 0.02197 | 795.85 |
| 6 | 45 | 0.01403 | 0.035 | 3540.5 | 0.912 | 0.02806 | 993.46 |
| 7 | 55 | 0.01362 | 0.034 | 4296.0 | 0.915 | 0.02723 | 1169.83 |
| 8 | 65 | 0.01566 | 0.037 | 4239.8 | 0.903 | 0.03132 | 1327.93 |
| 9 | 75 | 0.01356 | 0.034 | 5402.8 | 0.792 | 0.02711 | 1464.73 |
| 10 | 85 | 0.01505 | 0.037 | 5266.1 | 0.772 | 0.03011 | 1585.43 |
| 11 | 95 | 0.01336 | 0.034 | 6288.3 | 0.795 | 0.02671 | 1679.79 |
| 12 | 105 | 0.01413 | 0.035 | 6203.6 | 0.784 | 0.02825 | 1752.67 |
| 13 | 115 | 0.01233 | 0.032 | 7357.0 | 0.810 | 0.02467 | 1814.84 |
| 14 | 125 | 0.01282 | 0.033 | 7290.6 | 0.803 | 0.02563 | 1868.58 |
| 15 | 135 | 0.01115 | 0.030 | 8570.2 | 0.829 | 0.02230 | 1911.02 |
| 16 | 145 | 0.00865 | 0.026 | 11292.4 | 0.863 | 0.01729 | 1952.92 |

All values are as printed in the manual. The SHAKE16 rerun reproduces them to within ±0.1 in the
last printed digit; for example sublayer 9 G = 5402.7 vs 5402.8, and surface PGA 0.19040 vs 0.19037.

**Peak acceleration profile (g), Table B-2 (Option 6):**

| Location | PGA |
|---|---|
| Surface (outcrop) | **0.19037** |
| 5 ft | .19006 |
| 10 ft | .18876 |
| 20 ft | .18258 |
| 30 ft | .17208 |
| 40 ft | .15947 |
| 50 ft | .14288 |
| 60 ft | .12652 |
| 70 ft | .11050 |
| 80 ft | .09840 |
| 90 ft | .08999 |
| 100 ft | .08268 |
| 110 ft | .08559 |
| 120 ft | .08547 |
| 130 ft | .08198 |
| 140 ft | .07769 |
| 150 ft within | **.07617** |
| 150 ft outcrop | .10000 |

Times of peak: 11.28 s at the surface, 10.92 s at the base.

**Other published values:**

* Initial (elastic, 5 %) profile: T = 0.48 s from V_avg = 1253 ft/s; maximum amplification 13.80 at
  2.32 Hz.
* Final profile: T = 0.52 s from V_avg = 1153 ft/s; maximum amplification **20.47 at 2.11 Hz**.
  This is the SHAKE "layer 1 / base" amplification, interpreted as surface/within-base.
* Option 10, surface/rock-outcrop amplification: maximum **3.44 at 2.12 Hz (T = 0.47 s)**.

Selected values from SHAKE16 output (identical program):

* |surface/outcrop| = 1.3270 (1 Hz), 3.3181 (2 Hz), 3.4369 (2.125 Hz), 1.7386 (3 Hz), 2.0450 (5 Hz),
  1.6690 (10 Hz).
* 5 % surface spectrum, absolute acceleration (g):

| T (s) | 0.01 | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.47 | 0.50 | 0.60 | 0.80 | 1.0 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA (g) | 0.1898 | 0.2265 | 0.3977 | 0.3461 | 0.4811 | 0.7832 | 0.6366 | 0.4949 | 0.3812 | 0.1431 | 0.1351 | 0.0922 | 0.0631 | 0.0313 | 0.0139 |
| SD (cm) | 0.00047 | 0.01481 | 0.09603 | 0.34684 | 1.06696 | 3.11702 | 3.46976 | 3.04459 | 3.39854 | 2.24551 | 3.33190 | 5.09715 | 6.24650 | 6.94620 | 8.51406 |

Summary values:

* Over 0.1–2.5 s: maximum acceleration response 0.790 g; maximum velocity response 50.52 cm/s.
* Area under the acceleration spectrum 0.459; area under the velocity spectrum 67.06.

**Tolerances.** An equivalent-linear SOIL module using SHAKE91 conventions (SASSI complex modulus,
the same strain ratio and the same curves interpolated in log-strain) should match:

* PGA profile within 0.5 %.
* Effective strains, G and damping within 1 %.
* Spectra within 2 %.

Differences larger than this usually come from one of:

* the FFT length or zero-padding;
* the frequency cut-off at 25 Hz;
* the mass-from-unit-weight constant (32.2);
* interpolation of the curves (SHAKE91 interpolates linearly in log strain);
* the "outcrop" definition (2A_17).

### H.2 Darendeli (2001) curves (A — closed form; values D)

Darendeli (2001), UT Austin PhD; equations as in the Strata Technical Manual (Kottke & Rathje 2008,
PEER 2008/10) §2.1.3 ([pdf](https://peer.berkeley.edu/sites/default/files/web_peer810_albert_r._kottke_ellen_m._rathje.pdf)).
γ is in percent and σ'0 is in atm.
```
G/Gmax = 1/[1 + (γ/γ_r)^a],  a = 0.9190
γ_r(%) = (0.0352 + 0.0010·PI·OCR^0.3246) σ'0^0.3483
D_min(%) = (0.8005 + 0.0129·PI·OCR^−0.1069) σ'0^−0.2889 (1 + 0.2919 ln f)
D_Masing,a=1(%) = (100/π){4[γ − γ_r ln((γ+γ_r)/γ_r)]/[γ²/(γ+γ_r)] − 2}
D_Masing = c1 D_M1 + c2 D_M1² + c3 D_M1³ ;
   c1 = −1.1143a² + 1.8618a + 0.2533, c2 = 0.0805a² − 0.0710a − 0.0095, c3 = −0.0005a² + 0.0002a + 0.0003
D(%) = b (G/Gmax)^0.1 D_Masing + D_min,   b = 0.6329 − 0.0057 ln N   (N = 10, f = 1 Hz typical)
```
The Strata manual says "γ_r (not in percent)". With φ1 = 0.0352 the reference strain is in percent
(0.0352 % for clean sand at 1 atm), so treat that phrase as a typo **[VERIFY against Darendeli 2001
Table 8.x]**.

Exact property: G/Gmax = 0.5 at γ = γ_r.

| γ (%) | Sand, PI=0, OCR=1, 1 atm: G/Gmax | Sand: D % | Clay, PI=30, OCR=1, 1 atm: G/Gmax | Clay: D % |
|---|---|---|---|---|
| 0.0001 | 0.9955 | 0.839 | 0.9974 | 1.208 |
| 0.001 | 0.9635 | 1.175 | 0.9789 | 1.391 |
| 0.01 | 0.7607 | 3.959 | 0.8485 | 3.038 |
| 0.0352 (= γ_r sand) | 0.5000 | 8.655 | — | — |
| 0.1 | 0.2770 | 13.806 | 0.4030 | 11.131 |
| 0.3 | 0.1225 | 18.289 | — | — |
| 1.0 | 0.0441 | 20.736 | 0.0752 | 20.207 |

D_min: 0.8005 % for sand and 1.1875 % for PI = 30. γ_r for PI = 30 is 0.0652 %.

Other curve families should be digitised from the original papers and, ideally, cross-checked
against DEEPSOIL/Strata built-in tables before being used as "P" data. These include Vucetic &
Dobry (1991) PI-dependent curves, EPRI (1993) depth-dependent curves and Seed et al. (1986) gravel
curves.

---

## I. HOUSE element library — eigenvalue and static benchmarks

### I.1 Closed forms (E)

| Model | Exact result |
|---|---|
| Bar, fixed–free, axial | `f_n = (2n−1)/(4L) √(E/ρ)` |
| Shear column / soil column, fixed–free | `f_n = (2n−1) Vs/(4H)` (with discrete errors per A.3) |
| Euler–Bernoulli cantilever | `f_n = (β_nL)²/(2πL²) √(EI/m̄)`; β_nL = **1.875104, 4.694091, 7.854757, 10.995541, 14.137168** |
| Free–free or clamped–clamped beam | β_nL = 4.730041, 7.853205, 10.995608 |
| Simply supported beam | `f_n = (nπ)²/(2πL²) √(EI/m̄)` |
| Cantilever static | `δ_tip = PL³/(3EI)` (plus shear `PL/(κGA)` for Timoshenko) |
| Kirchhoff plate, simply supported a × b | `f_mn = (π/2)[(m/a)² + (n/b)²] √(D/(ρt))`, `D = Et³/[12(1−ν²)]` |
| Rigid block on springs | `f = (1/2π)√(K/m)` per decoupled DOF |

### I.2 NAFEMS free-vibration benchmarks (P)

NAFEMS, *The Standard NAFEMS Benchmarks*, TNSB Rev. 3 (1990), and NAFEMS R0015 (1987). Common
material: **E = 200 GPa, ν = 0.3, ρ = 8000 kg/m³**. Targets are as reproduced in the Abaqus
Benchmarks Guide §4.4 ([FV2](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf16.html),
[FV4](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf17.html),
[FV12](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf18.html),
[FV16](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf20.html),
[FV32](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf22.html),
[FV41](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf23.html),
[FV52](https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf25.html)) and Altair OptiStruct
OS-V (Test 21).

| Test | HOUSE element tested | Definition | Target frequencies (Hz) |
|---|---|---|---|
| FV2 pin-ended double cross (in-plane) | BEAM | 8 arms pinned at their ends. Arm geometry (5 m arms, 0.125 m square section, 4 elements per arm) **[VERIFY from TNSB figure]** | 11.336; 17.709 (×7, repeated); 45.345; 57.390 (×7) |
| FV4 cantilever with off-centre point masses | BEAM + MASS (offsets) | Geometry per TNSB figure **[VERIFY]** | 1.723, 1.727, 7.413, 9.972, 18.155, 26.957 |
| FV12 free thin square plate | SHELL (Kirchhoff) | 10 × 10 m, t = 0.05 m, unsupported (6 rigid modes; in-plane restrained in NAFEMS → 3 RBM) | RBM ×3; 1.622, 2.360, 2.922, 4.233, 4.233, 7.416 |
| FV16 cantilevered thin square plate | SHELL | 10 × 10 m, t = 0.05 m, clamped along one edge (y-axis) | 0.421, 1.029, 2.582, 3.306, 3.753, 6.555 |
| FV32 cantilevered tapered membrane | PLANE (plane stress) | Length 10 m, depth 5 m at the root tapering to 1 m at the tip **[VERIFY]**, t = 0.05 m, root clamped, out-of-plane restrained | 44.623, 130.03, 162.70, 246.05, 379.90, 391.44 |
| NAFEMS Test 21 (R0015) simply supported thick plate (Mindlin, closed form) | SHELL (thick) / SOLID | 10 × 10 × 1 m, simply supported, in-plane restrained | 45.897; 109.44 (×2); 167.89; 204.51 (×2); 256.50 (×2) |
| FV52 simply supported "solid" plate | SOLID | 10 × 10 × 1 m, uz = 0 along the 4 bottom edges; 3 in-plane RBM | Abaqus Guide: 44.092, 106.66 (×2), 156.23, 193.58, 200.13 (×2). Altair OS-V 0455: 45.897, 109.44 (×2), 167.89(?), 193.59, 206.19 (×2). **Conflicting published targets — [VERIFY with TNSB]** |
| FV41 free cylinder, axisymmetric | (axisymmetric; reference only) | wall 0.4 m, free | 243.53, 377.41, 394.11, 397.72, 405.28 |

The Abaqus guide shows typical accuracy of good elements: about 0.1–1 % (FV32 CPS8 within 0.6 %;
FV16 S8R5 within 1 % for modes 1–5). Mode 6 and the higher plate modes show 2–10 % errors with
coarse meshes. Set test tolerances per mesh: for example ≤ 2 % on the first 4 modes with the NAFEMS
mesh.

**ANSYS cross-check (recommended).** Build each benchmark once in ANSYS using the elements the
HOUSE converter supports: BEAM4/188, SHELL63/181, SOLID45/185, COMBIN14 and MASS21. Compare the
eigenvalues of the converted SASSI model with ANSYS to 0.1 %. This isolates converter and element
errors from mesh errors. For example, FV12 checks that SHELL63 and SASSI SHELL share the Kirchhoff
formulation.

---

## J. Recommended verification matrix (minimum set)

| # | Module | Benchmark | Type | Pass criterion |
|---|---|---|---|---|
| V1 | SITE | A.2 uniform layer TFs, both damping forms | E | 1e-8 relative |
| V2 | SITE | A.3 discrete column eigenvalues (consistent, lumped, 50/50) | D | 1e-6 |
| V3 | SOIL | H.1 SHAKE91 sample problem | P | PGA 0.5 %, strain 1 %, SA 2 % |
| V4 | POINT | D.3 Wong (1975) surface Green functions | P | ±0.005 abs |
| V5 | POINT | D.1 static limits and reciprocity | E | 1e-4 relative (discrete) |
| V6 | POINT | D.5 Rayleigh/Love speeds from thin-layer eigenproblem | E/D | O(h²) convergence, < 0.5 % at 10 sublayers/λ |
| V7 | ANALYS | B.1 disk static stiffness (fine mesh, extrapolated) | E | 3–5 % (welded vs relaxed) |
| V8 | ANALYS | C.2 Veletsos–Verbic k, c for a0 ≤ 1.5; C.1 asymptotes | A/E | ±10 % (k), ±10 % (c); asymptotes ±5 % |
| V9 | ANALYS | C.7 Tajirian–Tabatabaie layer case (resonance locations exact) | P/E | resonance A0 ±2 %; curves by digitised comparison |
| V10 | MOTION | E6 surface rigid foundation, vertical SH: H_u = 1 | E | < 1 % |
| V11 | MOTION | E6 embedded foundation vs Elsabee–Morray | A | qualitative ±20 % |
| V12 | HOUSE+ANALYS | F.2 3-DOF with code-computed impedances | E | 1e-8 |
| V13 | HOUSE+ANALYS | A.0 fixed-base hysteretic SDOF peak 1/(2β√(1−β²)) | E | 1e-6 |
| V14 | EQUAKE | G.1 NJ step/harmonic | E | 1e-5 |
| V15 | EQUAKE | G.2 SRP 3.7.1 Approach 2 checks on a generated RG 1.60 motion | P | all criteria pass |
| V16 | HOUSE | I.1 closed forms; I.2 NAFEMS FV12, FV16, FV32, Test 21 | E/P | ≤ 2 % first 4 modes (NAFEMS mesh) |
| V17 | ANALYS | E3/E4 flexible-volume vs subtraction vs ANSYS harmonic box | D | FV ≡ ANSYS box within 5 %; flag SM deviation above f_excavation |

---

## K. Open questions / items to confirm

1. **Original rigorous impedance tables** (Luco & Westmann 1971; Veletsos & Wei 1971; Wong & Luco
   1978/1985) were not obtainable from open sources. Obtain them, and also the Veletsos–Verbic
   ν = 0.45 vertical coefficients and the welded-vs-relaxed corrections, to turn V8 from an
   approximate-formula test into a published-number test.
2. **Kausel stratum vertical and torsional formulas** (1 + 1.28R/H; embedded vertical factor) are
   reproduced only from secondary memory. Confirm them in Kausel (1974, MIT PhD) or Kausel &
   Ushijima (1979, MIT R79-?).
3. **FV52 targets conflict** between the Abaqus guide and Altair. Confirm with NAFEMS TNSB Rev. 3.
   The FV2, FV4 and FV32 geometries also need the TNSB figures.
4. **Pekeris sign** of the pre-Rayleigh part (D.4): the structure is confirmed, the global sign is not.
5. **Tajirian & Tabatabaie (1985) curves:** digitise them (and the compliance normalisation) to give
   numeric targets for V9.
6. **Veletsos–Prasad κ_a vs ACS SASSI coherency parameter mapping** (E6) is needed before using it
   as an incoherency benchmark.
7. **Darendeli γ_r units** (percent vs fraction) as written in the Strata manual (H.2).
8. Confirm **which complex-modulus form ACS SASSI applies to springs and structural elements.** The
   SASSI2000 behaviour (A.0) supports the SASSI form for elements.

---

## L. References (with open-access links where found)

* Anderson, L. & Ostadan, F. (2014). *Validation of the SASSI2010 Subtraction Method Using Full
  Scale Independent Verification*. DOE NPH Workshop, Oct 2014 (slides).
* Apsel, R.J. (1979). *Dynamic Green's Functions for Layered Media and Applications to
  Boundary-Value Problems*. PhD diss., UCSD. https://nehrpsearch.nist.gov/static/files/NSF/PB81178204.pdf
* Bardet, J.P., Ichii, K. & Lin, C.H. (2000). *EERA — A Computer Program for Equivalent-linear
  Earthquake site Response Analyses*. USC. http://www.ce.memphis.edu/7137/PDFs/EERA2/EERAManual.pdf
* Costantino, C.J. (2009). *Consistent Site Response – SSI Calculations*. BNL Report
  N6112-051208 Rev. 1, for USNRC. https://www.nrc.gov/docs/ML0919/ML091980384.pdf
* DOE (2011). OE-3: 2011-02 *SASSI Software Problem* (Operating Experience Level 3), and *DOE
  Soil-Structure Interaction Report*, July 2011.
* Elsabee, F. & Morray, J.P. (1977). *Dynamic Behavior of Embedded Foundations*. MIT R77-33.
  https://nehrpsearch.nist.gov/static/files/NSF/PB286493.pdf
* Gazetas, G. (1991). Formulas and charts for impedances of surface and embedded foundations.
  *J. Geotech. Eng.* 117(9):1363–1381.
  http://ssi.civil.ntua.gr/downloads/journals/1991-ASCE_FORMULAS%20AND%20CHARTS%20FOR%20IMPEDANCES%20OF%20SURFACE%20AND%20EMBEDDED%20FOUNDATIONS.pdf
* IAEA (2013). *Review of Seismic Evaluation Methodologies for NPPs Based on a Benchmark Exercise*
  (KARISMA). IAEA-TECDOC-1722. https://www-pub.iaea.org/MTCD/Publications/PDF/TE-1722_web.pdf
* Idriss, I.M. & Sun, J.I. (1992). *User's Manual for SHAKE91*. UC Davis. Public-domain copy and
  test data: https://github.com/ocrickard/SHAKE16
* Jakub, M. & Roesset, J.M. (1977). *Dynamic Stiffness of Foundations: 2-D vs 3-D Solutions*. MIT
  R77-36. https://nehrpsearch.nist.gov/static/files/NSF/PB286504.pdf
* Kausel, E. (1981). *An Explicit Solution for the Green Functions for Dynamic Loads in Layered
  Media*. MIT R81-13. https://nehrpsearch.nist.gov/static/files/NSF/PB82147893.pdf
* Kausel, E. (2006). *Fundamental Solutions in Elastodynamics: A Compendium*. Cambridge Univ. Press.
* Kausel, E., Whitman, R.V., Morray, J.P. & Elsabee, F. (1978). The spring method for embedded
  foundations. *Nucl. Eng. Des.* 48:377–392.
* Kottke, A.R. & Rathje, E.M. (2008). *Technical Manual for Strata*. PEER 2008/10.
  https://peer.berkeley.edu/sites/default/files/web_peer810_albert_r._kottke_ellen_m._rathje.pdf
* Kramer, S.L. (1996). *Geotechnical Earthquake Engineering*. Prentice Hall (ch. 7).
* Luco, J.E., Frazier, G.A., Day, S.M. & Apsel, R.J. (1978). *Dynamic Response of
  Three-Dimensional Rigid Embedded Foundations*. Final report, NSF ENV 76-22632, UCSD.
  https://nehrpsearch.nist.gov/static/files/NSF/PB296617.pdf
* Luco, J.E. & Westmann, R.A. (1971). Dynamic response of circular footings. *J. Eng. Mech. Div.
  ASCE* 97(EM5):1381–1395.
* Lysmer, J., Tabatabaie, M., Tajirian, F., Vahdani, S. & Ostadan, F. (1981). *SASSI — A System
  for Analysis of Soil-Structure Interaction*. UCB/GT/81-02.
* Mertz, G.E., Costantino, M.C., Houston, T.W. & Maham, A.S. (2011). *SASSI Subtraction Method
  Effects at Various DOE Projects*. DOE NPH Workshop, Oct 25–26 2011.
* NIST (2012). *Soil-Structure Interaction for Building Structures*. NIST GCR 12-917-21.
  https://www.nehrp.gov/pdf/nistgcr12-917-21.pdf
* Nigam, N.C. & Jennings, P.C. (1969). Calculation of response spectra from strong-motion
  earthquake records. *BSSA* 59(2):909–922.
* Ostadan, F., Deng, N. & Roesset, J.M. (2004). Estimating total system damping for
  soil-structure interaction systems. *Proc. 3rd UJNR Workshop on SSI*, Menlo Park.
  https://www.pwri.go.jp/eng/ujnr/tc/a/ssi_w3/Contributions/Ostadan.pdf
* Pais, A. & Kausel, E. (1988). Approximate formulas for dynamic stiffnesses of rigid foundations.
  *Soil Dyn. Earthq. Eng.* 7(4):213–227.
* Pekeris, C.L. (1955). The seismic surface pulse. *PNAS* 41:469–480.
* Richards, P.G. (1979). Elementary solutions to Lamb's problem for a point source…. *BSSA*
  69(4):947–956. https://www.ldeo.columbia.edu/~richards/my_papers/Richards_BSSA1979_LambsProblem.pdf
* Richart, F.E., Hall, J.R. & Woods, R.D. (1970). *Vibrations of Soils and Foundations*.
  Prentice Hall.
* Shye, K.-Y. & Robinson, A.R. (1980). *Dynamic Soil-Structure Interaction*. UIUC SRS 484
  (contains the Veletsos–Verbic table). https://nehrpsearch.nist.gov/static/files/NSF/PB81144503.pdf
* Tajirian, F.F. & Tabatabaie, M. (1985). Vibration analysis of foundations on layered media.
  *ASCE Convention, Vibration Problems in Geotechnical Engineering*, Detroit.
* USNRC (2014). NUREG-0800 SRP 3.7.1 Rev. 4 *Seismic Design Parameters*.
  https://www.nrc.gov/docs/ML1419/ML14198A460.pdf
* USNRC (2014). RG 1.60 Rev. 2 *Design Response Spectra for Seismic Design of NPPs*.
  https://www.nrc.gov/docs/ML1321/ML13210A432.pdf
* Veletsos, A.S. & Meek, J.W. (1974). Dynamic behaviour of building-foundation systems.
  *EESD* 3:121–138.
* Veletsos, A.S. & Verbic, B. (1973). Vibration of viscoelastic foundations. *EESD* 2:87–102.
* Veletsos, A.S. & Wei, Y.T. (1971). Lateral and rocking vibration of footings. *J. SMFD ASCE*
  97(SM9):1227–1248.
* Wolf, J.P. (1994). *Foundation Vibration Analysis Using Simple Physical Models*. Prentice Hall.
  Wolf & Deeks (2004), Elsevier.
* Wong, H.L. (1975). *Dynamic Soil-Structure Interaction*. Caltech EERL 75-01.
* Xu, J., Miller, C., Costantino, C. & Hofmayer, C. (2006). *Assessment of Seismic Analysis
  Methodologies for Deeply Embedded NPP Structures*. NUREG/CR-6896.
  https://www.nrc.gov/docs/ML0608/ML060820521.pdf
* NAFEMS (1990). *The Standard NAFEMS Benchmarks*, TNSB Rev. 3; NAFEMS (1987) R0015 *Selected
  Benchmarks for Natural Frequency Analysis*. Targets via the Abaqus Benchmarks Guide §4.4
  (https://ceae-server.colorado.edu/v2016/books/bmk/ch04s04anf18.html etc.) and Altair OptiStruct
  OS-V 0400–0470 (https://help.altair.com/hwsolvers/os/topics/solvers/os/nafems_test_problem_21_r.htm).
