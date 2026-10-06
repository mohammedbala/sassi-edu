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

**Notation.** Matrices are bold in the text where useful and plain in the formula blocks; $^*$ marks
complex (damped) quantities; $^{\mathsf{T}}$ is the transpose (never the conjugate transpose: the SSI
matrices are complex *symmetric*); $i = \sqrt{-1}$; $\omega = 2\pi f$.

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
| Harmonic time factor | $\exp(+i\omega t)$ | R1 §0, D-CNV-01 |
| Horizontal propagation | $\exp(-ikx)$; outgoing waves have $\operatorname{Im} k < 0$, or $\operatorname{Re} k > 0$ when $k$ is real | R1 §0 |
| Cylindrical outgoing waves | Hankel functions of the second kind $H^{(2)}_\mu(k\rho)$ | R1 §3.2 |
| Vertical axis | **$z$ up**; layer interfaces numbered 1 (ground surface) downwards | R1 §0 |
| In-plane thin-layer variables | $\lbrace u_x,\ i\cdot u_z\rbrace$ ("i-scaled vertical", Kausel form): all layer matrices symmetric; physical $u_z = -i\cdot\phi_z$ | R1 §2.2 |
| FFT | forward $A_k = \sum a_n \exp(-2\pi i k n/N)$ (`numpy.fft.rfft`), inverse with $1/N$ (`irfft`) | D-CNV-02 |
| Frequency grid | $\Delta f = 1/(\Delta t\cdot\mathrm{NFFT})$ (or SITE `<fstep>`), SSI frequencies $f = n\cdot\Delta f$ with integer $n$ | requirements §4.0.1 |
| Transfer functions | seismic: $H = U/U_\text{cp}$ (total motion per unit control motion, dimensionless, the same for displacement and acceleration); vibration: displacement per unit load factor | D-CNV-06 |

With $\exp(+i\omega t)$ a time delay $\tau$ is the factor $\exp(-i\omega\tau)$, a viscous dashpot $c$
contributes $+i\omega c$ to the dynamic stiffness and a hysteretic material has a modulus with a
*positive* imaginary part. The transfer functions of ANALYS can therefore be multiplied directly by
`numpy.fft.rfft` of the input.

---

## 2. Complex-modulus damping

**Code:** `sassi/conventions.py` (`cfactor`, `complex_lame`), `sassi/elements/base.py`
(`material_from_M`, `material_from_layer`), `sassi/core/tlm.py` (`layer_moduli`),
`sassi/core/shake.py` (`complex_velocity`). **Verified by:** VP-01, VP-02a, VP-04, VP-39, VP-39b.

### 2.1 The SASSI form

Soil and structural damping are represented as **hysteretic** (rate-independent) damping, through a
complex modulus. The Berkeley codes (SHAKE after the 1973 Udaka-Lysmer revision, SHAKE91, TLUSH,
SASSI2000) use

```math
G^* = G\cdot c(\beta), \qquad c(\beta) = 1 - 2\beta^2 + 2i\beta\cdot\sqrt{1-\beta^2} \tag{2.1}
```

$c(\beta)$ has modulus 1: writing $c = \exp(i\delta)$ gives $\cos\delta = 1 - 2\beta^2$,
$\sin\delta = 2\beta\sqrt{1-\beta^2}$, i.e. $\delta = 2\cdot\arcsin\beta$. Damping rotates the modulus
in the complex plane without changing its magnitude. The complex wave velocity is

```math
V^* = \sqrt{G^*/\rho} = V\cdot\left(\sqrt{1-\beta^2} + i\beta\right), \qquad \lvert V^*\rvert = V \tag{2.2}
```

and a vertically propagating wave $\exp(-ik^*z)$ with
$k^* = \omega/V^* = (\omega/V)\left(\sqrt{1-\beta^2} - i\beta\right)$ decays as
$\exp(-\beta\omega z/V)$: per radian of phase the amplitude decreases by the factor
$\exp\left(-\beta/\sqrt{1-\beta^2}\right) \approx \exp(-\beta)$.

Evidence that SASSI uses this form (R1 §1.1): for the hysteretic SDOF transfer function the recovered
damping $1/(2\lvert H\rvert_{\max})$ is 4.99/9.95/14.83/19.60 % for $\beta = 5/10/15/20\,\%$ with
Eq. (2.1), against 4.98/9.81/14.37/18.57 % for $1 + 2i\beta$; Ostadan et al. (2004) published
5.0/9.9/14.7/19.5 % for SASSI2000. Only Eq. (2.1) reproduces the 19.5 % at $\beta = 20\,\%$.

The simpler form $c = 1 + 2i\beta$ (original SHAKE 1972) is available with `CMODFORM,1` for
benchmarks; its modulus $\sqrt{1+4\beta^2}$ slightly stiffens the material. Damping ratios
$\beta \ge 0.5$ are rejected (EDU-04): the real part $1 - 2\beta^2$ of Eq. (2.1) vanishes at
$\beta = 1/\sqrt{2}$, and such damping ratios have no physical meaning for soils or structures long
before that.

### 2.2 Which modulus carries which damping

The shear modulus carries the S-wave damping and the **constrained** modulus
$M = \lambda + 2G = \rho V_p^2$ carries the P-wave damping (D-CNV-04):

```math
G^* = G\cdot c(\beta_s), \qquad M^* = M\cdot c(\beta_p), \qquad \lambda^* = M^* - 2G^* \tag{2.3}
```

| Element / layer | Complex stiffness |
|---|---|
| soil layers (SITE, POINT), excavated soil, SOLID, PLANE | Lamé constants $(\lambda^*, G^*)$ from Eq. (2.3) |
| BEAMS | $E^*, G^*$ from $(M^*, G^*)$ ($\nu^* = (M^*-2G^*)/\left(2(M^*-G^*)\right)$, $E^* = 2G^*(1+\nu^*)$); $E^* = E\cdot c(\beta)$ when $\beta_p = \beta_s$ |
| SHELL, TSHELL | $E^* = E\cdot c(\beta)$, $\nu$ real (CHECK warns when $\beta_p \ne \beta_s$, EDU-27) |
| SPRING | $k^* = k\cdot c(\text{damp})$ |
| GENERAL | $K^* = \mathrm{MXR} + i\cdot\mathrm{MXI}$ as entered |

Masses are always real. Because every element stiffness is linear in the moduli, an element whose
$\beta_p = \beta_s$ has $K^* = c(\beta)\cdot K_0$.

### 2.3 The hysteretic SDOF and the meaning of β

A mass $m$ on a spring $k^* = k\cdot c(\beta)$ whose base moves with unit amplitude has the
total-motion transfer function

```math
\begin{aligned}
H(\omega) &= \frac{k^*}{k^* - m\omega^2}\\
\lvert H\rvert^2 &= \frac{1}{(1 - 2\beta^2 - r^2)^2 + 4\beta^2(1-\beta^2)}, \qquad r^2 = \frac{m\omega^2}{k}
\end{aligned}
\tag{2.4}
```

The minimum of the denominator is at $r^2 = 1 - 2\beta^2$, which gives the peak

```math
\lvert H\rvert_{\max} = \frac{1}{2\beta\cdot\sqrt{1-\beta^2}} \quad\text{at}\quad \omega = \omega_0\cdot\sqrt{1-2\beta^2} \tag{2.5}
```

(25.0050, 10.0125, 5.0252, 3.3715, 2.5516 for $\beta = 2, 5, 10, 15, 20\,\%$). VP-01 runs exactly
this system through SITE, POINT, HOUSE, ANALYS and MOTION and checks Eq. (2.5) to 1e-6; VP-39/39b do
the same with SPRING and GENERAL elements.

**Equivalent ANSYS input** (R1 §1.2): ANSYS structural damping multiplies the stiffness by
$(1 + ig)$. $E_\text{ANSYS} = E(1 - 2\beta^2)$ and $g = 2\beta\sqrt{1-\beta^2}/(1 - 2\beta^2)$ give
the same complex modulus.

---

## 3. Flexible-volume substructuring

**Code:** `sassi/core/ssi_solver.py`, `sassi/modules/analys.py`, `sassi/modules/house.py`.
**Verified by:** VP-16 (zero-SSI identity), VP-15, VP-E1, VP-S2, the impedance problems VP-10 to
VP-17, and examples 2 and 8 (FV against FI).

### 3.1 The idea

The soil-structure system is split into three parts (Lysmer et al. 1981; Tabatabaie 1982):

```math
\text{SSI system} = \text{free field} + \text{structure} - \text{excavated soil}
```

* the **free field**: the horizontally layered site *without* excavation, which SITE and POINT solve
  semi-analytically (thin-layer method); it is represented by its dynamic stiffness $X_{ff}$ at the
  *interaction nodes* and by its motion $U'_f$ there;
* the **structure**: the finite-element model of the building *including* its basement and any
  near-field soil (HOUSE matrices $K^*_s$, $M_s$);
* the **excavated soil**: a finite-element model of the soil volume that the basement replaces, with
  the free-field layer properties (HOUSE matrices $K^*_e$, $M_e$). It is subtracted, because the free
  field already contains that soil.

### 3.2 Derivation

Let $C(\omega) = K^* - \omega^2 M$ denote a dynamic stiffness. Consider the body "structure minus
excavated soil": its dynamic stiffness is $C^s - C^e$, assembled on one set of nodes (the structure
and the excavation share the nodes of the basement-soil interface). The rest of the world is the free
field with the excavated volume still in place but "cancelled" by the negative stiffness. The free
field has the free-field motion $U'$ when nothing else acts on it; imposing a different displacement
$U$ on its interaction nodes requires the forces $X_{ff}(U - U')$, where $X_{ff}$ is the free-field
dynamic stiffness (impedance) condensed to those nodes. By action and reaction the free field exerts
$-X_{ff}(U - U')$ on the body, and the body's equilibrium (no external force for a seismic analysis)
is

```math
(C^s - C^e)\,U = -X_{ff}(U - U') \quad\Rightarrow\quad \left[C^s - C^e + X_{ff}\right] U = X_{ff}\,U' \tag{3.1}
```

$U$ are **total** displacement amplitudes. With the degree-of-freedom sets of the manual, $s$
structure only, $i$ shared by structure and excavation (basement interface), $w$ excavation only,
and $f$ interaction DOFs, Eq. (3.1) is manual Eq. 2.1 (R1 §4.4):

```math
\begin{bmatrix}
C^s_{ss} & C^s_{si} & 0\\
C^s_{is} & C^s_{ii} - C^e_{ii} + X_{ii} & -C^e_{iw} + X_{iw}\\
0 & -C^e_{wi} + X_{wi} & -C^e_{ww} + X_{ww}
\end{bmatrix}
\begin{bmatrix} U_s\\ U_i\\ U_w \end{bmatrix}
=
\begin{bmatrix} F_s\\ X_{ii}\,U'_i + X_{iw}\,U'_w + F_i\\ X_{wi}\,U'_i + X_{ww}\,U'_w \end{bmatrix}
\tag{3.2}
```

where $X_{ab}$ is the block of $X_{ff}$ (zero when $a$ or $b$ is not an interaction set) and $F$ are
external forces (FORCE; $U' = 0$ for a pure vibration analysis). SASSI-EDU assembles it in the
general form (D-ANL-01)

```math
\begin{aligned}
C(\omega) &= A_s^{\mathsf{T}}(K^*_s - \omega^2 M_s)A_s - A_e^{\mathsf{T}}(K^*_e - \omega^2 M_e)A_e + A_f^{\mathsf{T}} X_{ff}\,A_f\\
b &= A_f^{\mathsf{T}} X_{ff}\,U'_f \quad \text{(seismic)} \qquad\text{or}\qquad b = P \quad \text{(vibration)}
\end{aligned}
\tag{3.3}
```

with Boolean gather matrices $A$ and the load vector $P$ of a vibration analysis from FILE9; the
structure and the excavated soil are on one global DOF map (HOUSE), and $A_f$ picks the translations
of the interaction nodes.

**Zero-SSI identity** (the first built-in check). If the "structure" is identical to the excavated
soil, $C^s = C^e$ and Eq. (3.1) gives $X_{ff}(U - U') = 0$, i.e. $U = U'$ at every node: the free
field is recovered exactly. VP-16 builds such a model (FV, SV, SH and P input) and finds $U = U'$ to
round-off (about 2e-15 relative, tolerance 1e-8).

### 3.3 Method variants: the choice of the interaction set

Eq. (3.1) is exact when **every** node of the excavated soil is an interaction node: the
**flexible-volume (FV)** method. The variants use fewer interaction nodes, which makes $X_{ff}$ (a
dense matrix) much smaller:

| Method | Interaction set $f$ | Consequence |
|---|---|---|
| FV (direct) | all excavated-soil nodes | exact (reference) |
| FI-FSIN, subtraction method (SM) | soil-foundation interface: lateral and bottom faces | the $w$ rows keep only $-C^e_{ww}$; the system is near-singular where $\det C^e_{ww}(\omega) \approx 0$, i.e. at the natural frequencies of the excavated volume with the interface nodes fixed: spurious peaks at and above the first excavated-volume frequency (DOE OE-3 2011-02; Mertz et al. 2011) |
| FI-EVBN, modified subtraction (MSM) | FSIN plus the ground-surface face of the excavation | the surface constraint raises the excavated-volume frequencies; much more accurate, anomalies possible above the new frequency |
| FFV (Ghiocel 2013) | EVBN plus internal horizontal levels every `skip` levels | the sub-volumes between interaction levels are small and stiff; close to FV at a fraction of the cost |

In SASSI-EDU the variants differ **only** by the interaction set (INT/INTGEN); the equations and the
code are the same. Example 2 shows the SM anomaly at 6 Hz on a 10 × 10 × 5 m box and its removal by
FI-EVBN; example 8 shows it at 15.5-16 Hz on a 24 × 24 m shear-wall building embedded 8 m, whose
interior structure has nodes of its own, so that the anomaly comes from the enclosed excavated soil alone
(its first natural frequency with the FI-FSIN nodes fixed is 15.45 Hz).
Screening frequencies used in practice (DOE/STP 2011): the soil-layer frequency $V_s/(4H)$
for embedment $H$ and the excavated-volume frequency $f_\text{EV}$.

---

## 4. The thin-layer method (SITE Mode 1)

**Code:** `sassi/core/tlm.py` (`column_matrices`, `rayleigh_modes`, `love_modes`,
`modal_flexibility`), `sassi/modules/site.py`. **Verified by:** VP-05 (Rayleigh phase velocity),
VP-06 (Love dispersion), VP-43 (monotone convergence), VP-02b, and through POINT VP-07 to VP-09.

### 4.1 Discretisation

The site is a stack of horizontal layers on a rigid base or a visco-elastic half-space. Each layer is
divided into thin **sublayers** in which the displacement varies **linearly** with depth, while the
horizontal dependence stays analytic, $\exp(i\omega t - ikx)$ (Lysmer and Waas 1972; Kausel 1981).
This is a finite-element discretisation in $z$ only, so the accuracy depends on the sublayer
thickness relative to the wavelength: the 1/5-wavelength rule $h \le V_s/(5 f_\text{cut})$ holds for
the **mixed mass** $M = \tfrac{1}{2} M_\text{consistent} + \tfrac{1}{2} M_\text{lumped}$ used by
SASSI-EDU in SITE, POINT and the excavated soil (D-CNV-05). At $h = \lambda/5$ the mixed mass keeps
the 1-D amplification within about 9 %, whereas the consistent mass errs by up to 17 % (R1 §2.1); the
mixed mass makes the shear-wave dispersion error of fourth order (R2 A.3, VP-03).

### 4.2 Layer matrices (Kausel form)

For a sublayer of thickness $h$ with Lamé constants $\lambda^*$, $G^*$ (Eq. 2.3) and density $\rho$,
in the DOF order $\lbrace u_{x1}, u_{y1}, \tilde u_{z1}, u_{x2}, u_{y2}, \tilde u_{z2}\rbrace$ (1 top,
2 bottom interface, $\tilde u_z = i\cdot u_z$, $z$ up), Kausel's Table 1 gives (R1 §2.2; dots are
zeros):

```math
\begin{aligned}
A_m &= \frac{h}{6}\begin{bmatrix}
2(\lambda+2G) & \cdot & \cdot & \lambda+2G & \cdot & \cdot\\
\cdot & 2G & \cdot & \cdot & G & \cdot\\
\cdot & \cdot & 2G & \cdot & \cdot & G\\
\lambda+2G & \cdot & \cdot & 2(\lambda+2G) & \cdot & \cdot\\
\cdot & G & \cdot & \cdot & 2G & \cdot\\
\cdot & \cdot & G & \cdot & \cdot & 2G
\end{bmatrix}\\[8pt]
B_m &= \frac{1}{2}\begin{bmatrix}
\cdot & \cdot & \lambda-G & \cdot & \cdot & -(\lambda+G)\\
\cdot & \cdot & \cdot & \cdot & \cdot & \cdot\\
\lambda-G & \cdot & \cdot & \lambda+G & \cdot & \cdot\\
\cdot & \cdot & \lambda+G & \cdot & \cdot & -(\lambda-G)\\
\cdot & \cdot & \cdot & \cdot & \cdot & \cdot\\
-(\lambda+G) & \cdot & \cdot & -(\lambda-G) & \cdot & \cdot
\end{bmatrix}\\[8pt]
G_m &= \frac{1}{h}\begin{bmatrix}
G & \cdot & \cdot & -G & \cdot & \cdot\\
\cdot & G & \cdot & \cdot & -G & \cdot\\
\cdot & \cdot & \lambda+2G & \cdot & \cdot & -(\lambda+2G)\\
-G & \cdot & \cdot & G & \cdot & \cdot\\
\cdot & -G & \cdot & \cdot & G & \cdot\\
\cdot & \cdot & -(\lambda+2G) & \cdot & \cdot & \lambda+2G
\end{bmatrix}
\end{aligned}
```

```math
\begin{aligned}
M_m^\text{cons} &= \frac{\rho h}{6}\begin{bmatrix} 2I & 1I\\ 1I & 2I \end{bmatrix}, \qquad
M_m^\text{lump} = \frac{\rho h}{2}\begin{bmatrix} I & 0\\ 0 & I \end{bmatrix}\\[4pt]
M_m &= \tfrac{1}{2} M_m^\text{cons} + \tfrac{1}{2} M_m^\text{lump}
= \rho h\begin{bmatrix} \tfrac{5}{12} I & \tfrac{1}{12} I\\[2pt] \tfrac{1}{12} I & \tfrac{5}{12} I \end{bmatrix}
\end{aligned}
```

The sublayer matrices are assembled by overlapping their blocks at shared interfaces, exactly as in
finite-element assembly, and give the interface load-displacement relation in the
frequency-wavenumber domain

```math
\bar P = \left(A k^2 + B k + G - \omega^2 M\right)\bar U \tag{4.1}
```

The direction sub-blocks used below are $A_x$ (with $\lambda+2G$), $A_z = A_y$ (with $G$),
$G_x = G_y$ (with $G$), $G_z$ (with $\lambda+2G$), the $x$-row/$z$-column coupling block $B_{xz}$, and
$C_\bullet = G_\bullet - \omega^2 M$. The i-scaled vertical variable makes every matrix complex
symmetric; the Waas/SASSI form (vertical displacement down, no $i$, a skew-symmetric $B$ entering as
$iBk$) is algebraically identical ($T^{-1} L_\text{Waas}\,T = L_\text{Kausel}$ with
$T = \operatorname{diag}(1, i, \ldots)$, R1 §2.3).

A rigid base deletes the bottom-interface DOFs; a half-space adds the dashpots of section 5.

### 4.3 Eigenproblems

Setting $\bar P = 0$ gives the wavenumbers $k$ of the free waves of the discrete layered medium at a
given frequency.

**Rayleigh (P-SV, in-plane).** The quadratic eigenproblem $(A k^2 + B k + C)\,\phi = 0$ is
linearised *without* doubling its size by the substitution $Z = \lbrace\phi_x;\ k\,\phi_z\rbrace$
(Kausel 1981, Eq. 16):

```math
\left(k^2\begin{bmatrix} A_x & 0\\ B_{xz}^{\mathsf{T}} & A_z \end{bmatrix}
+ \begin{bmatrix} C_x & B_{xz}\\ 0 & C_z \end{bmatrix}\right)
\begin{Bmatrix} \phi_x\\ k\,\phi_z \end{Bmatrix} = 0 \tag{4.2}
```

a generalised linear eigenproblem in $k^2$ of size $2N_f$ ($N_f$ free interfaces). SASSI-EDU solves
it as the standard eigenproblem of $-\bar A^{-1}\bar C$ ($\bar A$ is block lower triangular, so the
inverse is two symmetric solves). Each root $k_j = \sqrt{k_j^2}$ is taken on the branch with
**$\operatorname{Im} k_j < 0$**, or $\operatorname{Re} k_j > 0$ when $k_j$ is real (decaying or
outgoing in $+x$), and $\phi_z = Z_\text{lower}/k_j$.

**Love (SH, anti-plane).** $(A_y k^2 + C_y)\,\phi_y = 0$, $N_f$ modes.

**Normalisation** (Kausel 1981, Eq. 22a; the same as Waas'):

```math
\begin{aligned}
\phi_x^{\mathsf{T}} A_x \phi_x + \phi_z^{\mathsf{T}} A_z \phi_z + \phi_z^{\mathsf{T}} B_{xz}^{\mathsf{T}} \phi_x / k_j &= 1 &&\text{(Rayleigh)}\\
\phi_y^{\mathsf{T}} A_y \phi_y &= 1 &&\text{(Love)}
\end{aligned}
\tag{4.3}
```

With this normalisation the modes are orthogonal ($Y_i^{\mathsf{T}} \bar A Z_j = 0$ for $i \ne j$) and
the wavenumber-domain flexibility has the modal expansion (Kausel 1981, Eqs. 49-52)

```math
\begin{aligned}
&F_{xx} = \Phi_x D \Phi_x^{\mathsf{T}}, \qquad F_{xz} = k\,\Phi_x K^{-1} D \Phi_z^{\mathsf{T}}, \qquad
F_{zz} = \Phi_z D \Phi_z^{\mathsf{T}}, \qquad F_{yy} = \Phi_y D_L \Phi_y^{\mathsf{T}},\\
&D = \operatorname{diag}\left(\frac{1}{k^2 - k_j^2}\right), \qquad K = \operatorname{diag}(k_j)
\end{aligned}
\tag{4.4}
```

which SASSI-EDU checks against the direct inverse of Eq. (4.1) to round-off (unit tests of
`tlm.modal_flexibility`, R1 V1). The modes are sorted by decreasing $\operatorname{Re} k$ (shortest
wavelength first). FILE2 stores, per frequency, the generated sublayers and
$\lbrace k_j, \phi_x, \phi_z\rbrace$ ($2N_f$ Rayleigh modes) and $\lbrace k_l, \phi_y\rbrace$ ($N_f$
Love modes); one FILE2 serves POINT2 and POINT3, because the cylindrical and the plane-strain
eigenproblems are identical.

**What VP-05 and VP-06 show.** The fundamental Rayleigh velocity of a deep uniform stratum converges
at $O(h^2)$ to the exact $c_R/V_s$ (0.9194 at $\nu = 0.25$) within 0.5 % at 10 sublayers per
wavelength for $\nu \le 1/3$; for $\nu \ge 0.45$ linear thin layers lock volumetrically and need about
40 sublayers per wavelength (lead decision D-W1-02; CHECK warns for $\nu > 0.47$). The Love dispersion
of a layer over a half-space is reproduced within 0.5 %, but the variable-depth half-space cannot
represent the grazing field near a Love cut-off frequency (D-W1-03).

---

## 5. Half-space: variable depth and viscous boundary

**Code:** `sassi/core/tlm.py` (`halfspace_depth`, `halfspace_sublayers`, `build_column`,
`vertical_response`, `outcrop_response`). **Verified by:** VP-02b, VP-08, VP-10, VP-11 (radiation
damping), UT-19.

A visco-elastic half-space under the user layers is simulated with two techniques together
(manual §4.1.2 item 20, R1 §2.4):

1. **Variable depth.** SITE adds $n_l$ sublayers (SITE `<nl>`, 4-20, 0 = rigid base) with the
   half-space properties and a total thickness of one and a half shear wavelengths of the half-space,

   ```math
   H_\text{hs}(f) = 1.5\,\lambda_s = 1.5\,V_{s,\text{hs}}/f \tag{5.1}
   ```

   because the fundamental Rayleigh mode has decayed at that depth. The buffer therefore becomes
   deeper at low frequency and thinner at high frequency.
2. **Viscous boundary.** Lysmer-Kuhlemeyer dashpots per unit area are attached at the bottom of the
   buffer: $c_s = \rho V^*_s$ on the horizontal DOFs and $c_p = \rho V^*_p$ on the vertical one, with
   the *complex* velocities of Eq. (2.2), so that for vertically propagating waves the dashpot is the
   exact impedance of the damped half-space. With $\exp(+i\omega t)$ a dashpot adds $+i\omega c$ to
   the diagonal of $G$ at the base interface (in the i-scaled variables the vertical entry is
   unchanged because the $i$ factors cancel).

**Sublayer law.** The manual says only that the thicknesses increase with depth. SASSI-EDU offers
`EDUOPT,HSLAW,UNIFORM` ($h_i = H_\text{hs}/n_l$, the default), `GEOMETRIC` ($h_i = h_1 q^{i-1}$,
$h_1 = \min(h_\text{last}\cdot V_{s,\text{hs}}/V_{s,\text{last}},\ H_\text{hs}/n_l)$, $q$ found by
bisection so that $\sum h_i = H_\text{hs}$) and `LINEAR`; the non-uniform laws fall back to uniform
when a sublayer would exceed $\lambda_s/8$. VP-02b compared the outcrop transfer function of a uniform
layer on a half-space with the exact solution: the geometric law missed the 1 % criterion (1.10 % at
$n_l = 20$) and the uniform law passed (0.20 % / 0.95 % at $n_l$ = 20 / 10), so the default is
UNIFORM (lead decision D-W1-01, the fallback written into decision D-SIT-02). R1 V4 had found
≤ 0.7 % ($n_l = 20$) and ≤ 1.5 % ($n_l = 10$) with uniform sublayers, but up to 20 % with
$n_l = 5$.

**Within and outcrop motions.** The *within* transfer function between two interfaces above the
half-space does not depend on the half-space model; the half-space model matters for the motion
*below* the control point and for the *outcrop* motion. SASSI-EDU computes the discrete outcrop
motion at the top of the half-space as the surface motion of the half-space-only column (same
generated sublayers, dashpot and incident wave), R1 §2.7a.

---

## 6. Free-field motion and control-point normalisation (SITE Mode 2)

**Code:** `sassi/core/freefield.py`, `sassi/core/tlm.py` (`vertical_response`, `choose_mode`),
`sassi/modules/site.py`. **Verified by:** VP-02b, VP-15, VP-16, VP-24, VP-S2, example 4.

The seismic environment is a superposition of plane **wave fields** $w$: body waves (P, SV, SH) at
vertical or inclined incidence, and surface waves (Rayleigh, Love). For each wave and frequency SITE
computes the field at every user interface in the SITE axes $x'y'z'$ ($z'$ up, $x'$ in the vertical
plane of propagation, $y' = z' \times x'$).

### 6.1 Vertically propagating body waves

With $k = 0$ the layered column decouples into 1-D columns: SV ($x'$) and SH ($y'$) use the shear
matrices ($G_s$, built with $G^*$), P ($z'$) the compression matrices ($G_p$, built with
$M^* = \lambda^* + 2G^*$). For an incident (up-going) wave of displacement amplitude $E_b$ in the
half-space, the traction of the half-space on the column base is
$t = i\omega\rho V^*(u_I - u_R) = i\omega\rho V^*(2u_I - u)$, which gives (R1 §2.7a)

```math
\begin{aligned}
\left(G_\bullet - \omega^2 M + i\omega c\,e_N e_N^{\mathsf{T}}\right) u &= 2i\omega c\,E_b\,e_N &&\text{(dashpot base)}\\
u_N &= 1 &&\text{(rigid base)}
\end{aligned}
\tag{6.1}
```

with $E_b = 1$ for the dashpot base.

### 6.2 Normalisation to the control point

The control motion is specified at the **top of the control layer** `<cl>` in the direction `<cm>`
($x'$, $y'$ or $z'$), as the **within** (in-column) motion (D-SIT-07). Each wave field is divided by
its own motion there:

```math
\hat U_w(z) = \frac{u_w(z)}{u_{w,\text{cm}}(z_\text{cp})} \tag{6.2}
```

so the control point moves with unit amplitude; this ratio of displacements is also the ratio of
accelerations. With several waves the total field is

```math
U(z, x') = \sum_w r_w(f)\,\hat U_w(z)\,\exp\left(-i k_w (x' - x'_\text{cp})\right) \tag{6.3}
```

with participation ratios $r_w$ given at two frequency numbers (Frequency 1 and 2), interpolated
linearly in between, held constant outside, each in $(0, 1]$ and summing to 1 (Errors 53-55). Because
each field is unit-normalised and the ratios sum to one, the control point has unit motion. When the
design motion is an *outcrop* motion, the user runs SOIL first and uses the computed within motion
(R1 §2.7a; examples 4 and 6).

### 6.3 Inclined body waves and surface waves

* **Inclined P, SV, SH** (D-SIT-08): apparent wavenumber $k = \omega\sin\theta/V_\text{hs}$ ($V = V_s$
  for SV/SH, $V_p$ for P; $\theta$ from the vertical). The generated buffer is not used: the user
  column sits on the **exact** half-space stiffness at wavenumber $k$ (from P and SV potentials,
  Kausel and Roesset 1981), $K_\text{out}$ for the radiating down-going field and $K_\text{in}$ for
  the incident up-going field, and
  ```math
  \left(A k^2 + B k + G - \omega^2 M + K_\text{out}\,e_N e_N^{\mathsf{T}}\right) U = (K_\text{out} - K_\text{in})\,u_I \tag{6.4}
  ```
  with $u_I$ the displacement of the pure incident wave at the top of the half-space. For
  $\theta = 0$, $K_\text{out} = -K_\text{in} = i\omega\rho V^*$ and Eq. (6.4) reduces to Eq. (6.1).
* **Surface waves**: one mode of FILE2, chosen by WAVE `<opt>`: 1 "shortest wavelength" = the
  propagating mode with the largest $\operatorname{Re} k$, 2 "least decay" = the mode with the
  smallest $\lvert\operatorname{Im} k\rvert$ (ties within $10^{-8}\max\lvert k\rvert$ broken by the
  largest $\operatorname{Re} k$). A mode is *propagating* when
  $\lvert\arg k\rvert \le \arctan(0.5) + \arctan(\text{loss})$, the elastic sector widened by the
  material loss angle of the column, $\text{loss} = \tan(\delta/2)$ (D-SIT-09, D-W1-08). The field is
  $\phi_j(z)\exp\left(-ik_j(x' - x'_\text{cp})\right)$ with the physical vertical component
  $u_z = -i\phi_z$, normalised by Eq. (6.2).

The manual's definitions of "shortest wavelength" and "least decay" and the treatment of the
half-space for inclined waves are reconstructions (R1 §10 items 1 and 3).

### 6.4 Evaluation at the interaction nodes

FILE1 holds, per wave and frequency, the unit-normalised field at the user interfaces, the
wavenumber and the ratio. ANALYS evaluates it at an interaction node on interface $m$ at $(x, y)$:

```math
\begin{aligned}
x' &= (x - x_c)\cos a + (y - y_c)\sin a\\
U'(\text{node}) &= R_z(a)\sum_w r_w\,U_w(m)\,\exp(-i k_w x')
\end{aligned}
\tag{6.5}
```

where $a$ is the coordinate transformation angle (ANALYS `<ang>`, from $x'$ to the global $x$ axis,
counter-clockwise about $z$) and $(x_c, y_c)$ the control point, which has $x' = 0$. Vertical
incidence has $k = 0$: every node of an interface moves in phase. VP-24 checks the rotation on a
rotationally symmetric model ($\mathrm{ATF}_x = \cos a\cdot\mathrm{ATF}(a=0)$,
$\mathrm{ATF}_y = \sin a\cdot\mathrm{ATF}(a=0)$).

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
response to a load at depth $n$ depends only on the horizontal distance to the observation point:
POINT solves one load per **load interface** $n = 1 \ldots L+1$ ($L$ = POINT `<layer>`, the number
of embedment layers) on the axis of an axisymmetric model, and ANALYS translates and rotates the
solution to every pair of interaction nodes (section 8).

The load is expanded in Fourier harmonics about the vertical axis: a vertical load excites the
harmonic $\mu = 0$, a horizontal load $\mu = 1$, with the displacement pattern (symmetric case,
R1 §3.2)

```math
u_\rho = \tilde u_\rho\cos\mu\theta, \qquad u_\theta = -\tilde u_\theta\sin\mu\theta, \qquad u_z = \tilde u_z\cos\mu\theta \tag{7.1}
```

A true point load gives an infinite displacement under it. SASSI therefore spreads it over a small
**central zone** of radius $R_0$ (POINT `<rad>`), a finite-element core whose size calibrates the
diagonal (self) terms of the flexibility to the mesh of interaction nodes (section 7.5).

### 7.2 The core

One radial axisymmetric finite element (nodes on the axis and at $\rho = R_0$) per SITE sublayer,
including the generated half-space sublayers and the base dashpots (decision D-PNT-03; the
discretisation of the original "special cylindrical axisymmetric elements" is a reconstruction,
R1 §3.3). Bilinear interpolation $L_0 = 1 - \rho/R_0$, $L_1 = \rho/R_0$ in $\rho$ and the thin-layer
linear interpolation in $z$; strains

```math
\begin{aligned}
\varepsilon = \Big[\,&\frac{\partial\tilde u_\rho}{\partial\rho},\ \ \frac{\tilde u_\rho - \mu\tilde u_\theta}{\rho},\ \ \frac{\partial\tilde u_z}{\partial z},\ \ \frac{\partial\tilde u_\rho}{\partial z} + \frac{\partial\tilde u_z}{\partial\rho},\\
&\frac{\mu\tilde u_\rho}{\rho} + \frac{\partial\tilde u_\theta}{\partial\rho} - \frac{\tilde u_\theta}{\rho},\ \ \frac{\partial\tilde u_\theta}{\partial z} + \frac{\mu\tilde u_z}{\rho}\,\Big]^{\mathsf{T}}
\end{aligned}
\tag{7.2}
```

with the isotropic $D(\lambda^*, G^*)$, $K_c = \iint B^{\mathsf{T}} D B\,\rho\,d\rho\,dz$, 3 × 3
Gauss points ($\rho = 0$ is never evaluated), the mixed mass in $z$ (as in SITE) and the consistent
mass in $\rho$. Axis constraints: $\mu = 0 \to \tilde u_\rho = \tilde u_\theta = 0$;
$\mu = 1 \to \tilde u_\rho = \tilde u_\theta$ (one tied unknown), $\tilde u_z = 0$. A viscous base
adds $i\omega c\int L_a L_b\,\rho\,d\rho$ on the bottom interface.

### 7.3 The exterior and the transmitting boundary

For $\rho \ge R_0$ the field is a superposition of outgoing Rayleigh and Love modes of FILE2:

```math
\begin{aligned}
\tilde u(\rho) &= \Psi_\mu(\rho)\,\alpha, \qquad \alpha = \left[a_1 \ldots a_{2N_f};\ b_1 \ldots b_{N_f}\right]\\[6pt]
\Psi_\mu(\rho) &= \begin{bmatrix}
\Phi_x \operatorname{diag} H'_\mu(k^R\rho) & \Phi_y \operatorname{diag} \dfrac{\mu H_\mu(k^L\rho)}{k^L\rho}\\[8pt]
\Phi_x \operatorname{diag} \dfrac{\mu H_\mu(k^R\rho)}{k^R\rho} & \Phi_y \operatorname{diag} H'_\mu(k^L\rho)\\[8pt]
-\Phi_z \operatorname{diag} H_\mu(k^R\rho) & 0
\end{bmatrix}
\end{aligned}
\tag{7.3}
```

with $H_\mu = H^{(2)}_\mu$ and $H'_\mu = dH_\mu(x)/dx$ (the columns of Kausel's $C_\mu$ operator with
$J \to H^{(2)}$). The consistent generalised forces that the exterior exerts on the core at
$\rho = R_0$ (per radian) are

```math
\begin{aligned}
f_\rho &= R_0\left[E_{\lambda+2G}\,\frac{\partial\tilde u_\rho}{\partial\rho} + E_\lambda\,\frac{\tilde u_\rho - \mu\tilde u_\theta}{R_0} + Q_\lambda\,\tilde u_z\right]\\
f_\theta &= R_0\left[E_G\left(\frac{\mu\tilde u_\rho}{R_0} + \frac{\partial\tilde u_\theta}{\partial\rho} - \frac{\tilde u_\theta}{R_0}\right)\right]\\
f_z &= R_0\left[Q_G\,\tilde u_\rho + E_G\,\frac{\partial\tilde u_z}{\partial\rho}\right]
\end{aligned}
\tag{7.4}
```

with the thin-layer boundary integrals $E_c = \int c\,N N^{\mathsf{T}}\,dz$ (per sublayer
$c\,h/6\,[2\ 1;\ 1\ 2]$) and $Q_c = \int c\,N N'^{\mathsf{T}}\,dz$ (per sublayer
$c/2\,[1\ {-1};\ 1\ {-1}]$). Collecting them per mode column as $T_\mu$, the **transmitting
boundary** (consistent boundary of Waas, Lysmer and Kausel) is

```math
R_\mu = -T_\mu\,\Psi_\mu(R_0)^{-1} \tag{7.5}
```

a full $3N_f \times 3N_f$ matrix. Numerically, the Hankel functions are evaluated with the
exponentially scaled `hankel2e` and the basis is scaled to $R_0$ (column $j$ multiplied by
$\exp(ik_j R_0)$), which avoids overflow for $\rho \ge R_0$ and underflow of evanescent modes;
$R_\mu$ does not depend on that scaling.

### 7.4 Solution

For each load interface $n$ (unit load on the axis):

```math
\begin{aligned}
&\left[K_c - \omega^2 M_c + R_\mu\right]\tilde U = \tilde F, \qquad
\tilde F = P/\pi \ \ (\mu = 1), \qquad \tilde F = P/(2\pi) \ \ (\mu = 0)\\
&\alpha = \Psi_\mu(R_0)^{-1}\,\tilde U_\text{boundary}, \qquad \tilde u(\rho) = \Psi_\mu(\rho)\,\alpha \quad \text{for } \rho \ge R_0
\end{aligned}
\tag{7.6}
```

The load of $\mu = 1$ acts on the tied axis unknown. ($\pi$ and $2\pi$ are the angular integrals of
$\cos^2$, $\sin^2$ and 1.) FILE3 stores, per frequency, load interface and $\mu$, the amplitudes
$\alpha$, the core-axis displacements at interfaces $1 \ldots L+1$ and the mode rows at the
observation interfaces, so that ANALYS evaluates the field **exactly at any distance** without
tabulation (D-PNT-02).

**Exact reference.** For a point load in a thin-layer medium Kausel (1981, Eqs. 88-89) gives the
Green functions as Hankel series over the same modes, e.g. for a vertical unit load
$u_z = \tfrac{1}{4i}\sum_l \phi_z^{ml}\phi_z^{nl} H_0(k_l^R\rho)$. They are implemented in
`sassi/core/greens_tlm.py` and are the reference of VP-09: the POINT3 far field agrees within 2.35 %
at $3.3\,R_0$ (lead decision D-W1-04), 0.8 % at $6.7\,R_0$ and 0.25 % beyond $13\,R_0$ on the dominant
components (R1 V8). In the static limit on a deep graded stratum the Green functions match
Boussinesq and Cerruti within 0.1-1.7 % (VP-07, R1 V3).

### 7.5 The central-zone radius

With one radial element, every pair of interaction nodes at a distance $r \ge R_0$ lies in the exact
far field; $R_0$ only calibrates the near-singular diagonal terms. The rules of the manual are
$R_0 = 0.90\,h$ for square meshes, $0.85\,h$ for triangular (circular) meshes and $R_0 = h$ in 2D
($h$ = mesh size). The self-flexibility of the core is much softer than that of a uniform disk of
radius $R_0$ (equivalent disk radius 0.25-0.29 $R_0$, R1 V9), so a "tributary disk" model would not
reproduce SASSI's diagonal terms. For non-uniform meshes `RADIUS` computes
$r_e = \mathrm{Scale}\cdot\sqrt{A_\text{plan}}$ per excavated element and the average (D-PNT-05).

### 7.6 POINT2 (plane strain)

The 2D counterpart: a strip $\lvert x\rvert \le R_0$ of two plane-strain elements per sublayer, with
Waas-Lysmer transmitting boundaries at $x = \pm R_0$ built from the same modes (outgoing mode $j$:
$u_x = \phi_{xj}\exp\left(-isk_j(x - x_b)\right)$, $u_z = -is\,\phi_{zj}\exp(\ldots)$, $s = \pm 1$),
$R = -T\,\Psi^{-1}$; the far field is
$\sum\alpha_j\psi_j\exp\left(-ik_j(\lvert x\rvert - R_0)\right)$, the vertical response to a
horizontal load being odd in $x$. The reference is Kausel's line-load
Green function (R1 V6, V7); VP-T1 checks the FILE3 far field against an independent modal-sum
implementation of it: within 0.33 % at $\lvert x\rvert = 2R_0$ and 0.3 % farther out (lead decision
D-W3-12). HOUSE `<dim>` = 1 selects POINT2, and ANALYS solves 2D (plane-strain) models with it
(section 15).

---

## 8. Flexibility and impedance assembly

**Code:** `sassi/core/flexibility.py` (`flexibility_matrix`, `cylindrical_components`),
`sassi/core/ssi_solver.py` (`impedance_matrix`, `rigid_body_transform`, `global_impedance`).
**Verified by:** VP-07 (reciprocity), VP-10, VP-11, VP-13, VP-14, VP-17.

### 8.1 The 3 × 3 blocks

For interaction node $i$ on interface $m$ at $(x_i, y_i)$ and node $j$ on interface $n$, let
$r = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}$ and $\theta = \operatorname{atan2}(y_i - y_j,\ x_i - x_j)$
(from the load point $j$ to the observation point $i$). From POINT (load at $n$, observed at $m$)
take $\tilde u, \tilde v, \tilde w = (\tilde u_\rho, \tilde u_\theta, \tilde u_z)$ for $\mu = 1$ and
$\tilde p, \tilde q = (\tilde u_\rho, \tilde u_z)$ for $\mu = 0$. Rotating the $\cos\theta$ /
$-\sin\theta$ pattern of Eq. (7.1) to Cartesian axes gives (rows $u_x, u_y, u_z$ at $i$; columns
$P_x, P_y, P_z$ at $j$; R1 §4.1):

```math
F_{ij} = \begin{bmatrix}
\tilde u c^2 + \tilde v s^2 & (\tilde u - \tilde v)\,sc & \tilde p c\\
(\tilde u - \tilde v)\,sc & \tilde u s^2 + \tilde v c^2 & \tilde p s\\
\tilde w c & \tilde w s & \tilde q
\end{bmatrix}, \qquad c = \cos\theta, \quad s = \sin\theta \tag{8.1}
```

* $r = 0$ (same vertical line):
  $\operatorname{diag}(\tilde u_\text{axis}, \tilde u_\text{axis}, \tilde q_\text{axis})$ from the
  core-axis values;
* $0 < r < R_0$ (mesh finer than the central zone): linear in $r$ between the axis value and the
  value at $R_0$, consistent with the core shape functions;
* $r \ge R_0$: the exact expansion $\Psi_\mu(r)\,\alpha$. When a model has more than 6,000 distinct
  distances, the expansion is tabulated on a dense grid (geometric near $R_0$, spacing at most
  $\lambda_{\min}/48$) and interpolated by cubic splines of $r\cdot u(r)$.

**Reciprocity.** $F_{ij} = F_{ji}^{\mathsf{T}}$ holds exactly for the exact Green functions and to
about 1e-4 with the finite-element core (R1 V10, VP-07 checks
$\lVert F - F^{\mathsf{T}}\rVert/\lVert F\rVert < 10^{-3}$). SASSI-EDU symmetrises
$F \leftarrow \tfrac{1}{2}(F + F^{\mathsf{T}})$ (D-ANL-02), because SASSI treats the impedance as
symmetric.

**2D.** For POINT2 data the blocks are the 2 × 2 P-SV and 1 × 1 SH blocks with the odd coupling terms
multiplied by $\operatorname{sgn}(x_i - x_j)$; ANALYS inverts the in-plane (UX, UZ) part only
(section 15.1). For a half or quarter model (SYMM) the blocks of the image nodes are added before the
inversion, Eq. (15.2).

### 8.2 Impedance

```math
X_{ff}(\omega) = F_{ff}(\omega)^{-1} \tag{8.2}
```

which is manual Eq. 4.3, $K + i\omega D = (f + ig)^{-1}$: a full, complex symmetric
$3n_f \times 3n_f$ matrix, computed by a dense LU factorisation and inversion in place (LAPACK
getrf/getri; never Cholesky). Its size dominates the memory of large models: 16 bytes per entry,
$(3n_f)^2 \times 16$ bytes per matrix (14.4 GB for 10,000 interaction nodes).

### 8.3 Global (unconstrained) foundation impedance

With ANALYS `<impe>` > 0 the impedance is condensed to the six rigid-body DOFs of the interaction
nodes about the control point (D-ANL-07):

```math
\begin{aligned}
K_G(\omega) &= T^{\mathsf{T}} X_{ff}\,T, \qquad u_j = u_0 + \theta \times d_j, \qquad d_j = r_j - (x_c, y_c, z_c)\\[6pt]
T_j &= \begin{bmatrix}
1 & 0 & 0 & 0 & d_z & -d_y\\
0 & 1 & 0 & -d_z & 0 & d_x\\
0 & 0 & 1 & d_y & -d_x & 0
\end{bmatrix}
\end{aligned}
\tag{8.3}
```

and written as FOUNSTIF ($\operatorname{Re} K_G$), FOUNDASH ($\operatorname{Im} K_G/\omega$), FOUNDAMP
($\operatorname{Im} K_G/(2\lvert\operatorname{Re} K_G\rvert)$) and FOUNIMPD ($\lvert K_G\rvert$). This
is the impedance of the soil at the interaction nodes when they are forced to move as a
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

For every SSI frequency $f_q = n_q\,\Delta f$:

1. **Frequency survey**: the frequency numbers must exist in FILE1 (or FILE9) and in FILE3, otherwise
   the run stops before any solve.
2. **Flexibility** $F_{ff}$ from FILE3 (section 8.1) and **impedance** $X_{ff} = F_{ff}^{-1}$.
3. **Dynamic matrix** $C = (K^*_s - \omega^2 M_s) - (K^*_e - \omega^2 M_e)$ (sparse, HOUSE matrices)
   plus $A_f^{\mathsf{T}} X_{ff}\,A_f$.
4. **Load**: seismic $b_f = X_{ff}\,U'_f$ with the free field of Eq. (6.5) at each interaction node
   (incoherent motion, wave passage and multiple excitation multiply the load or the motion by the
   HOUSE factors of FILE77, Eq. 14.7); vibration $b$ = the FORCE vector,
   $P_k(\omega) = a_k\exp(-i\omega t_k)$ for a load factor $a_k$ and arrival time $t_k$ on the
   reference load history (a delay is a phase lag; moving loads are loads with increasing arrival
   times, VP-40).
5. **Solve** (section 9.2) and write the transfer function of every equation to FILE8 (fixed DOFs:
   $H = 0$).

### 9.2 Schur complement

$X_{ff}$ is dense on the interaction DOFs $f$, everything else is sparse. With the remaining DOFs $n$:

```math
\begin{bmatrix} C_{nn} & C_{nf}\\ C_{fn} & C_{ff} + X_{ff} \end{bmatrix}
\begin{bmatrix} U_n\\ U_f \end{bmatrix} = \begin{bmatrix} b_n\\ b_f \end{bmatrix}
```

```math
\begin{aligned}
S &= C_{ff} - C_{fn}\,C_{nn}^{-1}\,C_{nf}\\
U_f &= (S + X_{ff})^{-1}\left(b_f - C_{fn}\,C_{nn}^{-1}\,b_n\right)\\
U_n &= C_{nn}^{-1}\left(b_n - C_{nf}\,U_f\right)
\end{aligned}
\tag{9.1}
```

with $S$ the dense Schur complement. $C_{nn}$ is factorised by a sparse LU (SuperLU) and
$S + X_{ff}$ by a dense LU, both complex symmetric but solved with LU (never Cholesky or a Hermitian
solver). Every right-hand side (simultaneous cases, restarts) reuses both factorisations. If $C_{nn}$
alone is singular or ill-conditioned (estimated reciprocal condition number below 1e-12, for example
an undamped structure with clamped interaction DOFs at one of its fixed-base frequencies) the full
system with the dense block inserted is factorised by a sparse LU instead, with a warning. The cost
per frequency is dominated by the $O(n_f^3)$ dense operations.

### 9.3 Restarts

| ANALYS `<mode>` | Name | Stored by the initiation run | Re-used |
|---|---|---|---|
| 1 | New Structure | `COOXqqq` = $X_{ff}$ per frequency | $X_{ff}$; the structure is re-assembled and re-factorised |
| 2 | New Seismic Environment | `COOTKqqq` = the factorised system (sparse LU of $C_{nn}$, Schur data, dense LU of $S + X_{ff}$) | $X_{ff}$ and the factorisation; only $b$ changes |
| 3 | New Dynamic Loading | same as 2 | same as 2 |

Each record carries the frequency number and value and the FILE90 hashes of the interaction nodes,
the layering and (COOTK) the structure matrices; a restart refuses records that do not match, and
`COOXI`/`COOTKI` index the records by frequency number so that any subset of the saved frequencies
can be re-run. A restart gives the same FILE8 as a fresh initiation run to round-off (VP-22,
tolerance 1e-10).

### 9.4 Simultaneous cases and the low-frequency check

With `<simul>` = 1 the three coherent input directions (FILE1X/Y/Z, angle 0) are three right-hand
sides of one factorisation and give FILE8X/Y/Z; with `<simul>` = $N_l \ge 2$ in a vibration analysis
the $N_l$ load cases FILE9001 ... give FILE8001 ... VP-23 checks that each case equals a single run to
1e-12. An incoherent stochastic simulation with $N_s$ samples (`<simul>` = $N_s \le 50$) solves the
$3N_s$ right-hand sides FILE1X/Y/Z × FILE77001 ... FILE77Ns and writes `FILE8{3(s-1)+d}` (sample $s$,
direction $d = 1, 2, 3$; section 14.4). In 2D, `<simul>` = 1 analyses the in-plane X and Z cases only.

As the frequency tends to zero the whole system moves rigidly with the free field, so every
transfer function in the input direction tends to 1. ANALYS reports the deviation at the first SSI
frequency (G-19); VP-41 checks it on several models.

---

## 10. Transfer-function interpolation

**Code:** `sassi/core/interp.py` (`interpolate_tf`, `smooth_tf`, `phase_adjust`), used by MOTION and
STRESS. **Verified by:** VP-28 (exactness, ISRS from interpolated TFs), VP-32 (SRSS TF), VP-33 (phase
adjustment), VP-53 (CRITFREQ).

ANALYS solves at a few tens to a few hundred SSI frequencies; the convolution needs the transfer
function at every Fourier frequency $f_k = k\,\Delta f$ up to the cut-off.

### 10.1 The SASSI (Tajirian) interpolant

Over a short band, every SSI transfer function behaves like that of a **two-degree-of-freedom system
with hysteretic damping**, whose transfer function is a ratio of quadratics in $\omega^2$ (no odd
powers, because hysteretic damping does not depend on frequency; Tajirian 1981, TLUSH 1981):

```math
H(\omega) = \frac{C_1\omega^4 + C_2\omega^2 + C_3}{\omega^4 + C_4\omega^2 + C_5} \tag{10.1}
```

The five complex constants are fitted to five computed values $H_p = H(\omega_p)$ of a **window** of
five consecutive SSI frequencies:

```math
\begin{bmatrix} \omega_p^4 & \omega_p^2 & 1 & -\omega_p^2 H_p & -H_p \end{bmatrix}\cdot
\begin{Bmatrix} C_1 & C_2 & C_3 & C_4 & C_5 \end{Bmatrix}^{\mathsf{T}} = \omega_p^4 H_p, \qquad p = 1, \ldots, 5 \tag{10.2}
```

(the MHI transcription prints $+H_p$ in the last column; moving the denominator of Eq. 10.1 across
gives $-H_p$, and only that sign reproduces an exact 2-DOF transfer function, R1 §5.1). SASSI-EDU
uses the scaled variable $x = (\omega/\omega_{\max})^2$ of the window for conditioning. When the
window holds data of
fewer than two modes (an SDOF-like or flat transfer function) the 5 × 5 system has rank 4: every
member of the solution family interpolates, and the minimum-norm least-squares (SVD, relative
singular-value cut 1e-10) solution is used. The fit reproduces an exact 2-DOF hysteretic transfer
function to about 1e-14 (R1 V5; VP-28 checks 1e-10).

**Guard against spurious poles** (D-MOT-01 as amended by D-W1-10). When the fitted denominator nearly
vanishes on the window span ($\min\lvert D\rvert < 10^{-3}\max\lvert D\rvert$,
$D(x) = x^2 + C_4 x + C_5$) the window *may* hold a spurious near-real pole. The pole is accepted
when it is cancelled by a numerator zero, when it is a physical damped mode (for $\exp(+i\omega t)$ a
hysteretic mode has its pole at $x_p = (\omega_0/\omega_{\max})^2\,c(\beta)$, i.e.
$0 < \arg x_p < \pi$ with implied damping $\sin\left(\tfrac{1}{2}\arg x_p\right) \ge 10^{-4}$), or
when it is remote from the span; otherwise the complex cubic through the four nearest window points
replaces the rational form in that window. (The literal ratio test alone would destroy exact fits on
lightly damped poles.)

### 10.2 Window schemes (options 0-6)

For the interval $[f_j, f_{j+1}]$ the candidate windows are those of five consecutive points that
contain it: start index $m$ with $\max(1, j-3) \le m \le \min(j, N-4)$ (1-based).

| Option | Windows | Combination |
|---|---|---|
| 0 (SASSI2000) | all candidate windows | weighted mean, $w_m = \max\left(10^{-3},\ 1 - \lvert f - c_m\rvert/h_m\right)$ ($c_m$, $h_m$ centre and half-width of window $m$) |
| 1 (SASSI 1982) | starts 1, 5, 9, ... sharing end points; the last window is the last five points | single window |
| 2 | all candidate windows | arithmetic mean |
| 3 | the three candidates with the smallest $\lvert f - c_m\rvert$ | arithmetic mean |
| 4 | starts 2, 6, 10, ...; the first interval uses $m = 1$ | single window |
| 5 | starts 3, 7, 11, ...; the first two intervals use $m = 1$ | single window |
| 6 | independent not-a-knot cubic splines of $\operatorname{Re} H$ and $\operatorname{Im} H$ through the computed points and the $f = 0$ anchor | — |

With fewer than five SSI frequencies options 0-5 become option 6, with fewer than three complex linear
interpolation. Every scheme returns the computed values exactly at the SSI frequencies (TFI = TFU,
enforced bit-exactly). The weighting of option 0 and the shifts of options 4 and 5 are
reconstructions (R1 §10 item 4).

**Outside the computed range** (D-MOT-03): above the last SSI frequency $f_N$ the transfer function
is zero (the cut-off is a low-pass filter); below the first one a seismic transfer function is linear
between the rigid-body anchor $H(0)$ and $H_1$ ($H(0)$ = projection of the control direction on the
output DOF: 1 parallel to the input, $\cos a$ / $\sin a$ for rotated horizontals, 0 otherwise and for
rotations), a vibration transfer function is held at $H_1$.

### 10.3 Smoothing and phase adjustment

**Smoothing** (MOTION `<smo>` = $S$, options 0-5): in $[f_j, f_{j+1}]$, with the complex linear
reference $L(f)$ and the band $[A_{\min}, A_{\max}]$ of the two neighbouring amplitudes,

```math
\begin{aligned}
r &= \max\left(0,\ \frac{\lvert H\rvert}{A_{\max}} - 1\right) + \max\left(0,\ 1 - \frac{\lvert H\rvert}{\max(A_{\min},\ 10^{-12} A_{\max})}\right)\\
H_s &= L + \frac{H - L}{1 + S\,r}
\end{aligned}
\tag{10.3}
```

Values inside the band are unchanged; $S = 0$ is the identity. Genuine resonance peaks between
computed points are clipped too, hence $S = 0$ for coherent analyses.

**Phase adjustment** (MOTION `<pzadj>` = 1; decision D-MOT-05 as amended by D-W1-11):

```math
\begin{aligned}
\phi &= \operatorname{unwrap}(\arg H)\\
H' &= \lvert H\rvert\exp\left(i\left(\phi_0 + \rho\,(\phi - \phi_0)\right)\right), \qquad
\rho = \frac{1}{1 + S}\ \ \text{(options 0-5)}, \qquad \rho = 0\ \ \text{(option 6)}
\end{aligned}
\tag{10.4}
```

with the phase unwrapped along the grid from $f = 0$ and $\phi_0$ the phase at $f = 0$ (or the first
non-zero value). Referring the phase to $\phi_0$ keeps negative rigid-body anchors (a reversed input
direction reverses the response). With $\rho = 0$ the transfer function has zero differential phase,
the "upper bound" approach of incoherent analyses; VP-33 checks both cases exactly.

**SRSS transfer function** (MOTIONX `<srss>` with `SRSSTF.txt`, D-MOT-06): the modal FILE8 transfer
functions are interpolated and combined as $\lvert H\rvert = \sqrt{\sum\lvert H_k\rvert^2}$ with zero
phase, or with the phase of the coherent solution (VP-32).

---

## 11. Convolution, response spectra, relative displacements and stresses

**Code:** `sassi/core/signal.py`, `sassi/core/spectra.py`, `sassi/modules/motion.py`,
`sassi/modules/reldisp.py`, `sassi/core/stress_lib.py`, `sassi/modules/stress.py`.
**Verified by:** VP-29 (identity convolution), VP-30 (response spectra), VP-31 (baseline correction),
VP-34 (relative displacement), VP-35, VP-S1, VP-S2 (stress recovery), VP-40.

### 11.1 Convolution

The control acceleration $a(t)$ (in g, records rec1 ... rec2, scaled by `<mult>` or to `<max>`) is
zero-padded to NFFT points and transformed, $A = \operatorname{rfft}(a)$. The frequency step must
equal the FILE8 step (relative 1e-6). Then

```math
\begin{aligned}
&\text{seismic acceleration:} && a_k(t) = \operatorname{irfft}\left(H_k(f)\,A(f)\right)\\
&\text{quantities linear in displacement:} && U_g(f) = -g\,A(f)/\omega^2, \qquad r(t) = \operatorname{irfft}\left(H_r\,U_g\right)\\
&\text{vibration:} && u = \operatorname{irfft}(H F), \quad v = \operatorname{irfft}(i\omega H F), \quad a = \operatorname{irfft}(-\omega^2 H F)
\end{aligned}
\tag{11.1}
```

with $a_k(t)$ in g, the relative displacements, strains and stresses as the quantities linear in
displacement (the $f = 0$ term of $U_g$ set to 0), and $F$ the FFT of the reference load history.
Because FILE8 holds the dimensionless ratio $U/U_\text{cp}$, the same transfer function serves for
accelerations and displacements (D-CNV-06). The signal is periodic with period
$\mathrm{NFFT}\cdot\Delta t$; the zero-padded tail (quiet zone) must be long enough for the free
vibration to decay. The `.ACC` history is written for
$\min(\mathrm{NFFT},\ 1.2\cdot\mathrm{dur}/\Delta t)$ samples. VP-29 checks that $H \equiv 1$
returns the input to 1e-12.

### 11.2 Baseline correction (MOTION `<bl>` = 1)

Hudson-Housner time-domain correction (D-MOT-08): integrate the acceleration (trapezoidal rule) to
the velocity $v$; fit $v \approx c_1 t + c_2 t^2/2$ by least squares (the velocity produced by an
acceleration baseline error $c_1 + c_2 t$); subtract $c_1 + c_2 t$ from the acceleration and
integrate again. VP-31 removes a constant acceleration offset to 1e-6.

### 11.3 Response spectra (Nigam-Jennings)

For an oscillator of frequency $\omega$ and damping $\zeta$,
$\ddot x + 2\zeta\omega\dot x + \omega^2 x = -a(t)$, with the ground acceleration linear between
samples, the exact recurrence of Nigam and Jennings (1969) advances the state $s = (x, \dot x)$:

```math
s_{k+1} = A\,s_k + B_0\,a_k + B_1\,a_{k+1} \tag{11.2}
```

with $A$, $B_0$, $B_1$ functions of $\omega$, $\zeta$, $\Delta t$ (R2 G.1;
`spectra._nj_coefficients`). The spectral acceleration is the maximum absolute *absolute*
acceleration
$\lvert 2\zeta\omega\dot x + \omega^2 x\rvert$ (written to `.RS`); SV and SD are the maxima of the
relative velocity and displacement, $\mathrm{PSA} = \omega^2\,\mathrm{SD}$. Two refinements keep the
peak sampling error small: for oscillator frequencies above $0.1/\Delta t$ the record is
FFT-upsampled four times (D-MOT-07), and the piecewise-linear input is sub-divided so that every
oscillator period has at least 64 integration steps (D-W1-07; peak error
$\le 1 - \cos(\pi/64) = 0.12\,\%$). The spectra are computed over the full NFFT record at `<fstep>`
log-spaced frequencies from `<freq1>` to `<freq2>`. VP-30 checks the closed forms: step input
$\mathrm{PSA}/a_0 = 1 + \exp\left(-\zeta\pi/\sqrt{1-\zeta^2}\right) = 1.854468$ at $\zeta = 5\,\%$,
harmonic input at resonance $A/(2\zeta)$, impulse $\mathrm{PSV}/\Delta V = 0.92669$.

### 11.4 Relative displacements (RELDISP)

From the complex `.TFI` of a node and of a reference (another node's `.TFI`, or the free field: unit
amplitude and zero phase in the input direction, D-RDP-04):

```math
D(f) = \left(H_\text{node}(f) - H_\text{ref}(f)\right)\cdot U_g(f), \qquad d(t) = \operatorname{irfft}(D) \tag{11.3}
```

This avoids the drift of double integration. Above the cut-off the node TFs are zero while the
free-field reference is still 1, so a free-field-referenced relative displacement contains the small
ground displacement above $f_N$: exactly $\operatorname{irfft}\left((H - 1)\,U_g\right)$ (VP-34). The
`.TFD` file is the complex relative-displacement transfer function per unit control acceleration
(length per g).

### 11.5 Stresses and forces (STRESS)

HOUSE stores per output element a complex **recovery operator** $S$, so that a component is
$S_c\cdot u_e$ (SOLID/PLANE centroid stresses $D^* B u$ in global axes, SHELL membrane stresses and
moments per unit length $D^*_m B_m u$, $D^*_b\,\kappa$ in local axes, BEAMS local end forces
$k_L T u$ exerted on the element, SPRING $k^*(u_J - u_I)$). STRESS then

1. forms the **stress transfer functions** at the SSI frequencies,
   $\mathrm{STF}_c(f_j) = S_c\left(U_e(f_j) - r_e\right)$, where $r_e$ is the rigid-body part (the
   control motion projected on the element DOFs; $S$ annihilates it in exact arithmetic, and
   subtracting it keeps the low-frequency values accurate);
2. **interpolates the STF** (not the nodal TFs) with the MOTION schemes; below the first SSI frequency
   a seismic STF goes linearly to 0 at $f = 0$ (a rigid-body motion has no stress);
3. **convolves** with the control displacement spectrum (Eq. 11.1);
4. computes the non-linear derived quantities in the time domain, e.g. the octahedral shear stress

   ```math
   \tau_\text{oct} = \tfrac{1}{3}\sqrt{(\sigma_{xx}-\sigma_{yy})^2 + (\sigma_{yy}-\sigma_{zz})^2 + (\sigma_{zz}-\sigma_{xx})^2 + 6(\tau_{xy}^2 + \tau_{yz}^2 + \tau_{xz}^2)}
   ```

VP-35 checks that a rigid-body motion gives zero stress, $\tau_\text{oct} = (\sqrt{6}/3)\,\tau$ for
pure shear and the equilibrium of beam end forces; VP-S1 the base moment of a cantilever stick
against an independent calculation; VP-S2 the shear stress of a soil column under vertical SV against
$G^* \times \text{free-field strain}$.

---

## 12. The SHAKE equivalent-linear method (SOIL)

**Code:** `sassi/core/shake.py`, `sassi/modules/soil.py`. **Verified by:** VP-02a (closed-form layer
transfer functions, 1e-8), VP-04 (SHAKE91 sample problem), example 4.

### 12.1 Wave solution and recursion

In sublayer $m$ (local depth $z_m$ measured down from its top) a vertically propagating shear wave is
the sum of an up-going wave $E_m$ and a down-going wave $F_m$ (Schnabel, Lysmer and Seed 1972;
SHAKE91):

```math
\begin{aligned}
u_m &= E_m\exp\left(i(\omega t + k^*_m z_m)\right) + F_m\exp\left(i(\omega t - k^*_m z_m)\right)\\
k^*_m &= \omega/V^*_m, \qquad V^*_m = \sqrt{G^*_m/\rho_m}
\end{aligned}
\tag{12.1}
```

Continuity of displacement and shear stress $\tau = G^*\,\partial u/\partial z$ at each interface
gives

```math
\begin{aligned}
E_{m+1} &= \tfrac{1}{2}E_m(1 + \alpha_m)\exp(ik^*_m h_m) + \tfrac{1}{2}F_m(1 - \alpha_m)\exp(-ik^*_m h_m)\\
F_{m+1} &= \tfrac{1}{2}E_m(1 - \alpha_m)\exp(ik^*_m h_m) + \tfrac{1}{2}F_m(1 + \alpha_m)\exp(-ik^*_m h_m)\\
\alpha_m &= \frac{\rho_m V^*_m}{\rho_{m+1}V^*_{m+1}}
\end{aligned}
\tag{12.2}
```

starting from the free-surface condition $E_1 = F_1\ (= 1)$. The **within** motion at the top of
sublayer $m$ is $E_m + F_m$, the **outcrop** motion $2E_m$. Transfer functions are ratios of these
quantities; the input is applied at the top of the control sublayer as outcrop or within motion, and
the last SPRO entry is the half-space (D-SOL-02). VP-02a checks the closed form
$1/\left[\cos k^*H + i\alpha^*\sin k^*H\right]$ of a uniform layer to 1e-8 for both damping forms.

### 12.2 Strains and the equivalent-linear iteration

The shear strain at mid-height of each sublayer, from the acceleration-wave amplitudes, is

```math
\gamma_m(\omega) = \frac{E_m\exp(ik^*h/2) - F_m\exp(-ik^*h/2)}{i\omega V^*_m} \tag{12.3}
```

transformed to the time domain to obtain $\gamma_{\max} = \max\lvert\gamma(t)\rvert$. The **effective
strain** is $\gamma_\text{eff} = R_\gamma\,\gamma_{\max}$ (SOIL `<ratio>`, typically 0.5-0.7;
SHAKE91 suggests $(M - 1)/10$ for magnitude $M$). New properties are read from the strain-dependent
curves:

```math
G = G_{\max}\cdot(G/G_{\max})(\gamma_\text{eff}), \qquad \beta = D(\gamma_\text{eff}) \tag{12.4}
```

interpolated **linearly in $\log_{10}(\gamma)$** with constant extrapolation (SHAKE91's rule;
D-SOL-03). SOIL runs exactly `<iter>` iterations (SHAKE91 has no convergence test; 8 are recommended)
and lists the changes of $G$ and $\beta$ per iteration; the half-space is not iterated. For vertical
input (SOILX
`<indir>` = 1) the P-wave velocity and damping are used without iterations. Fourier components above
`EDUOPT,SOILCUTOFF` are removed (D-SOL-08; 25 Hz reproduces the SHAKE91 sample problem).

### 12.3 Hand-off to SSI

FILE88 holds per sublayer the thickness, the effective strain, $G$, $V_s$, $\beta_s$, $V_p$ and
$\beta_p$. With `SITEX,1` SITE replaces the TOPL layer properties by these, position by position (the
layer counts must match, EDU-07). For $V_p$ and $\beta_p$ the default policy keeps Poisson's ratio
(`EDUOPT,VPPOLICY,NU`; `VP` keeps $V_p$ for saturated soils) and sets $\beta_p = \beta_s$ (D-SOL-06).
The same complex-modulus form must be used in SOIL
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

HOUSE assembles four frequency-independent sparse matrices on one global DOF map: $K^*_s$ (complex)
and $M_s$ of the structure (including near-field soil) and $K^*_e$, $M_e$ of the excavated soil.
Because damping is hysteretic, none depends on the frequency; ANALYS combines them at each frequency.

### 13.1 The element library

| Element | Stiffness | Mass |
|---|---|---|
| SOLID | 8-node isoparametric hexahedron, $D^* = \lambda^* D_\lambda + G^* D_G$; Gauss rule 2×2×2 / 3×3×3 / 4×4×4 for EINT 0/1/2; for structural solids 9 Wilson incompatible modes ($1-\xi^2$, $1-\eta^2$, $1-\zeta^2$ per direction) with the Taylor centroid-Jacobian correction, condensed statically $K = K_{uu} - K_{ua}K_{aa}^{-1}K_{au}$ (passes the patch test); prisms and pyramids by repeated nodes; positive Jacobian required (EDU-05) | ½ lumped (row sum) + ½ consistent |
| PLANE | 4-node plane strain in the X-Z plane, unit thickness, 2×2 Gauss, 4 incompatible modes for structural elements; either node orientation accepted | ½ + ½ |
| BEAMS | 2-node 3D Timoshenko frame, $\phi = 12EI/(G A_s L^2)$ ($A_s = 0$: Euler-Bernoulli); bending in the 1-2 plane with $I_3$ and $A_{s2}$, in the 1-3 plane with $I_2$ and $A_{s3}$; end releases by static condensation of the local matrices; exact for end loads | consistent (translations from the Timoshenko interpolation, axial $\rho A L/6\,[2\ 1;\ 1\ 2]$, torsion $\rho(I_2+I_3)L/6\,[2\ 1;\ 1\ 2]$, no rotary inertia) |
| SHELL | flat facet: membrane = plane-stress Q4 with 4 incompatible modes (CST for triangles); plate = DKQ / DKT discrete-Kirchhoff elements (Batoz et al. 1980, 1982); **zero drilling stiffness**; warped quads projected on the mean plane | lumped $\rho t A/n$ on translations, no rotary inertia |
| TSHELL | flat facet, Mindlin-Reissner plate (section 13.2): membrane as SHELL; transverse shear by the MITC4 assumed natural strains (Bathe and Dvorkin 1985) for both EINT values; bending EINT 1 = 2×2 Gauss, EINT 0 = one point with hourglass stabilisation (coefficient 0.1); triangles: CST membrane, MITC3 shear with a condensed rotation bubble and a stabilised curl part; shear factor 5/6; drilling penalty 1e-4 × the smallest membrane diagonal × area | lumped $\rho t A/n$ on translations and the Mindlin rotary inertia $\rho t^3 A/(12n)$ on the two bending rotations; none on the drilling rotation |
| SPRING | six uncoupled global springs $k^*_d = k_d\,c(\text{damp})$ between I and J; zero length allowed | none |
| GENERAL | user 12×12 $K_R + iK_I$ and mass (÷g in weight units); with a K node, local axes as BEAMS and $T = \operatorname{blockdiag}(\Lambda, \Lambda, \Lambda, \Lambda)$ | user |
| nodal masses | MT/MR on the diagonal (÷g with MUNITS = 1), ignored on fixed DOFs | |

**Beam stiffness in one bending plane** (deflection $v$, rotation $\theta = dv/dx$; Przemieniecki
1968):

```math
k = \frac{E I}{(1 + \phi)\,L^3}\begin{bmatrix}
12 & 6L & -12 & 6L\\
6L & (4+\phi)L^2 & -6L & (2-\phi)L^2\\
-12 & -6L & 12 & -6L\\
6L & (2-\phi)L^2 & -6L & (4+\phi)L^2
\end{bmatrix} \tag{13.1}
```

with the sign of the coupling terms flipped in the 1-3 plane ($\theta_2 = -du_3/dx$). Local axes:
$e_1$ = I → J, $e_2$ toward K in the plane IJK, $e_3 = e_1 \times e_2$; BEAMS end forces are the
forces exerted on the element.

**Excavated soil.** Excavated SOLID/PLANE elements take the L-table layer of their MSET index
(D-ELM-12): $G = \rho V_s^2$, $M = \rho V_p^2$, $\rho = \text{weight}/g$, with the complex moduli of
Eq. (2.3); never incompatible modes. ETYPE 0 elements are excavated when all nodes are at or below the
ground elevation and the centroid is strictly below it (D-HOU-01). VP-H2 checks the excavated mass
$\sum\rho V$ and the stiffness identity with the layer properties.

**DOF management** (D-ELM-09). The active DOFs of a node are the union of the DOFs of its elements
minus the fixed ones; DOFs that no element defines (rotations of solid-only nodes) are removed
automatically, DOFs that are defined but not stiffened (shell drilling rotations) are kept and need
FIXROT/FIXSHLROT. Interaction DOFs are the translations of the interaction nodes.

**What the element VPs show.** VP-37: Euler-Bernoulli and Timoshenko cantilevers, end releases,
simply supported plates and bars, patch tests (exact constant stress). VP-38: NAFEMS free-vibration
benchmarks FV12, FV16 (SHELL bending; the lumped translational mass converges from below at
$O(h^2)$, hence the 32 × 32 mesh of lead decision D-W1-05) and FV32 (membrane). VP-03: the dispersion
of discrete soil columns with consistent, lumped and mixed mass (the mixed mass gives
$\omega_h/\omega = 0.99451$ at five elements per wavelength). The TSHELL problems are described in
section 13.2.

### 13.2 The TSHELL thick shell (Mindlin-Reissner)

The thin SHELL element (Kirchhoff) neglects transverse shear deformation. For thick basemats, the
shear walls of nuclear buildings or the higher modes of any plate the Mindlin-Reissner theory is the
better model: the normal to the mid-surface stays straight but not normal, so the rotations $\beta$
of the normal are independent of the slope of $w$ and the transverse shear strain carries the shear
forces. In the local axes $x'y'z'$ of the facet (defined as for SHELL), with right-hand rotations
$\theta'$:

```math
\begin{gathered}
u(z) = u' + z\,\beta_x, \qquad v(z) = v' + z\,\beta_y, \qquad \beta_x = \theta'_y, \qquad \beta_y = -\theta'_x\\[8pt]
\begin{aligned}
&\text{membrane} && \varepsilon = \left[u_{,x}\quad v_{,y}\quad u_{,y} + v_{,x}\right], && N = D_m\,\varepsilon\\
&\text{bending} && \kappa = \left[\beta_{x,x}\quad \beta_{y,y}\quad \beta_{x,y} + \beta_{y,x}\right], && M = D_b\,\kappa\\
&\text{shear} && \gamma = \left[w_{,x} + \beta_x\quad w_{,y} + \beta_y\right], && Q = D_s\,\gamma
\end{aligned}\\[8pt]
D_m = \frac{E t}{1-\nu^2}\begin{bmatrix} 1 & \nu & 0\\ \nu & 1 & 0\\ 0 & 0 & (1-\nu)/2 \end{bmatrix}, \qquad
D_b = \frac{E t^3}{12(1-\nu^2)}\begin{bmatrix} 1 & \nu & 0\\ \nu & 1 & 0\\ 0 & 0 & (1-\nu)/2 \end{bmatrix}\\[8pt]
D_s = \kappa_s G t\,I, \qquad \kappa_s = 5/6
\end{gathered}
\tag{13.2}
```

so $M_{xx}$ is positive with tension on the $+z'$ face and the plate equilibrium reads
$Q_x = M_{xx,x} + M_{xy,y}$. In the Kirchhoff limit $\gamma \to 0$ gives $\beta = -\nabla w$, the
SHELL convention.

**Shear locking and the MITC4 field.** With bilinear $w$ and $\beta$ a thin plate cannot make
$\gamma = 0$ everywhere and the element becomes far too stiff (with plain 2×2 integration of the
shear, a cantilever plate 1000 times longer than thick, 8 × 2 elements, deflects about 6500 times too
little; VP-TS1). The quadrilateral therefore uses the assumed natural strains of Bathe and Dvorkin
(1985): the covariant shear strains $\gamma_\xi = w_{,\xi} + \beta\cdot x_{,\xi}$ and
$\gamma_\eta = w_{,\eta} + \beta\cdot x_{,\eta}$ are sampled at the edge mid-points A (0, 1),
C (0, -1), B (1, 0), D (-1, 0) and interpolated,

```math
\begin{aligned}
\tilde\gamma_\xi &= \tfrac{1}{2}(1 + \eta)\,\gamma_\xi(A) + \tfrac{1}{2}(1 - \eta)\,\gamma_\xi(C), \qquad
\tilde\gamma_\eta = \tfrac{1}{2}(1 + \xi)\,\gamma_\eta(B) + \tfrac{1}{2}(1 - \xi)\,\gamma_\eta(D)\\
\gamma &= J^{-1}\left[\tilde\gamma_\xi,\ \tilde\gamma_\eta\right]
\end{aligned}
\tag{13.3}
```

integrated with 2×2 Gauss points for both EINT values. At the centre $\tilde\gamma$ is the one-point
(tying average) value of the manual's "1-point transverse shear", and the linear terms act as an
assumed-strain stabilisation of it: a literal one-point shear would leave the $w = \pm 1$ hourglass
mode without stiffness. The MITC4 field neither locks nor has a spurious mode.

**Bending.** EINT 1 (*selective*): 2×2 Gauss integration. EINT 0 (*reduced*, the manual's default):
one point at the centre plus a physical hourglass stabilisation,

```math
K_b = A\,B_0^{\mathsf{T}} D_b B_0 + \varepsilon_\text{hg}\sum_g w_g\,\lvert J_g\rvert\,(B_g - B_0)^{\mathsf{T}} D_b\,(B_g - B_0), \qquad \varepsilon_\text{hg} = 0.1 \tag{13.4}
```

$B_g - B_0$ vanishes for every linear rotation field, so constant curvature is exact on any
quadrilateral (patch test) and only the two rotation hourglass modes get $\varepsilon_\text{hg}$ times
their 2×2 stiffness. **Membrane** (both EINT values): the SHELL membrane, a Q4 with the four
Wilson-Taylor incompatible modes condensed.

**Triangles.** CST membrane. Plate: $w$ linear, rotations linear plus a cubic bubble
$f_4 = 27\,r s\,(1 - r - s)$ whose two DOFs are condensed inside the element; transverse shear from
the MITC3 field of the corner DOFs (Lee and Bathe 2004),
$\tilde\gamma = a + c\,(-(y - y_c),\ x - x_c)$, whose tangential component on each side equals the
side average of the displacement-based strain; the bubble enters through its element mean,
$\gamma_0 = a + (27/60)\,\beta_b$. The shear energy is
$A\,\gamma_0^{\mathsf{T}} D_s\,\gamma_0 + w_\text{lin}\int (c\,\mathrm{rot})^{\mathsf{T}} D_s\,(c\,\mathrm{rot})\,dA$
with the Lyly-Stenberg-Vihinen (1993) weight $w_\text{lin} = t^2/(t^2 + \alpha h^2)$, $\alpha = 0.2$
and $h$ the longest side, on the linear ("curl") part only. The plain MITC3 triangle locks: a simply
supported plate with $t/a = 1/1000$ on a 16 × 16 mesh deflects 46 % of Navier's value.

**Drilling rotation.** "A small rotational stiffness is automatically added in HOUSE" (manual): the
nodal drilling rotations are tied to the in-plane rotation of the membrane at the centre,
$\omega = (v_{,x} - u_{,y})/2$, by the penalty $k_d\sum_i(\theta'_{z,i} - \omega)^2$ with
$k_d = 10^{-4}\min(\operatorname{diag} K_m)\,A$. It vanishes for a rigid in-plane rotation, so the
free element keeps its six rigid-body modes, and it is too small to alter the membrane response.
TSHELL models therefore need no FIXROT.

**Mass, damping, recovery.** Lumped $\rho t A/n$ on each translation and $\rho t^3 A/(12n)$ on the
two bending rotations (required to reproduce the Mindlin frequencies of NAFEMS Test 21, which include
rotary inertia); no drilling inertia (an inertia on a penalty DOF would create spurious low-frequency
modes). $E^* = E\,c(\beta_s)$ with $\nu$ real, so $K^* = c(\beta_s)\,K_0$. STRESS recovers at the
centre the membrane forces NXX NYY NXY, the transverse shear forces QXZ QYZ (MITC field) and the
moments MXX MYY MXY per unit length. With `THSHLSTR,1` it adds the face stresses of D-TSH-01 from the
maxima of the eight components, which occur at different times with unknown signs:

```math
\begin{aligned}
\mathrm{SXX} &= s_N\,\frac{\mathrm{NXX}}{t} + s_M\,\frac{6\,\mathrm{MXX}}{t^2}, \qquad
\mathrm{SYY} = s_N\,\frac{\mathrm{NYY}}{t} + s_M\,\frac{6\,\mathrm{MYY}}{t^2}\\
(s_N, s_M) &= {+}{+},\ {-}{-},\ {+}{-},\ {-}{+}\\
\mathrm{TXY} &= \frac{\mathrm{NXY}}{t} + \frac{6\,\mathrm{MXY}}{t^2}\ \ ({+}{+}\ \text{only})\\
\mathrm{TXZ} &= \frac{1.5\,\mathrm{QXZ}}{t}, \qquad \mathrm{TYZ} = \frac{1.5\,\mathrm{QYZ}}{t}
\end{aligned}
\tag{13.5}
```

with TXZ and TYZ the mid-surface maxima, and the principal stresses and the plane-stress strains of
every permutation.

**Verification.** VP-TS1: membrane, bending and transverse-shear patch tests exact (1e-10) on
distorted quadrilaterals and triangles; no locking (cantilever plate $t/L = 1/1000$ within 1 % of
Kirchhoff on 8 × 2 elements and 0.02 % on 32 × 4; simply supported plates $t/a$ = 1/1000 and 1/10
with both triangulation patterns within 2 % (16 × 16) and 0.5 % (32 × 32) of Navier, observed order
2); exactly six zero-energy modes; $K^* = c(\beta)\,K_0$. VP-38T: NAFEMS Test 21 (simply supported
thick plate 10 × 10 × 1 m, 45.897, 109.44 (×2), 167.89 Hz ...) within 2 % on 16 × 16 and 32 × 32
meshes with an observed order of 2; the NAFEMS 8 × 8 rows are informative (lead decision D-W3-01: the
lumped mass converges from below at $O(h^2)$, -5.4 % on mode 4 with 8 × 8); the thin limit FV12/FV16
within 0.5 % of the Kirchhoff SHELL on the same mesh. VP-54T: the THSHLSTR face stresses $N/t$ and
$\pm 6M/t^2$ of pure membrane and pure bending states.

### 13.3 The node-numbering optimizer (HOUSEX `<optimize>` = 1)

A profile (skyline) solver, the solver of the original SASSI, stores every row of the matrix from its
first non-zero to the diagonal, so its cost grows with the profile $\sum_i(i - \mathrm{first}_i)$ and
the bandwidth. The optimizer (D-HOU-03) renumbers the nodes by reverse Cuthill-McKee (Cuthill and
McKee 1969; George and Liu 1981): the breadth-first search starts from **all interaction nodes at
once**, so they get the last numbers, ascending in their bottom-up table order (EDU-21 by
construction), and the dense $X_{ff}$ block that ANALYS adds stays in the last rows of the profile.
Components without interaction nodes are
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

Two points of the free-field surface a horizontal distance $D$ apart do not move alike. Part of the
difference is a **delay** (the wave sweeps across the foundation at an apparent velocity: *wave
passage*), the rest is random (*incoherence*, from scattering). The **unlagged (plane-wave) coherency**
$\gamma(f, D)$ is the correlation coefficient of the Fourier amplitudes of the two motions at
frequency $f$ after the delay has been removed: 1 at $D = 0$ or $f = 0$, decaying with frequency and
distance. SASSI applies the two effects separately: HOUSE builds complex **incoherency factors**
$s_i(\omega)$ per SSI frequency, motion component (X, Y, Z) and interaction node, from the coherency
matrix of the interaction nodes and from the wave-passage delays, and writes them to FILE77; ANALYS
multiplies the coherent free field of SITE by them:

```math
u_{i,\text{inc}}(\omega) = s_i(\omega)\,u_{i,\text{coh}}(\omega) \tag{14.1}
```

### 14.2 Coherency models

`WPASS <cohf>` selects the model (D-INC-01). The coefficients of the Abrahamson models are **data**,
transcribed from the primary reports into `sassi/data/coherency/*.json` with the citation (report,
equation, table, page) and points read from the published figures, which the unit tests reproduce;
they are never typed from memory (D-INC-04). A model without a verified file refuses to run with a
"coefficients not available" error that names the alternative.

| `<cohf>` | Model | $\gamma(f, D)$ | Status and source |
|---|---|---|---|
| 1 | Luco and Wong (1986) | $\exp\left[-(\gamma_c\,\omega D/V_s)^2\right]$, $\omega = 2\pi f$ in rad/s (D-INC-05); $\gamma_c$ = INCOH `<gammax>`, `<gammay>`, `<gammaz>` per component (≥ 0.1, Error 58), $V_s$ = INCOH `<alpha>`, the mean shear-wave velocity (Error 59); $D$ = plain horizontal distance | available; form confirmed from the EPRI/NRC CLASSI incoherency presentation (NRC ADAMS ML072620217) |
| 2 | Abrahamson (1993), all soil types | - | refused: the seminar paper could not be retrieved |
| 3 | Abrahamson (2005), all sites, surface foundations | Eq. (14.2) with $a_2 f_c$ in the second factor | available: EPRI 1012968, Eq. 3-1, Tables 3-4 and 3-5 (also EPRI 1015110 Eq. 5-2) |
| 4 | Abrahamson (2006), all sites, embedded foundations | - | refused: EPRI 1014101 could not be retrieved (EPRI 1015110 recommends the hard-rock model for embedded foundations) |
| 5 | Abrahamson (2007), hard rock | Eq. (14.2) | available: EPRI 1015110 chapter 6, Eq. 6-1, Tables 6-1 and 6-2 |
| 6 | Abrahamson (2007), soil sites, surface foundations | Eq. (14.2) | available: EPRI 1015110 chapter 7, Eq. 7-1, Tables 7-1 and 7-2 |
| 7 | user tables | `COHXUSER`, `COHYUSER`, `COHZUSER` on the grids `FREQCOH` × `DISTCOH`: bilinear in $(f, D)$, clamped at the table ends, $\gamma(f, 0) = 1$ (D-INC-11) | available |

The Abrahamson models share the plane-wave form ($f$ in Hz, separation $\xi$ in metres; feet are
converted when the gravity is > 20):

```math
\gamma_\text{pw}(f, \xi) = \left[1 + \left(\frac{f\tanh(a_3\xi)}{a_1 f_c(\xi)}\right)^{n_1}\right]^{-1/2}\cdot\left[1 + \left(\frac{f\tanh(a_3\xi)}{a_2\,[f_c(\xi)]}\right)^{n_2}\right]^{-1/2} \tag{14.2}
```

with separate horizontal (X, Y) and vertical (Z) coefficient sets; $a_1$, $a_2$, $a_3$, $n_1$, $n_2$
and $f_c$ are constants or functions of $\xi$ given in the files (for example
$f_c = 27.9 - 4.82\ln(\xi+1) + 1.24\left[\ln(\xi+1) - 3.6\right]^2$ for the horizontal hard-rock
model), and the factor $[f_c]$ is present in the 2005 model only.
Beyond the end of the formula range (150 m for the 2007 models, whose coefficient expressions become
non-physical further out) the coherency is held at its end value, the more coherent and conservative
choice, and HOUSE warns when separations leave the validated range 0-150 m.

Models 2-7 use the **directional distance** of the manual, with $\Delta X'$, $\Delta Y'$ the
separations along and across Line D (the horizontal line at `WPASS <ang>` from the X axis) and
$\alpha$ = INCOH `<alpha>` the directionality factor:

```math
D_{ij} = \sqrt{2\left(\alpha\,\Delta X'^2 + (1 - \alpha)\,\Delta Y'^2\right)} \tag{14.3}
```

$\alpha = 0.5$ gives the Euclidean distance; $\alpha = 0.1$ weights the $Y'$ separations three times
more than the $X'$ separations (VP-27 checks $\sqrt{2\alpha}$, $\sqrt{2(1 - \alpha)}$ and the ratio
1/3). These models are applied only with the wave-passage option on (HOUSE `<wpass>` = 1, EDU-26
otherwise); $V_\text{app} = 10^9$ suppresses the delays.

### 14.3 Incoherent spatial modes

For every SSI frequency and component HOUSE forms the coherency matrix of the $N$ interaction nodes,
$\Sigma_{ij} = \gamma(f, D_{ij})$ (horizontal projections, unit diagonal, real symmetric), and
decomposes it spectrally (LAPACK syevd):

```math
\Sigma = \Phi\Lambda\Phi^{\mathsf{T}}, \qquad \lambda_1 \ge \lambda_2 \ge \ldots \ge 0, \qquad \sum_k\lambda_k = \operatorname{trace}\Sigma = N \tag{14.4}
```

Column $\phi_k$ is the $k$-th **incoherent spatial mode** and $\lambda_k$ its variance; negative
round-off eigenvalues are clipped to 0. The listing reports the trace check
$\lvert\sum\lambda - N\rvert/N < 10^{-8}$ (manual Eq. 6.1; a warning when it fails), the contribution
$100\,\lambda_1/N$ of the first mode and the number of modes that carry 90 % of the variance
(Eq. 6.4). When the X and Y matrices are identical one decomposition serves both.

Eigenvectors are defined up to their sign. D-INC-03 orients every mode so that $\sum_i\phi_{ik} \ge 0$
(`EDUOPT,INCOHSIGN,ADJUST`, the default; `RAW` keeps the eigensolver's signs); a mode whose sum
vanishes (an antisymmetric mode of a symmetric layout) is oriented by its first significant component.
**Repeated eigenvalues** (mirror-image modes of a symmetric layout) leave the eigenvectors undetermined
within their eigenspace, and the deterministic sums below depend on that choice. With ADJUST the
eigenspace gets a canonical basis that depends on the subspace only (D-W3-06): first the projection of
the uniform field $\mathbf{1}$ (the coherent motion), then pivoted orthogonalisation in node order,
using only the projector $P = QQ^{\mathsf{T}}$. The factors are then reproducible and independent of
round-off and of the eigensolver. Coincident plan positions of nodes on different levels give equal
rows of $\Sigma$; EDU-17 warns
when projections from different levels are closer than 0.1 × the mesh size, and `EDUOPT,INCOHMERGE,1`
merges coincident positions (D-INC-09).

### 14.4 Synthesis of the factors

INCOH `<nmodes>` selects the modes (0 = all, $k > 0$ = the $k$ largest, $-k$ = mode $k$ only) and
HOUSEX `<supmode>` the superposition (0 Linear, 1 Quadratic):

```math
\begin{aligned}
&\text{stochastic simulation (SS):} && s^{(r)} = \sum_k\sqrt{\lambda_k}\,\phi_k\exp\left(i\theta_k^{(r)}\right)\\
&\text{algebraic sum (AS):} && s = \sum_k\sqrt{\lambda_k}\,\phi_k\\
&\text{single mode } k \text{ (SRSS):} && s^{(k)} = \sqrt{\lambda_k}\,\phi_k
\end{aligned}
\tag{14.5}
```

with the phases $\theta_k^{(r)}$ uniform in $[-\mathrm{RandPhz}, +\mathrm{RandPhz}]$ for SS, zero
modal phases for AS and INCOH `<nmodes>` = $-k$ for a single mode. The simulation is stochastic when
a seed is non-zero and RandPhz > 0 (D-INC-02). With RandPhz = 180° the phases cover the circle,
$E\left[\exp\left(i(\theta_k - \theta_l)\right)\right] = \delta_{kl}$, and every sample has the target
coherency on average: $E\left[s\,s^{\mathsf{H}}\right] = \Phi\Lambda\Phi^{\mathsf{T}} = \Sigma$. VP-I2
runs 50 samples (FILE77001 ... FILE77050) and checks that the mean of $s\,s^{\mathsf{H}}$ over $N_s$
samples converges to $\Sigma$ with the exact variance of the estimator, i.e. an RMS error proportional
to $1/\sqrt{N_s}$ (observed log-log slope -0.499). The random streams are reproducible (D-INC-07): X
and Y use the two children of `SeedSequence(HSeed)`, Z `SeedSequence(VSeed)`, one PCG64 generator per
sample drawing one phase per mode at every frequency in ascending order, so sample $r$ does not
depend on how many samples are computed together.

As $f \to 0$, $\Sigma \to \mathbf{1}\mathbf{1}^{\mathsf{T}}$, $\lambda_1 = N$,
$\phi_1 = \mathbf{1}/\sqrt{N}$ and the AS factor is $s = \mathbf{1}$: the coherent motion and a
zero-frequency ATF of 1.00 are recovered (the manual's check, VP-27). The single modes $k > 1$ vanish
as $f \to 0$, so MOTION interpolates their transfer functions to $H(0) = 0$.

| Approach | HOUSE input | Runs | Results |
|---|---|---|---|
| SS, "Simulation Mean" (the reference approach) | HSeed, VSeed ≠ 0, RandPhz > 0; HOUSEX `<nsim>` = $N_s \le 50$; Linear | one HOUSE run; ANALYS `<simul>` = $N_s$: FILE1X/Y/Z × FILE77sss → FILE8{3(s-1)+d} | MOTION per sample; mean of the ISRS (`AVERAGE`) |
| AS (deterministic) | zero seeds or RandPhz = 0; Linear | one run | FILE8 used directly; EDU-13: valid for rigid foundations only |
| SRSS TF (Quadratic) | `<nmodes>` = $-k$, `<supmode>` = 1 | one HOUSE + ANALYS run per mode | MOTIONX `<srss>` with `SRSSTF.txt`: $\lvert H\rvert = \sqrt{\sum_k\lvert H_k\rvert^2}$, zero phase or the coherent phase (D-MOT-06) |
| SRSS FRS (Linear, single modes) | `<nmodes>` = $-k$, `<supmode>` = 0 | one run per mode | SRSS of the modal ISRS (line mathematics `SRSS`) |

For deeply embedded foundations the per-level approach decomposes the coherency level by level and
`BUILDFILE77` assembles the per-level FILE77 files into one; its accuracy rests on the consistent mode
signs of D-INC-03.

### 14.5 Wave passage and multiple excitation

Wave passage (HOUSE `<wpass>` = 1) delays the motion along Line D, measured from the ANALYS control
point (D-INC-06):

```math
\tau_i = \frac{(x_i - x_c)\cos a + (y_i - y_c)\sin a}{V_\text{app}}, \qquad s_i \leftarrow s_i\exp(-i\omega\tau_i) \tag{14.6}
```

with $a$ = WPASS `<ang>` and $V_\text{app}$ = WPASS `<appv>` (> 0, Error 113). With
$\exp(+i\omega t)$ a later arrival is a phase lag; for coherent motion the factors are pure delays,
$\lvert s\rvert = 1$ (VP-27 checks the phase difference $-\omega L/V_\text{app}$ of two nodes $L$
apart). **Multiple excitation** (HOUSE `<me>` = 1, which needs `<wpass>` = 1) multiplies, after
incoherency and wave passage, the factors of the interaction nodes of zone $k$ (ME: nodes `nfirst`
... `nlast`) by its spectral amplification ratio $\mathrm{SAR}_k(\omega)$, one value per SSI frequency
(AMP; complex with HOUSE `<cmplxspec>` = 1; modulus in $[0, 10]$; $\mathrm{SAR} = 1$ outside the
zones; Errors 115-119; D-INC-10).

### 14.6 Application in ANALYS: free-field load (FFL) or free-field motion (FFM)

With $S = \operatorname{diag}(s)$ (the factor of the control-motion direction multiplies the three
translations of a node):

```math
\begin{aligned}
&\text{FFL:} && b_f = A_f^{\mathsf{T}}\left(s \odot X_{ff}\,U'_f\right) = A_f^{\mathsf{T}} S\,X_{ff}\,U'_f\\
&\text{FFM:} && b_f = A_f^{\mathsf{T}} X_{ff}\left(s \odot U'_f\right) = A_f^{\mathsf{T}} X_{ff}\,S\,U'_f
\end{aligned}
\tag{14.7}
```

with FFL selected by ANALYSX `<ffm>` = 0 (the default) and FFM by ANALYSX `<ffm>` = 1. FFL multiplies
the coherent free-field load of every interaction node by its factor (the form validated by EPRI);
FFM multiplies the motion before the impedance acts on it and is meant for surface foundations (a
warning with embedded interaction nodes, D-INC-12). With $s = \mathbf{1}$ both reproduce the coherent
load exactly (VP-26: difference 0).

**FFL versus FFM (lead decision D-W3-05).** The two loads differ by the commutator,
$b_\text{FFL} - b_\text{FFM} = A_f^{\mathsf{T}}[S, X_{ff}]\,U'_f$. For a rigid massless surface
foundation under vertically propagating waves the free field is a rigid-body translation,
$U'_f = T\,t_d$ ($T$ of Eq. 8.3, $t_d$ the unit vector of the input direction), and the foundation
response is $q = K_G^{-1}T^{\mathsf{T}}b$ with $K_G = T^{\mathsf{T}}X_{ff}\,T$, so

```math
K_G\,(q_\text{FFL} - q_\text{FFM}) = T^{\mathsf{T}}[S, X_{ff}]\,T\,t_d, \qquad \left(T^{\mathsf{T}}[S, X_{ff}]\,T\right)^{\mathsf{T}} = -T^{\mathsf{T}}[S, X_{ff}]\,T \tag{14.8}
```

because $S$ and $X_{ff}$ are both (complex) symmetric. An antisymmetric matrix has a zero diagonal,
hence: the generalised force in the input direction, $t_d^{\mathsf{T}}K_G\,q$, is the same for FFL and
FFM; so is the response of a foundation DOF that $K_G$ does not couple to the others (the vertical
translation of a doubly symmetric mat); and coherent factors ($S = I$) make the commutator vanish.
The horizontal translations, coupled to rocking through $K_G$, and the rotations do differ: the
manual's statement that
FFL and FFM give identical results for a surface rigid foundation on rock holds exactly only in these
cases. VP-26 (20 m mat on rock, hard-rock coherency, AS) checks the identities (≤ 3e-10) and Eq. (14.8)
evaluated independently from FILE11 and FILE77 (residual 4e-10); the foundation ATF differences, about
1e-3 (X) and 3e-2 (Y), are reported as informative. VP-I1 compares the ANALYS ATF of a rigid mat under
Luco-Wong coherency with the rigid-foundation average $K_G^{-1}T^{\mathsf{T}}b$ computed independently
from FILE3, FILE1 and FILE77, for FFL and FFM (3e-10, tolerance 1e-8): the ATF is 1 at $f \to 0$ and
decreases with frequency (0.86 at 7.3 Hz).

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
  by $\operatorname{sgn}(x_i - x_j)$ (section 8.1); the impedance is the inverse of the in-plane block,
  stored in the 3D layout with zero UY rows and columns. The free field is evaluated at $(x, y_c, z)$
  (the $y$ coordinate is ignored) and the coordinate transformation angle must be 0 or 180°. ANALYS
  `<simul>` = 1 analyses the X and Z cases (FILE8X, FILE8Z).
* MOTION, RELDISP and STRESS (PLANE: SXX SZZ TXZ) work on FILE8 and FILE4 unchanged.

Not available in 2D: the global impedance (3D only) and incoherency, wave passage and multiple
excitation (manual §6.5.4). The manual does not recommend 2D SSI for design: a plane-strain model
radiates energy along the whole length of the structure and overestimates radiation damping, which is
unconservative.

**Verification.** VP-T1: the POINT2 far field against Kausel's exact line-load Green functions of the
same layered column (R1 V7: 0.33 % at $\lvert x\rvert = 2R_0$, 0.3 % beyond; lead decision D-W3-12).
VP-T3: the 2D zero-SSI identity, an embedded PLANE excavation whose structure equals the soil, gives
$U = U'_f$ at every node for vertical SV and P and an inclined SV wave (5e-15, tolerance 1e-8).
VP-42: a rigid strip of half-width $B$ on a layer of thickness $H$ over a rigid base, loaded through
FORCE (three load cases), against Jakub and Roesset (1977), $\nu = 0.30$:

```math
\frac{K_x}{G} = 1.175\left(1 + 2.15\,\frac{B}{H}\right), \qquad \frac{K_\phi}{G B^2} = 2.394\left(1 + 0.17\,\frac{B}{H}\right) \qquad \left(\frac{1}{8} \le \frac{B}{H} \le \frac{1}{2}\right) \tag{15.1}
```

The FE stiffness converges as $O(h)$ (singular edge stresses of a rigid punch), so the criterion uses
the Richardson value $2K(32) - K(16)$ (observed mesh-difference ratio about 2): within 5 % (observed
0.3 % to 4.5 %); the chain FORCE + ANALYS equals the direct rigid-strip impedance
$T^{\mathsf{T}}X_{ff}\,T$ of FILE3, and $K_x$ decreases as $H$ grows (it vanishes for a half-space in
plane strain).

### 15.2 Symmetry and antisymmetry planes (SYMM)

A structure symmetric about a vertical plane, on a horizontally layered site, under a loading that is
symmetric or antisymmetric about that plane, has a response with the same property, so a half model
(one plane) or a quarter model (two orthogonal planes, parallel to XZ and YZ; in 2D a line parallel to Z)
can be analysed. With $n$ the plane normal, the reflection $P = I - 2nn^{\mathsf{T}}$ and the sign
$s = +1$ (symmetric, SYMM type 0) or $-1$ (antisymmetric, type 1), the field obeys
$u(x') = s\,P\,u(x)$ for the mirror image $x'$ of $x$. HOUSE fixes on the plane nodes the DOFs that
this makes zero: for symmetric
loading the normal translation and the two in-plane rotations, for antisymmetric loading the two
in-plane translations and the normal rotation (rotations are pseudo-vectors, D-ANL-12). A horizontal X
input is antisymmetric about a plane normal to X and symmetric about a plane normal to Y; a vertical
input is symmetric about both.

**The soil: image superposition (D-W3-13).** ANALYS needs the flexibility of the *reduced* interaction
set. The interaction forces of the full model have the symmetry of the loading,
$f(j_S) = s_S P_S\,f(j)$, so the displacement at a reduced-model node $i$ is the sum over the $2^k$
combinations $S$ of the $k$ planes (the empty combination is the node itself):

```math
\begin{aligned}
F_\text{red}(i, j) &= \sum_S s_S\,F(i, j_S)\,P_S\\
&= F(i, j) + s\,F(i, j')\,P \qquad \text{(one plane)}
\end{aligned}
\tag{15.2}
```

with $j_S$ the image of $j$ in the planes of $S$ and $s_S$, $P_S$ the products of their signs and
reflections. (The earlier text $F(i,j) \pm P\,F(i,j')\,P$ of D-ANL-12 was wrong; for a layered site
$F(i, j')\,P = P\,F(i', j)$, so the reflection belongs on the load side.) A node on a plane carries
half (on two planes a quarter) of the full-model force, exactly as the structure, the excavated soil,
the masses and the loads of a half model carry half of that node's full-model values. For the DOFs
the symmetry constrains, $(I + sP)$ is zero: their rows and columns of $F_\text{red}$ vanish and they
are removed **before** the inversion $X_\text{red} = F_\text{red}^{-1}$; $F_\text{red}$ is
symmetrised like $F_{ff}$. ANALYS checks at every frequency that the free field has the declared
symmetry, $U'(\text{image}) = s\,P\,U'$, and stops otherwise.

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

1. **HOUSE**, iteration 0 (`.liq` ≠ 1): $G = \mathrm{GFAC}\cdot G_\text{layer}$,
   $\beta_s = \mathrm{DFAC}\cdot\beta_{s,\text{layer}}$,
   $\beta_p = \mathrm{DFAC}\cdot\beta_{p,\text{layer}}$; later (`.liq` = 1):
   $G = G_{\max}\cdot(G/G_{\max})(\gamma_\text{eff})$,
   $\beta_s = \beta_p = D(\gamma_\text{eff})$ (curve ICURVE of FILE73); one internal material per
   element (density and Poisson's ratio of the element material); writes FILE78 (the properties
   used).
2. **ANALYS**: New Structure restart (`<mode>` 1): the impedance $X_{ff}$ of the initiation run
   (`COOXqqq`) is reused.
3. **STRESS** (`<iter>` = 1, once per input direction): strain histories of the nonlinear elements,
   $\gamma_\text{eff} = \mathrm{ESF}\cdot\max_t\lvert\gamma(t)\rvert$ → FILE74 (one per direction:
   FILE74X, FILE74Y, FILE74Z).
4. **COMBXYZSTRAIN**: $\gamma_\text{eff} = \sqrt{\gamma_X^2 + \gamma_Y^2 + \gamma_Z^2}$ per element →
   FILE74.
5. Converged when

   ```math
   \max\,\lvert\Delta G/G\rvert < 2\,\% \quad\text{and}\quad \max\,\lvert\Delta\beta\rvert < 0.5\,\%\ \text{(absolute)} \tag{16.1}
   ```

   at most 8 iterations.

**What the material and the .pin factors mean** (requirements §4.4 item 7). The M-table material of a
nonlinear element is the **low-strain** near-field soil: density, Poisson's ratio and
$G_{\max} = \rho V_s^2$, the reference of the $G/G_{\max}$ curve (never strain-compatible values,
which the curve would soften a second time). $G_\text{layer}$, $\beta_\text{layer}$ are the free field
at the element: the TOPL layer that contains the element's mid-depth, with its L properties.
$\mathrm{GFAC} = \mathrm{DFAC} = 1$ ("same shear modulus as in free-field") starts the near field
from the free-field state; a looser backfill starts lower ($\mathrm{GFAC} < 1$). GFAC and DFAC choose
the starting point only: as in SHAKE, the converged strain-compatible state is set by $G_{\max}$, the
curves and the motion. When $G$ changes, Poisson's ratio is kept (the constrained modulus follows
$G$) and the P-wave damping equals the shear damping, as in SOIL (D-SOL-06). Curve interpolation is
that of SOIL: linear in $\log_{10}(\gamma)$, constant outside the tabulated range.

**Strain measures** (ISTR; D-NLS-03): 0 = the largest engineering shear-strain component
$\max(\lvert\gamma_{xy}\rvert, \lvert\gamma_{xz}\rvert, \lvert\gamma_{yz}\rvert)$
($\lvert\gamma_{xz}\rvert$ in 2D), the SHAKE measure; 1 = the octahedral shear strain

```math
\tfrac{2}{3}\sqrt{(\varepsilon_x-\varepsilon_y)^2 + (\varepsilon_y-\varepsilon_z)^2 + (\varepsilon_z-\varepsilon_x)^2 + 1.5\,(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)}
```

(the maximum shear strain $\sqrt{(\varepsilon_x-\varepsilon_z)^2 + \gamma_{xz}^2}$ in 2D). For simple
shear $\gamma$, ISTR 0 gives $\gamma$ and ISTR 1 gives $\sqrt{2/3}\,\gamma$.

**Convergence measure.** Row $k$ of `NLSOIL_CONVERGENCE.TXT` compares the properties *used* in
iteration $k$ (FILE78) with the strain-compatible properties of its response (FILE74):
$\Delta G/G = (G_\text{new} - G_\text{used})/G_\text{new}$ (relative to the new value, as in SHAKE)
and $\Delta\beta = \beta_\text{new} - \beta_\text{used}$ in percentage points.
`NLSSIITER,<variables>,<maxit>` runs the command lists until (16.1) holds; `NLSSIRESET` starts a new
analysis (`<model>.liq` = 0). The files carry the FILE4 hash of the properties they belong to, so a
stale strain file of another iteration is refused.

Each iteration is a New Structure restart (2-4 times faster than an initiation run). **VP-N1**: a
laterally uniform column of near-field soil inside an FV excavation of a uniform deposit, with the free
field and the excavated soil at the strain-compatible properties of a SOIL run of the same column,
motion and curves. Without a structure there is no secondary nonlinearity, so (a) the iterations started
from the low-strain properties converge monotonically ($\max\lvert\Delta G/G\rvert$ = 318, 42, 10,
2.6, 0.65 %) to the SHAKE profile, within 5 % on $G/G_{\max}$ and 10 % on the damping of every layer
(observed 1.3 % and 0.5 %), and (b) started from the free field ($\mathrm{GFAC} = \mathrm{DFAC} = 1$)
they are converged at once. **VP-N2**: an iteration
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
(BL, BR, TR, TL in the panel axes $e_h$ horizontal and $e_v = Z$) are used, with $u = d\cdot e_h$,
$w = d\cdot e_v$, the width $L$ and the height $H$:

```math
\begin{aligned}
\gamma &= \frac{1}{2}\,\frac{(u_{TL} - u_{BL}) + (u_{TR} - u_{BR})}{H} + \frac{1}{2}\,\frac{(w_{BR} - w_{BL}) + (w_{TR} - w_{TL})}{L}\\
\kappa &= \frac{(w_{TR} - w_{TL}) - (w_{BR} - w_{BL})}{L H}\\
\varepsilon_v &= \frac{1}{2}\,\frac{(w_{TL} - w_{BL}) + (w_{TR} - w_{BR})}{H}
\end{aligned}
\tag{17.1}
```

with $\gamma$ the shear strain (Disp. Opt 1), $\kappa$ the curvature (Disp. Opt 2) and
$\varepsilon_v$ the axial strain (informative). A rigid translation or in-plane rotation gives zero
for all three. A nonlinear spring uses its elongation $u_J - u_I$ along one DOF. The histories come
from RELDISP (`.THD`, free-field reference,
which cancels in every difference) for the X, Y and Z inputs, added time step by time step by
COMB_XYZ_THD (D-NON-08), because the nonlinear behaviour must be driven by the simultaneous
three-component input (manual §1.5.4).

**Equivalent-linear properties.** For every element (requirements §4.15):

```math
\begin{aligned}
x_\text{eq} &= \mathrm{EDF}\cdot\max_t\lvert x(t)\rvert\\
E_\text{new} &= E_\text{el}\cdot\frac{K_\text{sec}(x_\text{eq})}{K_\text{el}}, \qquad K_\text{sec} = \frac{F_\text{bb}(x_\text{eq})}{x_\text{eq}}, \qquad K_\text{el} = \frac{Y_1}{X_1}\\
\xi_h &= \frac{E_D}{4\pi E_S}, \qquad E_S = \frac{x_\text{eq}\,F(x_\text{eq})}{2}\\
\xi &= \min\left(\text{cutoff},\ \text{scale}\cdot\xi_h + [\xi_\text{el}]\right)
\end{aligned}
\tag{17.2}
```

with $\mathrm{EDF} \approx 0.8$ (0.7-0.9), $K_\text{el}$ the first BBC slope (D-NON-11), $\xi_h$
from the stabilised loop at $x_\text{eq}$, and scale 0 → 1, cutoff 0 → none (D-NON-04).
$F_\text{bb}$ is the backbone curve (BBC): odd, piecewise linear through the origin and the user
points (point 1 = the cracking point, the end of the elastic range), held constant beyond the last
point. For a panel the modulus of its own material is scaled with Poisson's ratio kept, so shear,
axial and bending stiffness degrade together (each panel needs a material of its own, D-NON-12); for a
spring the spring constant of its DOF is scaled. $E_D$ is the area of the third of three symmetric
cycles at $\pm x_\text{eq}$ through the hysteresis model. The ductility
$\mu = \max\lvert x\rvert/x_\text{cr}$ and the force-reduction factor
$F_\mu = K_\text{el}\max\lvert x\rvert/\lvert F(t^*)\rvert$ at the time $t^*$ of $\max\lvert x\rvert$
(the inelastic absorption factor of ASCE 43-05) are reported.

**Hysteresis models** (BBC type / Force Opt, D-NON-03):

* **GMR, General Masing Rule** (springs): virgin loading on the BBC; a branch from the reversal point
  $(x_r, F_r)$ is $F = F_r + 2F_\text{bb}\left((x - x_r)/2\right)$ (Masing 1926), with the extended
  memory rules (Pyke 1979; Kramer 1996): a branch that reaches the previous reversal continues on the
  branch it left, one that reaches the largest past excursion joins the backbone. Symmetric loops are
  closed with $E_D = 8\int_0^x F\,dx - 4x\,F(x)$; for an elastic-perfectly-plastic BBC
  $\xi_h = 2(x - x_y)/(\pi x)$ and $K_\text{sec} = F_y/x$ (VP-45, to round-off).
* **CMS, Cheng-Mertz Shear** (panels; low-rise RC walls): the S1 model transcribed rule by rule from
  subroutine HYST04 of INRESB-3D-SUP (Cheng and Mertz 1989): elastic until cracking, loading on the
  multi-linear backbone; unloading in three force bands with the degrading stiffnesses
  $\mathrm{S1} = \mathrm{SI}\,\min\left(1.4675\,(\mathrm{DC}/\mathrm{DMAX})^{0.345},\ 1\right)$,
  $\mathrm{S2} = \mathrm{SI}\,\min\left(0.7761\,(\mathrm{DC}/\mathrm{DMAX})^{0.5195},\ 1\right)$,
  $\mathrm{S3} = \mathrm{SI}\,\min\left(0.0707\,(\mathrm{DC}/\mathrm{DMAX})^{1.369},\ 1\right)$
  ($\mathrm{SI} = \mathrm{PC}/\mathrm{DC}$ the initial stiffness); pinched reloading below
  $0.75\,\mathrm{PC}$ with the slip stiffness
  $\mathrm{SR} = \mathrm{SI}\,\min\left((\mathrm{DC}/\mathrm{DMAX})^{1.02},\ 1\right)$; reloading
  toward $0.95\,\mathrm{PMAX}$ on the first unloading branch and then toward
  $(1.04\,\mathrm{DMAX},\ \mathrm{PMAX})$; up to ten stored small loops. The original is
  load-incremental; here every rule is a straight branch followed exactly in displacement control.
* **TAK, Takeda** (experimental, `EDUOPT,NONEXT,1`): trilinear backbone through the cracking and yield
  points, flat beyond; unloading after yield with $K_u = K_{cy}\,(D_y/D_{\max})^{0.4}$,
  $K_{cy} = (P_y + P_c)/(D_y + D_c)$ (Takeda, Sozen and Nielsen 1970; Otani 1974). The manual allows
  TAK only with the elastic damping, so $\xi = \xi_\text{el}$.
* CMB (Cheng-Mertz Bending) is not included in the manual's version and is refused.

**Iteration and convergence.** NONLINEAR writes the new properties into `<model>_new.hou`; the next SSI
analysis is a New Structure restart. The state machine of the manual: without `PANEL.NON` /
`SPRING.NON` the run is the *elastic run* and creates the file; with `.NON` = 1 it is an *iteration*
that compares the properties of the analysis just made (`*_EQL_Matl_Prop.txt`) with the new ones.
Convergence (D-NON-06): $\max\lvert E_\text{new} - E_\text{old}\rvert/E_\text{old} < 2\,\%$ (relative
to the previous value) and $\max\lvert\xi_\text{new} - \xi_\text{old}\rvert < 0.5\,\%$ (absolute), at
most 10 iterations. When the wall shear is set by the inertia forces (force control) the fixed point is
where the backbone force equals the demand $F_d$; one iteration $x \leftarrow F_d\,x/F_\text{bb}(x)$
reduces the distance to it only by the factor $1 - x\,F_\text{bb}'/F_\text{bb}$, close to 1 just after
cracking, where the backbone is nearly flat: hence the manual's advice to use smooth BBCs.

**Shear capacities and backbones** (SHEAR, BBCGEN; spec 10 §3.11). ACI 318-08, Wood (1990), Barda et
al. (1977) and Gulec and Whittaker (2009) capacities with $\sqrt{f'_c}$ in psi (the formulas are
listed in
[OPTION_NON.md §5](../user/OPTION_NON.md#5-shear-capacities-and-backbone-generation-shear-bbcgen)).
BBCGEN (D-NON-09) writes a 22-point CMS curve: point 1 = cracking
$(V_\text{cr}/(G A_W),\ V_\text{cr})$ with $V_\text{cr} = 3\sqrt{f'_c}\,A_W$ (ASCE 4-17) or
$\mathrm{CFL}\cdot V_u$; points 2-21 equally spaced in strain up to the yield point $(0.004,\ V_u)$ on
$V = V_\text{cr} + (V_u - V_\text{cr})\left[1 - (1 - \xi)^2\right]$,
$\xi = (\gamma - \gamma_\text{cr})/(\gamma_y - \gamma_\text{cr})$ (zero slope at yield); point 22 =
failure $(0.02,\ 1.02\,V_u)$. Its first slope is $G A_W$ by construction.

**Verification.** VP-45: GMR loops of an elastic-perfectly-plastic BBC at 1.25 ... 50 $x_y$ (closed,
$\xi_h$ and $K_\text{sec}$ to 1e-6, observed 2e-16) and the NONLINEAR state machine. VP-46: the SHEAR
numbers of the spec 10
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
extended Masing rules for unloading and reloading (after a reversal at $(\gamma_R, \tau_R)$:
$\tau = \tau_R + 2F\left((\gamma - \gamma_R)/2\right)$; a branch that reaches the previous reversal
continues on the curve of the previous cycle; beyond the largest past strain the backbone is
followed):

```math
\tau = \frac{G_0\,\gamma}{1 + \beta\left(\lvert\gamma\rvert/\gamma_r\right)^s}\ \left[+\ \eta\,\frac{d\gamma}{dt}\right] \tag{18.1}
```

NLSLAYER gives $\beta$, $s$, $\gamma_r$ (percent) and the viscosity $\eta$ per SOIL sublayer, or
`curvefit` = 1 fits the backbone to the sublayer's $G/G_{\max}$ curve. $\beta$ and $\gamma_r$ enter
only through $\beta/\gamma_r^s$, so the fit sets $\beta = 1$ and returns $\gamma_r$ and $s$, with $s$
restricted to $[0.2, 1]$ so that the stress never decreases (SN-3). The hysteretic damping of a
Masing loop of amplitude $\gamma_a$ is

```math
\begin{aligned}
D &= \frac{2}{\pi}\left[\frac{2\int_0^{\gamma_a} F\,d\gamma}{\gamma_a\,F(\gamma_a)} - 1\right]\\
\int_0^{\gamma_a} F\,d\gamma &= \frac{G_0\,\gamma_a^2}{2}\cdot{}_2F_1\left(1, \frac{2}{s}; 1 + \frac{2}{s}; -z\right), \qquad z = \beta\left(\frac{\gamma_a}{\gamma_r}\right)^s
\end{aligned}
\tag{18.2}
```

(for $s = 1$ the closed form of the hyperbolic model,
$D = (2/\pi)\left[2(1+z)\left(z - \ln(1+z)\right)/z^2 - 1\right]$). It is zero at small strains, so a
small-strain viscous damping $\xi_{\min}$ (the DYNP damping at the smallest strain of the curve) is
added (SN-4) as NLDampType 1, a classical modal damping matrix
$C = M\,\Phi\operatorname{diag}(2\xi_n\omega_n)\,\Phi^{\mathsf{T}}M$ of the small-strain column on a
fixed base with strain-energy weighted modal ratios (frequency independent at the modal frequencies,
SN-5); 2, element dashpots $\eta/h$; or 3, Rayleigh $C = \alpha M + \beta K_0$ (matched at $f_1$ and
$5f_1$ when the multipliers are 0).

**Base and equation of motion.** BedInt 0: rigid base moving with the input (a *within* motion).
BedInt 1: elastic half-space by the Joyner and Chen (1975) dashpot $c = \rho_r V_r$ per unit area
with the *outcrop* motion as input ("currently disabled" in ACS SASSI V3, implemented here, SN-7).
Relative to the input frame:

```math
M\ddot u + C\dot u + R(u) + c_b\,\dot u_\text{base} = -M\,\mathbf{1}\,a_\text{in}(t) \tag{18.3}
```

integrated by the implicit Newmark average-acceleration method ($\gamma = 1/2$, $\beta = 1/4$) with
Newton-Raphson iterations and the consistent tangent of the Masing branches (force and displacement
norms ≤ 1e-6, at most 25 iterations, sub-steps bisected up to $2^{12}$ times; SN-10). Sub-increments:
NSTimeSunInc fixed steps, or (0) flexible steps of at most 1/400 s, repeated when a strain increment
exceeds 0.005 % (SN-9). Each sublayer is divided into an odd number of elements of thickness
$h \le V_s/(10 f_{\max})$, $f_{\max} = \min(25\,\text{Hz}, \text{Nyquist})$ (SN-11). Outputs:
`ACCxxx.TH`, `SNxxx.TH`, `SSxxx.TH` and FILE88 with the secant $G$ and the damping
$\xi_{\min}$ + Masing damping at the effective strain, so SITE can use the properties of the model
that was integrated (SN-12). Every soil sublayer needs an NLSLAYER set (Error 125), the half-space
does not.

**Verification.** VP-SN1: in the linear limit ($\gamma_r = 10^6\,\%$) the time-domain column
reproduces the frequency-domain linear SOIL solution of the same profile and motion (elastic base,
outcrop input:
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
isotropic elastic solid obeys
$\rho\ddot u = (K + 4G/3)\,\nabla(\nabla\cdot u) - G\,\nabla\times(\nabla\times u)$; for $G \to 0$ the
only stress is the pressure $p = -K\,\nabla\cdot u$, and

```math
\rho\,\ddot u = -\nabla p, \qquad \ddot p = c^2\,\nabla^2 p, \qquad c = \sqrt{K/\rho} \approx 1480\,\text{m/s} \tag{19.1}
```

the linear acoustic equation of an inviscid compressible fluid, written in displacements. Far below
the first compression frequency of the water column, $c/(4H)$ (74 Hz for $H = 5\,\text{m}$), the
water is incompressible, $\nabla^2 p = 0$, with $\partial p/\partial n = -\rho\,a_n$ on the walls and
the floor (imposed by the interface springs along the wall normals, frictionless along the walls) and
$p = 0$ at the free surface: the potential-flow problem of Westergaard (1933) and Housner (1963). The
water pushed by the walls adds the **impulsive mass**; for a rigid rectangular tank of length $2L$ in
the direction of motion and depth $H$ the exact series is

```math
\frac{m_i}{m} = \frac{2}{L H^2}\sum_n\frac{\tanh(\lambda_n L)}{\lambda_n^3}, \qquad \lambda_n = \frac{(2n + 1)\pi}{2H} \tag{19.2}
```

(exactly ½ for $L = H$; Housner's closed form $\tanh(\sqrt{3}\,L/H)/(\sqrt{3}\,L/H)$ is 6-11 %
conservative).

**Why $G = 10^{-8} K$ (D-W3-07).** D-WAT-01 proposed $\nu = 0.49$, i.e. $G = 0.0201\,K$. Such a
"water" is an elastic solid whose shear modes lie in the seismic band (about 10 Hz for 5 m of water);
below them it moves rigidly with the tank and gives the **total** mass instead of the impulsive mass
($F/(m a) = 1.00$ at 0.2-2 Hz for $L = H$, where the exact value is 0.50). The smaller $G$, the lower
the spurious shear modes: with $G = 10^{-8} K$ ($V_s = 10^{-4} V_p$) they lie below about 0.05 Hz for
metre-sized pools and the water behaves as an inviscid fluid from 0.2 Hz up; $K/G = 10^8$ is far from
the double-precision limit, and at every non-zero frequency the inertia controls the
near-zero-stiffness shear deformations. The water SOLIDs must use **incompatible modes** (FILLPOOL and
MERGEPOOL set `MOPT,0`): a trilinear hexahedron with 2×2×2 Gauss points locks volumetrically for a
nearly incompressible material (with `MOPT,1` $F/(m a)$ is again about 1.0 below 1 Hz and wrong at
higher frequencies). Oblique walls use the diagonal of the coupled interface spring matrix
$\mathrm{Stiff}\,P_N + \mathrm{stiff2}\,(I - P_N)$, because zero-length GENERAL elements are rejected
(Error 8, D-CHK-07). Sloshing (the convective mass, which needs the restoring force $\rho g$ at the
free surface) is not modelled; FILLPOOL prints Housner's estimate of it.

**Verification.** VP-54W: REFINEMODEL counts, areas and volumes; FILLPOOL water volume
$a\cdot b\cdot(H - m\cdot h)$, springs, offsets, the water material ($K = 2.2\,\text{GPa}$,
$G/K = 10^{-8}$), LISTPOOLINTER and MERGEPOOL. VP-W1: a rigid tank 10 m × 2 m filled to 5 m on a
practically rigid site, through SITE → POINT → HOUSE → ANALYS, with SHELL walls and with SOLID walls on
a slab: the water mass equals $\rho V$ (1e-10), the interface forces balance the inertia of the water
(1e-6), the tank follows the control motion (1e-3), and the base shear ratio $F/(m a)$ = 0.507-0.509
at 2-8 Hz equals the exact impulsive mass ratio 0.5 within 5 % (the +1.5 to +1.9 %
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
reference frame and the motion relative to it, $u = u_r + \iota\,u_g$, with $\iota\,u_g$ a rigid-body
translation (the free-field ground motion, or the motion of a reference node). A rigid-body
translation strains nothing ($K\iota = 0$), so with the interface DOFs $b$ prescribed the structural
DOFs $s$ obey

```math
M_{ss}\,\ddot u_{r,s} + C\,\dot u_{r,s} + K_{ss}\,u_{r,s} = -M_{ss}\,\iota_s\,a_g(t) - K_{sb}\,u_{r,b}(t) \tag{20.1}
```

* **Dynamic second step** (LOADGENDYN, method REL): an ANSYS transient analysis with `ACEL` =
  $a_g(t)$ (the acceleration of the reference frame; ANSYS applies $-M\cdot\mathrm{ACEL}$ to every
  mass) and the interface
  displacements relative to the ground prescribed with D (RELDISP `.THD` with the free-field reference,
  started at rest: the zero-mean offset of the periodic FFT solution is removed). LOADGEN writes them as
  APDL TABLE arrays. Method ACC fixes the interface and drives it with the absolute acceleration of a
  reference node (valid for a rigid foundation).
* **Equivalent static second step** (LOADGEN): freeze Eq. (20.1) at a critical time $t^*$ and drop
  the damping and velocity terms, $K_{ss}\,u_s = -M_{ss}\,a_s(t^*) - K_{sb}\,u_b(t^*)$: the nodal
  inertia forces of the **absolute** accelerations, $F = -m\,a(t^*)$ (MOTION `.ACC`), plus the
  interface displacements at $t^*$.
  Relative and total interface displacements differ by a rigid translation, which produces no stress.
  Without the displacements ("Acceleration") the interface is fixed: the classical fixed-base equivalent
  static analysis. The critical times are the largest local maxima of the base shear (or of a moment or a
  node response).

The masses are the row sums of the HOUSE structure mass matrix per direction,
$m_{i,d} = \sum_j M_{(i,d),(j,d)}$ (exact for lumped masses, the total mass per direction preserved
for consistent ones), or load nodes of the ANSYS model slaved to master nodes of the SSI model,
$a_k = a_M + \alpha_M \times (x_k - x_M)$ (a coarse stick → a refined ANSYS floor). Rayleigh damping
can be given by a damping ratio $\zeta$ at two frequencies:
$\alpha = 2\zeta\omega_1\omega_2/(\omega_1 + \omega_2)$, $\beta = 2\zeta/(\omega_1 + \omega_2)$.

**What is approximate.** For a linear structure with the same mass, stiffness and damping as the SSI
model the decomposition is exact. The dynamic step replaces the frequency-independent hysteretic
damping of SASSI by Rayleigh damping, $\zeta(\omega) = \alpha/(2\omega) + \beta\omega/2$, exact at two
frequencies only ("a significant
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

1. **Target**: a design spectrum (Hz, SA in g) at damping $\zeta$, interpolated log-log; constant
   spectral displacement below its first frequency, constant SA above its last.
2. **Initial motion**: a stationary Gaussian process with random phases (PCG64, seed `<rand>`) and
   Fourier amplitudes from the SIMQKE power-spectral-density estimate (Gasparini and Vanmarcke 1976),
   multiplied by an envelope with a stationary strong-motion part; or the Fourier phases of a seed
   record (`<accopt>` = 1).
3. **Frequency-domain matching** (Levy-Wilkinson): Fourier amplitudes multiplied by
   $\mathrm{SA}_\text{target}/\mathrm{SA}_\text{computed}$ (Nigam-Jennings spectra at 100 points per
   decade), phases kept; with random
   phases each iteration also clips the stationary process at the target zero-period acceleration
   (iterative clipping and filtering), which lowers the crest factor so that the PGA can equal the
   target ZPA.
4. **Time-domain refinement** (P1): Al Atik and Abrahamson (2010) tapered-cosine wavelets added at the
   response-peak times to move the peaks to the target.
5. **Drift control and baseline correction** in every step (D-EQK-03, D-W1-09): a zero-phase high-pass
   filter on the frequency-domain updates, zero-mean acceleration, and a degree-5 displacement
   polynomial fitted by constrained least squares (zero final velocity and displacement).
6. **Checks** (D-EQK-04): the manual's criteria and SRP 3.7.1 Rev. 4 Option 1 Approach 2, reported
   separately: duration $\ge 20\,\text{s}$; $\Delta t \le 0.005\,\text{s}$ (manual) and Nyquist
   $\ge 50\,\text{Hz}$ (SRP); spectrum at $\ge 100$ points per decade never more than 10 % below or
   30 % above the target, at most 9 adjacent points below; strong-motion duration (Arias 5-75 %)
   $\ge 6\,\text{s}$; cross-correlation $\lvert\rho\rvert \le 0.16$ between components;
   $\mathrm{PSD} \ge 80\,\%$ of the target PSD over 0.3-24 Hz when a target PSD is given.

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
