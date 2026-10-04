# R1: Theory of the Original SASSI Method (implementation-grade note)

**Purpose.** This note rebuilds the equations of the Berkeley SASSI programs (Lysmer et al. 1981; SASSI2000,
1999) and their SHAKE front end. It covers the modules SITE, POINT, ANALYS, MOTION and SOIL, the
damping convention, and the frequency-interpolation scheme. The level of detail is enough to code each step.
It complements the ACS SASSI manual specs `02_theory_modules.md` and `05a`–`05d`. Those specs record what
the ACS manual *says*. This note records the *underlying formulation*, which the manual does not give.

**Status of sources.** The two primary documents are the SASSI report UCB/GT/81-02 and the SASSI2000
Theory Manual. Both are held behind the paid NISEE e-library and **were not accessible** for this note
(the access attempt returned a membership page). The formulation is therefore rebuilt from:

* openly available primary papers by the SASSI authors and their close collaborators: Tajirian &
  Tabatabaie 1985 [3], Kausel 1981 [4], TLUSH 1981 [7], the SHAKE 1972/1973 report and code [8, 9], and
  Ostadan et al. 2004 [10];
* NRC/DOE documents that transcribe SASSI2000 equations [5, 6, 13–17];
* the ACS SASSI manual [18];
* numerical checks run for this note. The scripts are in `docs/spec/R1_checks/` and are re-run with
  `python run_all.py`.

**Tag legend (used on every subsection).**

| Tag | Meaning |
|---|---|
| **[S#]** | Confirmed by source # (see §11). Equations marked [S] are copied from that source, with only the notation adapted as stated. |
| **[R]** | Reconstruction or engineering decision by the author of this note. It is not stated in any source that was accessed. |
| **[V]** | Checked numerically in this study (`R1_checks/`). §9 gives the test and its result. |

---

## 0. Conventions (read first)

| Item | Convention used in this note | Status |
|---|---|---|
| Time factor | $e^{+i\omega t}$ | [S8, S9, S4]. SHAKE, Kausel and Waas all use $e^{i\omega t}$. The SASSI damping sign (Im $G^*>0$, §1) also requires it. [5] Eq. 2-3 prints $e^{-i\omega t}$, which is taken to be a transcription slip. |
| Horizontal propagation | $e^{-ikx}$. Outgoing waves have $\mathrm{Im}\,k<0$, or $\mathrm{Re}\,k>0$ when $k$ is real. | [S4] p.15 ("Following Waas, we choose…") |
| Outgoing cylindrical waves | $H^{(2)}_\mu(k\rho)$ (second-kind Hankel) | [S4] Eq. 88–89 |
| Vertical axis | **Kausel form: $z$ up.** Interfaces are numbered 1 (ground surface) to $N$ downward. The in-plane vector is $\{u_x,\, i u_z\}$, the "i-scaled vertical", which makes all layer matrices symmetric. | [S4] Eq. 1, Fig. 1 |
| Waas/SASSI form | $z$ (here called $v$) **down**, no $i$ factor. $B$ is skew-symmetric and enters as $iBk$. The two forms are algebraically identical (§2.3). | [S4] Appendix, [S5] Eq. 2-8 |
| SITE local axes | $z'$ vertical up. $x'$ lies in the vertical plane of propagation. $y'$ completes the right-handed set (SH/Love motion). | [S18] ch. 6 SITE options |
| Damping | $G^* = G\,(1-2\beta^2+2i\beta\sqrt{1-\beta^2})$ (§1) | [S7, S8, S9]. Inferred for SASSI from [10] and [V]. |
| Units | Any consistent set. $\rho$ is mass density. SASSI input uses unit weight and $g$. | — |

---

## 1. Complex-modulus (hysteretic) damping convention  — task item (4)

**1.1 The form used by the Berkeley codes.**

$$
G^* = G\left(1-2\beta^2+2i\beta\sqrt{1-\beta^2}\right),\qquad
V^* = \sqrt{G^*/\rho} = V\left(\sqrt{1-\beta^2}+i\beta\right),\ |V^*|=V .
$$

* **[S8]** The original SHAKE (Schnabel, Lysmer & Seed 1972, Eq. 7) used $G^*=G(1+2i\beta)$. The
  Udaka & Lysmer revision of **September 1973**, bound into the same report (p. S-2), says the purpose
  of the change was "to redefine the complex modulus from G* = G(1+2iβ) to G* = G(1−2β²+i2β√(1−β²))".
  The same revision appears as a comment in the 1973 SHAKE source listing.
* **[S9]** The SHAKE91 source (subroutine `CXSOIL`, public domain) codes exactly
  `GREAL = GL*(1-2*BL**2)`, `GIMAG = 2*BL*GL*SQRT(1-BL**2)`.
* **[S7]** TLUSH (Kagawa, Mejia, Seed & Lysmer, UCB/EERC-81/14, 1981) is from Lysmer's group in the same
  year as SASSI and acknowledges Tajirian. Its Eq. (3) is $G^*=G(1-2\beta^2+2i\beta\sqrt{1-\beta^2})$.
* **Evidence for SASSI itself [S10 + V].** Ostadan, Deng & Roesset (2004), Table 1, ran a fixed-base
  SDOF in SASSI2000 with β = 5/10/15/20 %. They back-calculated "1/(2·Umax)" = 5.0/9.9/14.7/19.5 %. For
  the total-acceleration transfer function $k^*/(k^*-m\omega^2)$ this note computes
  1/(2|H|max) = 4.99/9.95/14.83/19.60 % for the Udaka–Lysmer form, and 4.98/9.81/14.37/18.57 % for
  $1+2i\beta$ (`test9_damping.py`). Only the Udaka–Lysmer form reproduces the published 19.5 % at
  β = 20 %. **Conclusion:** SASSI2000 uses the Udaka–Lysmer form. Confidence is high, but the evidence
  is indirect.
* **[R]** Apply the factor $c(\beta)=1-2\beta^2+2i\beta\sqrt{1-\beta^2}$ separately to the shear
  modulus (with $\beta_s$) and to the constrained modulus $M=\lambda+2G=\rho V_p^2$ (with $\beta_p$). The
  SASSI layer card carries Vp, Vs, βp and βs separately ([S18] `L` command). Then
  $\lambda^*=M^*-2G^*$. Use the same factor for HOUSE element moduli and springs. Keep $1+2i\beta$ only as
  a selectable benchmark option.

**1.2 Properties to keep in mind.**

* The loss angle is $\delta$ with $\sin\delta = 2\beta\sqrt{1-\beta^2}$ and $\cos\delta = 1-2\beta^2$, so
  $G^*=G e^{i\delta}$ with $\delta = 2\arcsin\beta$. The magnitude stays at $G$. With $1+2i\beta$ the
  magnitude becomes $G\sqrt{1+4\beta^2}$, a slight stiffening.
* For a vertically propagating wave, $k^*=\omega/V^* = (\omega/V)(\sqrt{1-\beta^2}-i\beta)$. The amplitude
  decays as $e^{-\beta\omega z/V}$ and the phase velocity is $V/\sqrt{1-\beta^2}$.
* **Mapping to ANSYS harmonic analysis [R]** (for benchmarking against ANSYS). ANSYS structural damping
  gives $K(1+ig)$. To reproduce the SASSI form, set $E_{\rm ANSYS}=E(1-2\beta^2)$ and
  $g = 2\beta\sqrt{1-\beta^2}/(1-2\beta^2)$. Anderson & Ostadan [15] used "ANSYS Full Harmonic Response"
  with frequency-independent damping as an independent SASSI benchmark.

---

## 2. SITE module: thin-layer (discrete-layer) formulation  — task item (1)

### 2.1 Discretization **[S5, S18, S4]**

* The site is horizontal layers on a rigid base or on a viscoelastic halfspace ([S5] §2.2.2). Within a
  layer the displacement **varies linearly with depth** ([S5]: "The displacements in each soil layer are
  assumed to vary linearly with the thickness"; [S4] Eq. 12).
* **All interaction nodes must lie on layer interfaces** ([S5] §2.2.2.1; [S18] 1.5.1). The HOUSE
  excavated-soil mesh planes must therefore coincide with the SITE interfaces (the embedment layers).
* Layer-thickness rule: $h \le \lambda_{\min}/5$ with $\lambda_{\min}=V_s/f_{\rm cut}$ [S18 §4 item 12].
  The 1/5 rule is tied to the **mixed (½ lumped + ½ consistent) mass**. Lysmer & Kuhlemeyer (1969)
  report better than 10 % accuracy with $h\le\lambda/8$ for lumped *or* consistent mass, and
  $h\le\lambda/5$ for the mixed mass ([S5] §4.2, citing Lysmer & Kuhlemeyer 1969 and Lysmer et al.
  1999b). SASSI2000 forms mass matrices as ½ lumped + ½ consistent ([S5] Eq. 2-2 text). The ACS manual
  says the same for structural/solid elements [S18 1.5.1 item 12].
* **[R]** Use the mixed mass for the SITE layers as well, i.e. $M_{\rm layer}=\tfrac12M_{\rm cons}+\tfrac12M_{\rm lump}$.
  [V]: at $h=\lambda/5$ the mixed mass keeps the 1-D amplification within about 9 %, while the
  consistent mass errs by up to 17 % (`test4b.py`, last two columns). The same mass choice must be used
  in POINT (§3) and in the HOUSE excavated soil.

### 2.2 Layer matrices (exact)  **[S4] Table 1**, notation adapted

For layer $m$ of thickness $h$, with Lamé $\lambda^*$, shear modulus $G^*$ (complex, §1) and density
$\rho$, take the dof order **$\{u_{x1},u_{y1},\tilde u_{z1},u_{x2},u_{y2},\tilde u_{z2}\}$**. Here
$1$ = top interface, $2$ = bottom interface, and $\tilde u_z = i\,u_z$ with $z$ up. Dots are zeros.

$$
A_m=\frac h6\begin{bmatrix}
2(\lambda+2G)&\cdot&\cdot&\lambda+2G&\cdot&\cdot\\
\cdot&2G&\cdot&\cdot&G&\cdot\\
\cdot&\cdot&2G&\cdot&\cdot&G\\
\lambda+2G&\cdot&\cdot&2(\lambda+2G)&\cdot&\cdot\\
\cdot&G&\cdot&\cdot&2G&\cdot\\
\cdot&\cdot&G&\cdot&\cdot&2G\end{bmatrix},\qquad
B_m=\frac12\begin{bmatrix}
\cdot&\cdot&\lambda-G&\cdot&\cdot&-(\lambda+G)\\
\cdot&\cdot&\cdot&\cdot&\cdot&\cdot\\
\lambda-G&\cdot&\cdot&\lambda+G&\cdot&\cdot\\
\cdot&\cdot&\lambda+G&\cdot&\cdot&-(\lambda-G)\\
\cdot&\cdot&\cdot&\cdot&\cdot&\cdot\\
-(\lambda+G)&\cdot&\cdot&-(\lambda-G)&\cdot&\cdot\end{bmatrix}
$$

$$
G_m=\frac1h\begin{bmatrix}
G&\cdot&\cdot&-G&\cdot&\cdot\\
\cdot&G&\cdot&\cdot&-G&\cdot\\
\cdot&\cdot&\lambda+2G&\cdot&\cdot&-(\lambda+2G)\\
-G&\cdot&\cdot&G&\cdot&\cdot\\
\cdot&-G&\cdot&\cdot&G&\cdot\\
\cdot&\cdot&-(\lambda+2G)&\cdot&\cdot&\lambda+2G\end{bmatrix},\qquad
M_m^{\rm cons}=\frac{\rho h}6\begin{bmatrix}2I_3&I_3\\ I_3&2I_3\end{bmatrix},\quad
M_m^{\rm lump}=\frac{\rho h}2\begin{bmatrix}I_3&0\\0&I_3\end{bmatrix}.
$$

The global matrices are $A=\{A_m\}$, $B=\{B_m\}$, $G=\{G_m\}$ and $M=\{M_m\}$, assembled by overlapping
the layer blocks at shared interfaces exactly as in FE assembly [S4 Fig. 3]. The interface load–
displacement relation in the frequency–wavenumber domain is ([S4] Eq. 11, 13a)

$$
\bar P = \big(A k^2 + B k + G - \omega^2 M\big)\,\bar U .
$$

For a rigid base, delete the bottom-interface dofs. For a halfspace, see §2.4 (dashpots are added to
$G$).

*Derivation check [R][V].* The same matrices follow from the weak form with $u=N(z)U\,e^{-ikx}$ and
linear $N$. The $k^2$ terms give $\int N^TN$ weighted by $(\lambda+2G, G, G)$. The $k^0$ terms give
$\int N'^TN'$ weighted by $(G,G,\lambda+2G)$. The $k^1$ terms give the $B$ coupling. Modal-versus-direct
inversion agrees to $10^{-14}$ (`tlm.py`).

**Sub-blocks used below** (order by direction, as in [S4] Eq. 15). $A_x$ uses $(\lambda+2G)$; $A_z$ and
$A_y$ use $G$. $G_x$ and $G_y$ use $G$; $G_z$ uses $\lambda+2G$. $B_{xz}$ is the $x$-row/$z$-column block
of $B$. $C_\bullet = G_\bullet-\omega^2M$.

**SH/Love (antiplane, 1 dof/node)** [S4, S5 Eq. 2-7/2-9]: $\big(A_y k^2 + G_y - \omega^2 M\big)\phi_y = 0$,
with $A_y=\frac{Gh}6\begin{bmatrix}2&1\\1&2\end{bmatrix}$ and $G_y=\frac Gh\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$
per layer.

### 2.3 Equivalence with the Waas/SASSI form **[S4 Appendix, S5 Eq. 2-8] [V]**

Waas (1972), and hence SASSI SITE, write ([S5] Eq. 2-8)

$$
\big([A]k^2 + i[B_W]k + [G] - \omega^2[M]\big)\{V\}=0,
$$

with vertical displacement $v$ positive **down** and no $i$ factor. In the order
$\{u_1,v_1,u_2,v_2\}$, $A$, $G$ and $M$ are as above restricted to the in-plane dofs, and

$$
B_W=\frac12\begin{bmatrix}
\cdot&-(\lambda-G)&\cdot&\lambda+G\\
\lambda-G&\cdot&\lambda+G&\cdot\\
\cdot&-(\lambda+G)&\cdot&\lambda-G\\
-(\lambda+G)&\cdot&-(\lambda-G)&\cdot\end{bmatrix}\quad(\text{skew-symmetric}).
$$

With $v = i\,\tilde u_z$ (so that $v_{\rm down}=-u_{z,\rm up}$) and $T={\rm diag}(1,i,1,i)$:
$T^{-1}\,L_W\,T = L_{\rm Kausel}$ holds exactly (`test10_waas_kausel.py`, residual $0.0$). An
implementation can use either form. This note uses Kausel's form, because the explicit Green's functions
(§3.2) are published in it.

### 2.4 Halfspace: variable-depth method + viscous base **[S18, S5, S3, S12]**

* **Confirmed.**
  - SITE adds $n$ extra sublayers with the **halfspace properties** below the user layers. Their
    **total thickness is $1.5\lambda_s$**, with $\lambda_s = V_{s,\rm hs}/f$, so the buffer depth varies
    with frequency. The sublayers increase in thickness with depth and become thinner as frequency
    rises. **Viscous dashpots** in the horizontal and vertical directions are attached at the base of
    the extended system ([S18] §4 items 8 and 20; "Rigid Base Rock vs. Halfspace Condition").
  - The choice of 1.5λ comes from the decay of the fundamental Rayleigh mode with depth ([S18], [S5]
    §4.2).
  - The recommended $n$ is 10–20, at most 20; $n=0$ means a rigid base ([S18]). SASSI2000 requires at
    least 4 sublayers. Their thickness depends on the total thickness, the number of sublayers and the
    thickness of the soil layer above the halfspace ([S5] §2.2.2.2.2).
  - Dashpots per unit area (Lysmer–Kuhlemeyer): $c_p=\rho V_p$ and $c_s=\rho V_s$ ([S5] Eq. 2-13/2-14).
  - Ghiocel [S12] calls this arrangement "buffer layers plus Lysmer–Kuhlemeyer viscous boundaries".
  - The halfspace is used both in the SITE eigenproblem (Mode 1) and in POINT ([S3]: "a halfspace can be
    simulated by using the variable depth and viscous boundary methods").
* **[R] Implementation.**
  1. With $e^{i\omega t}$, a dashpot $c$ adds $+i\omega c$ to the diagonal of $G$ at the bottom
     interface: $c_s$ on the $x$ and $y$ dofs and $c_p$ on the $z$ dof. In i-scaled variables the
     $z$ entry is unchanged, because the $i$ factors cancel. Use $V^*$ from the complex moduli. For 1-D
     vertical waves, $\rho V^*$ is then the *exact* impedance of a damped halfspace. Using real $V$
     instead is a negligible difference below 1.5λ of damped material.
  2. Sublayer thicknesses: a geometric series $h_i = h_1 q^{i-1}$, $i=1..n$, with
     $\sum h_i = 1.5\lambda_s$, $h_1 = \min\{h_{\rm last\ soil}\cdot V_{s,\rm hs}/V_{s,\rm last},\ 1.5\lambda_s/n\}$,
     and $q\ge1$ solved by bisection. Cap every $h_i\le\lambda_s/8$. If the cap is violated, fall back to
     uniform $h_i = 1.5\lambda_s/n$.
     [V] (`test4b.py`, `test4c.py`; vertical SH; outcrop TF compared with exact SHAKE at 0.5–15 Hz):
     - Uniform sublayers: ≤0.7 % amplitude error for $n=20$, ≤1.5 % for $n=10$, up to 20 % for $n=5$.
     - Geometric series started from the thin last soil sublayer: ≤4.6 % for $n=20$, but up to 26 % for
       $n=10$ and up to 140 % for $n=5$, because the deepest sublayers become too thick.
     - Geometric series started from the velocity-scaled sublayer: ≤7 % for $n=10$.

     This is consistent with the ACS advice to use 20 sublayers and the SASSI2000 statement that
     "10 extra half-space layers are sufficient" for most cases [S5 §5.4].

### 2.5 Mode 1: Rayleigh and Love eigenproblems **[S4] Eq. 14–22, [S5] Eq. 2-8/2-9, [S18]**

Set $\bar P=0$. With $C=G-\omega^2M$ (including the dashpots), solve, per frequency:

* **Rayleigh (generalized P-SV)** $(A k^2 + Bk + C)\phi=0$, a quadratic eigenproblem of size $2N_f$, where
  $N_f$ is the number of free interfaces. It is linearized *without* doubling the size ([S4] Eq. 16):

$$
\left(k^2\underbrace{\begin{bmatrix}A_x&0\\B_{xz}^T&A_z\end{bmatrix}}_{\bar A}+\underbrace{\begin{bmatrix}C_x&B_{xz}\\0&C_z\end{bmatrix}}_{\bar C}\right)
\underbrace{\begin{Bmatrix}\phi_x\\k\phi_z\end{Bmatrix}}_{Z}=0 .
$$

  This is a generalized linear eigenproblem in $k^2$ with $2N_f$ eigenvalues. Take
  $k_j=\sqrt{k_j^2}$ on the branch with $\mathrm{Im}\,k_j<0$, or $\mathrm{Re}\,k_j>0$ if $k_j$ is real
  (decaying or outgoing in $+x$) [S4 p.15]. Recover $\phi_z = Z_{\rm lower}/k_j$.
* **Normalization** [S4 Eq. 22a]. With $Y_j=\{k_j\phi_{xj};\phi_{zj}\}$ and $Z_j$ as above:
  $Y_j^T\bar A Z_j = k_j$, i.e.
  $\phi_x^TA_x\phi_x + \phi_z^TA_z\phi_z + k_j^{-1}\phi_z^TB_{xz}^T\phi_x = 1$.
  Kausel shows this is the same normalization Waas used [S4 Appendix]. Orthogonality is
  $Y_i^T\bar A Z_j=0$ for $i\neq j$ [S4 Eq. 21].
* **Love** $(A_yk^2+C_y)\phi_y=0$ gives $N_f$ modes, normalized by $\phi_y^TA_y\phi_y=1$.
* **Matrix sizes** ([S5]): for $n$ layers, Rayleigh is $2n\times2n$ and Love is $n\times n$. The
  body-wave systems (§2.7) include the halfspace interface, giving $2(n+1)$ and $(n+1)$.
* **Modal flexibility in the wavenumber domain** [S4 Eq. 49–52] [V, machine precision], with
  $D=\mathrm{diag}\,(k^2-k_j^2)^{-1}$ and $K_R=\mathrm{diag}(k_j)$:

$$
F_{xx}=\Phi_xD_R\Phi_x^T,\quad F_{xz}=k\,\Phi_xK_R^{-1}D_R\Phi_z^T=F_{zx}^T,\quad F_{zz}=\Phi_zD_R\Phi_z^T,\quad F_{yy}=\Phi_yD_L\Phi_y^T .
$$

  Two useful identities are $\Phi_xK_R^{-1}\Phi_z^T=0$ [S4 Eq. 34] and the factorization above
  (`tlm.py`, error $10^{-14}$).
* **What SITE stores in FILE2** [S18] (contents [R]): for each frequency, the generated sublayer
  thicknesses and properties, $\{k^R_j,\phi_x^j,\phi_z^j\}_{j=1}^{2N_f}$ and
  $\{k^L_l,\phi_y^l\}_{l=1}^{N_f}$. One FILE2 serves 2-D and 3-D problems, because the eigenproblems are
  identical [S18; S4 p.13: "the stiffness matrix for cylindrical coordinates is identical to that of the
  plane strain case"].
* [V] For a uniform deep stratum with ν = 0.25 at $h=\lambda/20$, the fundamental mode gives
  $c_R/V_s$ = 0.9210 (consistent mass), 0.9194 (lumped) and 0.9202 (mixed). The exact value is 0.9194
  (`test2.py`).

### 2.6 Plane-strain transmitting boundary (Waas–Lysmer) **[S5 Eq. 2-11/2-12, S4 p.23] [V]**

For a vertical boundary at $x=x_b$ with the semi-infinite layered region at $x>x_b$, the boundary nodal
forces are $\{P\}=[R]\{U\}$ with $[R]=i[A][V][K][V]^{-1}+[D]$ ([S5] Eq. 2-12). Here $D$ is a constant
matrix built from $\lambda$ and $G$. A form that can be coded directly and has been verified
(physical variables, $z$ up) is given below.

* Outgoing mode $j$: $u_x=\phi_{xj}e^{-ik_j(x-x_b)}$ and $u_z=-i\,s\,\phi_{zj}e^{-isk_j(x-x_b)}$. Here
  $s=+1$ for a right boundary and $s=-1$ for a left boundary (outgoing toward $-x$).
* Consistent boundary forces on the interior, per unit $y$, use $E_c=\int cN^TN\,dz$ (per layer
  $c\,h/6\,[\,2\ 1;1\ 2\,]$) and $Q_c=\int cN^TN'\,dz$ (per layer, $z$ up, top node first,
  $c/2\,[\,1\ {-1};1\ {-1}\,]$):

$$
f_x = s\left(E_{\lambda+2G}\,\partial_xu_x + Q_\lambda u_z\right),\qquad f_z = s\left(Q_G u_x + E_G\,\partial_x u_z\right),\qquad \partial_x\to -isk_j .
$$

* Stack mode columns: $\Psi=[u_x;u_z]$ and $\mathcal T=[f_x;f_z]$. Then $R=-\mathcal T\Psi^{-1}$ is added
  to the interior FE matrix on the boundary dofs. In Kausel variables this is $iA\Phi K\Phi^{-1}$ plus a
  constant $D$ term, i.e. the [S5] form.
* [V] A plane-strain FE strip ($-R_0..R_0$, 2 elements) with two such boundaries reproduces Kausel's
  line-load Green function within 0.3 % for $|x|\ge2R_0$ (`point2.py`).

### 2.7 Mode 2: free-field motion and control-point normalization **[S5, S18] + [R]**

SASSI uses one free-field model for all interaction nodes. Its motion is a superposition of plane
**body waves** (P, SV, SH, vertical or inclined) and **surface waves** (Rayleigh, Love) ([S5] §2.2.2;
[S18] SITE options). For each wave type $w$ and frequency, SITE builds a unit-normalized field
$\hat U^{(w)}(z,x')$. It then superimposes them with participation ratios $r_w(\omega)$. The ratios are
given at two frequencies, interpolated linearly, and must sum to 1 at every frequency [S18]. Eq. 2-10 of
[S5] gives the horizontal dependence as $\{U(x)\}=\delta\{U\}e^{-ikx}$, where $\delta$ is a "mode
participation factor obtained from the input control motion".

**(a) Vertically propagating SV ($x'$), SH ($y'$), P ($z'$) — the standard nuclear case.** Set $k=0$ in
§2.2; the problem then decouples into 1-D columns. Shear (SV/SH) uses $G_s$ (built with $G^*$) and
$c=c_s$. P uses $G_p$ (built with $\lambda^*+2G^*$) and $c=c_p$. For an incident (upgoing) wave of
displacement amplitude $E_b$ at the dashpot level, the column equation is [R; derivation below; V]:

$$
\big(G_\bullet-\omega^2M+i\omega c\,e_Ne_N^T\big)\,u = 2\,i\omega c\,E_b\,e_N .
$$

Derivation: the halfspace traction on the column is $t=i\omega\rho V^*(u_I-u_R)=i\omega\rho V^*(2u_I-u)$.
In general form, $(K_{\rm col}+K_{\rm out})u=(K_{\rm out}+K_{\rm in})u_I$, which reduces to the line
above when $K_{\rm out}=K_{\rm in}=i\omega\rho V^*$. A rigid base instead prescribes $u_N=1$.

* **Normalization to the control point [S18][R].** The control motion is specified at the **top of the
  control layer $L_{cp}$**, in a direction $d$. The free-field transfer function is
  $\hat U(z)=u(z)/u_d(z_{cp})$, so the control point moves with unit amplitude. This is a ratio of
  displacements, so the same function serves for accelerations.
* **Within vs outcrop.** SITE takes the control motion as an in-column ("within") motion at a layer
  interface [S18 SITE options]. The outcrop flag exists in SOIL [S18 `SOIL`/`SACC` commands].
  Recommended workflow [R]: when the design motion is an outcrop motion, run SOIL first to get the
  within motion at the control point, or define the control point at the surface. If an outcrop control
  point is needed inside SITE, compute it as "surface of the truncated column". [V]
  (`test4b.py`): the discrete outcrop motion at the halfspace top equals the surface motion of the
  halfspace-only column (same sublayers and dashpot, same load). It matches SHAKE's $2E$ within 4.6 % for
  $n=20$. *Caution*: away from the halfspace, SHAKE's "outcrop = 2E" differs from the true outcrop of a
  truncated profile [S17 §2.1–2.2].
* [V] The within TF (surface over top of halfspace) **does not depend on the halfspace model**: results
  for $n$ = 5, 10 and 20 are identical. It converges to SHAKE as $h$ decreases. With consistent mass the
  maximum error near resonant peaks is 20 % at $h=\lambda/10$ and 5.6 % at $h=\lambda/20$ (`test4.py`).
  The halfspace model matters for motions **below** the control point.

**(b) Inclined body waves [R].** The source excerpts give only the form: [S5] Eq. 2-6 is
$([A]k^2+[\bar B]k+[G]-\omega^2[M])\{U\}=\{0;P_b\}$, with "$P_b$ and $k$ determined by the incidence
angle and wave type". Reconstruction:

* $k=\omega\sin\theta_w/V_{w,\rm hs}$, with $\theta_w$ the angle from the vertical $z'$ axis [S18] and
  $V_{w,\rm hs}$ the halfspace P or S velocity.
* Solve $(Ak^2+Bk+C+K_{\rm out}(k))\,U=(K_{\rm out}+K_{\rm in})\,u_I$ at the base. Here $K_{\rm out}$ is
  either the dashpot $i\omega\,{\rm diag}(c)$ (approximate for oblique incidence) or the exact halfspace
  stiffness of Kausel & Roesset (1981). $u_I$ is the analytic incident P/SV/SH plane wave in the
  halfspace material (amplitude 1), evaluated at the base.
* The field is $\hat U(z)\,e^{-ik(x'-x'_{cp})}$, normalized to the control point as in (a).
* Do not use the variable-depth dashpot base for large incidence angles; prefer the exact halfspace
  stiffness.

**(c) Surface waves [S5 Eq. 2-10, S18] + [R].**

* The field is $\hat U(z,x')=\phi_j(z)\,e^{-ik_j(x'-x'_{cp})}/\phi_{j,d}(z_{cp})$. For Kausel mode
  shapes the physical vertical component is $u_z=-i\,\phi_z$.
* Mode selection: "**Shortest wavelength**" = the propagating mode with the largest $\mathrm{Re}\,k_j$.
  "**Least decay**" = the mode with the smallest $|\mathrm{Im}\,k_j|$. Both options are named in [S18];
  their exact definitions are [R].

**(d) Transfer to the model axes.** ANALYS maps $x'y'z'$ to the structure's $xyz$ ([S18]: "The
transformation … will be done by the ANALYS module"). [R]: $U_{xyz}=R_z(a)\,U_{x'y'z'}$, with the
phase evaluated at $x'_{\rm node}=x\cos a+y\sin a$.

**(e) FILE1 [S18] (contents [R]).** For each frequency: the normalized free-field vector at every
interface (3 complex components), the wavenumber and type of each participating wave, and the ratios.

---

## 3. POINT module: soil flexibility from point-load solutions  — task item (2)

### 3.1 What POINT solves **[S3, S5, S18, S11]**

* "The basic problem in determining the dynamic flexibility matrix is to find the response of a layered
  halfspace to a harmonic point load. Each column of the flexibility matrix is formed by applying a unit
  point load at the interaction degree-of-freedom associated with that column and by computing the
  resulting displacements at all the interaction nodes" [S3].
* The model is axisymmetric: "a central core of special cylindrical axisymmetric finite elements
  connected at the perimeter to a semi-infinite layered zone which is represented by axisymmetric
  transmitting boundaries". Displacements are obtained "both at the central nodes and at any point
  outside the cylindrical elements … These displacements which are computed in a cylindrical coordinate
  system are transformed to the global cartesian coordinate system" [S3].
* Unit loads are applied successively at the interaction nodes on the **centre line** of a single
  column of loaded nodes. "The solution … for a single column of the interaction nodes can be used as the
  compliance matrix for the remaining columns" [S5 §2.2.3], i.e. by horizontal translation [S11].
* Loads are needed at interfaces $1..L+1$, where $L$ is "Last Layer Number in Near Field Zone" [S18].
  The soil column is "a FE model with plane strain elements for 2D models, or axisymmetrical elements for
  3D models" [S18 §4.1.2 item 17].
* SC-SASSI [S11] notes that the traditional formulation is not defined directly under a concentrated
  load, where the continuum displacement is infinite. The FE core is what gives a finite, mesh-dependent
  self-flexibility.

### 3.2 Exact TLM Green's functions (reference and far-field oracle)  **[S4] §2.6, Eq. 61, 70, 88, 89** [V]

Cylindrical Fourier decomposition [S4 Eq. 5–8]:

* symmetric case: $u_\rho=\tilde u_\rho\cos\mu\theta$, $u_\theta=-\tilde u_\theta\sin\mu\theta$,
  $u_z=\tilde u_z\cos\mu\theta$;
* vertical loads excite $\mu=0$; horizontal ($x$) loads excite $\mu=1$.

For a **unit point load at interface $n$**, the displacement at interface $m$ at radius $\rho$ is
([S4] Eq. 88–89, $z$ up, $H\equiv H^{(2)}$):

*Horizontal load $P_x=1$ ($\mu=1$):*

$$
\tilde u_\rho=\frac1{4i}\Big[\sum_{l=1}^{2N}\phi_x^{ml}\phi_x^{nl}\,\frac{d}{d\rho}\frac{H_1(k^R_l\rho)}{k^R_l}+\frac1\rho\sum_{l=1}^{N}\phi_y^{ml}\phi_y^{nl}\frac{H_1(k^L_l\rho)}{k^L_l}\Big],\quad
\tilde u_\theta=\frac1{4i}\Big[\frac1\rho\sum_{l}^{2N}\phi_x^{ml}\phi_x^{nl}\frac{H_1(k^R_l\rho)}{k^R_l}+\sum_{l}^{N}\phi_y^{ml}\phi_y^{nl}\frac{d}{d\rho}\frac{H_1(k^L_l\rho)}{k^L_l}\Big],
$$
$$
\tilde u_z=-\frac1{4i}\sum_{l}^{2N}\phi_z^{ml}\phi_x^{nl}H_1(k^R_l\rho).
$$

*Vertical load $P_z=1$ (upward, $\mu=0$):*
$\ u_\rho=\frac1{4i}\sum_l\phi_x^{ml}\phi_z^{nl}H_1(k^R_l\rho)$, $\ u_\theta=0$,
$\ u_z=\frac1{4i}\sum_l\phi_z^{ml}\phi_z^{nl}H_0(k^R_l\rho)$.

With the normalization of §2.5, these are **true physical displacements** ($u_z$ up). [V] In the static
limit on a deep graded stratum they match Boussinesq and Cerruti within 0.1–1.7 % at
$\rho$ = 1–5 m, including all signs (`test3.py`). Example: at the surface, ahead of a horizontal push
($\theta=0$), the ground moves *down*, with $u_z=-(1-2\nu)P/(4\pi G\rho)$ upward-positive. This follows
from Boussinesq by reciprocity.

**Uniform disk load** of radius $R$ and intensity $q$ [S4 Eq. 61, 70, Table 2]. Example: vertical disk,
$u_z=qR\sum_l\phi_z^{ml}\phi_z^{nl}I^R_{1l}$. The closed-form integrals are

$$
I_{1l}=\int_0^\infty\frac{J_0(k\rho)J_1(kR)}{k^2-k_l^2}dk=\begin{cases}\frac{\pi}{2ik_l}J_0(k_l\rho)H_1(k_lR)-\frac{1}{Rk_l^2},&\rho\le R\\ \frac{\pi}{2ik_l}J_1(k_lR)H_0(k_l\rho),&\rho\ge R\end{cases},\qquad
I_{3l}=\int_0^\infty\frac{J_1(k\rho)J_1(kR)}{k(k^2-k_l^2)}dk=\begin{cases}\frac{\pi}{2ik_l^2}J_1(k_l\rho)H_1(k_lR)-\frac{\rho}{2Rk_l^2}\\ \frac{\pi}{2ik_l^2}J_1(k_lR)H_1(k_l\rho)-\frac{R}{2\rho k_l^2}\end{cases},
$$

together with $I_2=-dI_1/d\rho$ and $I_4$ [S4 Table 2]. (The scan of the $\rho\ge R$ branch of $I_3$
shows a leading "−"; continuity at $\rho=R$ and the point-load limit both require "+".) Average
displacement under the disk:
$u_{z,\rm avg}=2q\sum\phi_z\phi_zI^R_{3l}|_{\rho=R}$ and
$u_{x,\rm avg}=q[\sum\phi_x\phi_xI^R_{3l}+\sum\phi_y\phi_yI^L_{3l}]_{\rho=R}$ [S4 Eq. 65, 72].
Kausel–Peek (1982, BSSA) is the journal version of [S4].

### 3.3 SASSI central zone: axisymmetric FE core + cylindrical consistent boundary  [S3, S12, S18] + **[R] [V]**

The source excerpts confirm the architecture: an FE core of radius $R_0$ ("Radius of Central Zone") plus
the "Kausel–Waas axisymmetric consistent boundaries" [S12]. They do not give the internal
discretization of the core. **Reconstruction:** use **one element radially** (nodes on the axis and at
$\rho=R_0$), with the same vertical discretization as SITE (one element per sublayer). This is the only
discretization that the single input $R_0$ fully defines. It also explains the manual's rules
(§3.4): every off-diagonal pair with $r\ge h>R_0$ is then evaluated in the exact far field.

**(i) Core element** (harmonic $\mu$, per radian). Nodal unknowns are $\tilde u_\rho,\tilde u_\theta,\tilde u_z$.
Interpolation is bilinear: $L_0=1-\rho/R_0$, $L_1=\rho/R_0$ in $\rho$, and TLM-linear $N_a(z)$ in $z$.
The strain vector is

$$
\tilde\varepsilon=\Big[\partial_\rho\tilde u_\rho,\ \tfrac{\tilde u_\rho-\mu\tilde u_\theta}{\rho},\ \partial_z\tilde u_z,\ \partial_z\tilde u_\rho+\partial_\rho\tilde u_z,\ \tfrac{\mu\tilde u_\rho}{\rho}+\partial_\rho\tilde u_\theta-\tfrac{\tilde u_\theta}{\rho},\ \partial_z\tilde u_\theta+\tfrac{\mu\tilde u_z}{\rho}\Big]^T,
$$

with $D$ = isotropic $(\lambda^*,G^*)$. Then $K_c=\iint B^TDB\,\rho\,d\rho\,dz$ and
$M_c=\iint\rho_mN^TN\rho\,d\rho\,dz$, integrated with 3×3 Gauss points so that $\rho=0$ is never
evaluated. Use the same mass lumping in $z$ as in SITE. **Axis constraints:** for $\mu=0$,
$\tilde u_\rho=\tilde u_\theta=0$; for $\mu=1$, $\tilde u_\rho=\tilde u_\theta$ (the same unknown) and
$\tilde u_z=0$. For a dashpot base, add $i\omega c\int_0^{R_0}L_aL_b\rho\,d\rho$ on the bottom-interface
dofs.

**(ii) Exterior field.** For $\rho\ge R_0$ it is a superposition of outgoing modes, written
$\tilde u(\rho)=\Psi_\mu(\rho)\,\alpha$ with $\alpha=[a_{1..2N_f};b_{1..N_f}]$:

$$
\Psi_\mu(\rho)=\begin{bmatrix}
\Phi_x\,{\rm diag}\,H_\mu'(k^R\rho) & \Phi_y\,{\rm diag}\,\frac{\mu H_\mu(k^L\rho)}{k^L\rho}\\[2pt]
\Phi_x\,{\rm diag}\,\frac{\mu H_\mu(k^R\rho)}{k^R\rho} & \Phi_y\,{\rm diag}\,H_\mu'(k^L\rho)\\[2pt]
-\Phi_z\,{\rm diag}\,H_\mu(k^R\rho) & 0\end{bmatrix},\qquad H_\mu' = \frac{dH_\mu(x)}{dx}.
$$

These are the columns of Kausel's $C_\mu$ operator [S4 Eq. 7] with $J\to H^{(2)}$. Derivatives use
$H_\mu''(x)=-H_\mu'(x)/x-(1-\mu^2/x^2)H_\mu(x)$, giving
$\partial_\rho[H_\mu'(k\rho)]=kH_\mu''$ and
$\partial_\rho[\mu H_\mu(k\rho)/(k\rho)]=\mu k\,[H_\mu'/x-H_\mu/x^2]$.

**(iii) Boundary stiffness.** The consistent generalized forces (per radian) exerted on the core at
$\rho=R_0$ by mode columns are

$$
\tilde f_\rho=R_0\big[E_{\lambda+2G}\,\partial_\rho\tilde u_\rho+E_\lambda(\tilde u_\rho-\mu\tilde u_\theta)/R_0+Q_\lambda\tilde u_z\big],\quad
\tilde f_\theta=R_0\big[E_G(\mu\tilde u_\rho/R_0+\partial_\rho\tilde u_\theta-\tilde u_\theta/R_0)\big],\quad
\tilde f_z=R_0\big[Q_G\tilde u_\rho+E_G\,\partial_\rho\tilde u_z\big],
$$

with $E_c$ and $Q_c$ as in §2.6. Collect them as $\mathcal T_\mu$ (one column per mode) and set
$\boxed{R_\mu=-\mathcal T_\mu\,\Psi_\mu(R_0)^{-1}}$. This is a full $3N_f\times3N_f$ matrix. For
$\mu=0$ it splits into $(\rho,z)$ Rayleigh and $\theta$ (torsion, Love) blocks.

**(iv) Solve.** For each load interface $n$:

$$
\big[K_c-\omega^2M_c+R_\mu\big]\tilde U=\tilde F,\qquad \tilde F=\frac{P}{\pi}\ (\mu=1,\ \text{on the tied axis dof}),\quad \tilde F=\frac{P}{2\pi}\ (\mu=0).
$$

The factors $\pi$ and $2\pi$ are the angular integrals of $\cos^2$, $\sin^2$ and 1. Then
$\alpha=\Psi_\mu(R_0)^{-1}\tilde U_{\rm boundary}$ and $\tilde u(\rho)=\Psi_\mu(\rho)\alpha$ for all
$\rho\ge R_0$ and all interfaces.

**(v) Verification [V]** (`point3.py`; 28 layers, f = 5 Hz, $R_0$ = 0.9 m). The far field was compared
with the exact Kausel point-load functions (§3.2) for both $\mu$ and several load and observation
depths:

| Distance | Agreement |
|---|---|
| $\rho=2R_0$ | within 1–6 % (expected: the core is a distributed load, not a point load) |
| $3.3R_0$ | ≤2.3 % |
| $6.7R_0$ | ≤0.8 %, except two components (1.9 % and 4.5 %) whose magnitude is 10–100 times below the dominant component |
| $13$–$28R_0$ | ≤0.25 % |

The assembled flexibility is reciprocal to $\sim10^{-4}$ (`test8.py`).

**(vi) Self-flexibility (diagonal).** This is the core **axis** displacement at the loaded interface.
[V] (`test6.py`): it is much softer than a uniform disk of radius $R_0$. The equivalent uniform-disk
radius is only about $0.25$–$0.29R_0$. A "tributary-area disk" implementation would therefore **not**
reproduce SASSI's diagonal terms. Implement the FE core.

### 3.4 Central-zone radius rules **[S18]** (figures on manual p.128; see `05a` §8.3)

| Mesh | Rule |
|---|---|
| 3-D square excavation mesh, element size $h$ | $R_0=0.90h$ |
| 3-D triangular (circular foundation) mesh | $R_0=0.85h$ |
| 2-D (plane strain) | $R_0=h$ |

* For non-uniform meshes, use the average element radius and run sensitivity studies (max/avg/min)
  [S18 1.5.1 item 9].
* [R] Rationale: $R_0<h$ guarantees that all neighbour pairs lie in the exact far field. $R_0$ then only
  calibrates the near-singular diagonal terms to the nodal flexibility of a mesh of size $h$.

### 3.5 POINT2 (plane strain) [S5 Fig. 2-2, S18] + [R][V]

The core is a strip $|x|\le R_0$ (2 elements, nodes at $-R_0, 0, R_0$) with Waas boundaries (§2.6) at
$\pm R_0$. Unit line loads are applied at the centre node. The far field is
$U(x)=\sum\alpha_j\phi_je^{-ik_j(|x|-R_0)}$, with the vertical component odd in $x$ for horizontal loads.
The exact reference is Kausel's line-load Green function ([S4] Eq. 47, 48; [V] `test7.py`):

$$
\begin{Bmatrix}U_x\\U_z\end{Bmatrix}=\frac1{2i}\begin{bmatrix}\Phi_xE^R_{|x|}K_R^{-1}\Phi_x^T&\pm\Phi_xE^R_{|x|}K_R^{-1}\Phi_z^T\\ \pm\Phi_zK_R^{-1}E^R_{|x|}\Phi_x^T&\Phi_zE^R_{|x|}K_R^{-1}\Phi_z^T\end{bmatrix},\qquad U_y=\frac1{2i}\Phi_yE^L_{|x|}K_L^{-1}\Phi_y^T .
$$

Here $E_{|x|}={\rm diag}(e^{-ik_j|x|})$, and the upper sign applies for $x\ge0$. These are i-scaled
quantities: the true vertical displacement due to a horizontal load is $\mp\frac12\Phi_zK^{-1}E\Phi_x^T$
[S4 p.23].

### 3.6 FILE3 content and how ANALYS evaluates it [S18] + [R]

* Confirmed: FILE3 holds the "point load solution". Its size grows with the number of embedment layers
  $L$ [S18].
* [R] **Store, per frequency, per load interface $n=1..L+1$, and for $\mu\in\{0,1\}$:** the mode
  amplitudes $\alpha_{n,\mu}$ ($3N_f$ complex), the axis displacements at interfaces $1..L+1$, and
  $R_0$. ANALYS then evaluates $\tilde u(\rho)=\Psi_\mu(\rho)\alpha$ **exactly at any distance**. No
  tabulation in $r$ is needed. If tabulation is wanted for speed, store $\tilde u(\rho_k)$ on a
  log-spaced grid and interpolate $\log|\cdot|$ and the unwrapped phase in $\log\rho$.

---

## 4. ANALYS: impedance and the flexible-volume equation  — task items (2) and (3)

### 4.1 Building $F_{ff}$ from POINT [S3, S5, S18] + [R][V]

Interaction node $i$ lies at $(x_i,y_i)$ on interface $m$; node $j$ lies at $(x_j,y_j)$ on interface $n$.
Let $r=|\mathbf x_i-\mathbf x_j|$ and $\theta={\rm atan2}(y_i-y_j,x_i-x_j)$, measured from load point $j$
to observation point $i$. From POINT (load at $n$, observe at $m$), take $\tilde u,\tilde v,\tilde w$
= $(\tilde u_\rho,\tilde u_\theta,\tilde u_z)$ for $\mu=1$, and $\tilde p,\tilde q$ =
$(u_\rho,u_z)$ for $\mu=0$. Rotating the $\cos\theta$ / $-\sin\theta$ pattern gives the 3×3 block
(rows = $u_{x,y,z}$ at $i$, columns = unit $P_{x,y,z}$ at $j$):

$$
F_{ij}=\begin{bmatrix}\tilde u c^2+\tilde v s^2&(\tilde u-\tilde v)sc&\tilde p\,c\\ (\tilde u-\tilde v)sc&\tilde u s^2+\tilde v c^2&\tilde p\,s\\ \tilde w\,c&\tilde w\,s&\tilde q\end{bmatrix},\qquad c=\cos\theta,\ s=\sin\theta .
$$

* **Same vertical line** ($r=0$): $F_{ij}={\rm diag}(\tilde u_{\rm ax},\tilde u_{\rm ax},\tilde q_{\rm ax})$
  from the core axis values. The coupling terms vanish by symmetry.
* $0<r<R_0$ (mesh finer than the rule): interpolate linearly in $r$ between the axis value and
  $\tilde u(R_0)$, consistent with the core shape functions.
* Reciprocity: $F_{ij}=F_{ji}^T$ analytically for the exact Green's functions. With the FE core it holds
  to about $10^{-4}$ [V]. [R]: symmetrize $F\leftarrow\frac12(F+F^T)$, because SASSI treats $X_f$ as
  symmetric [S3].
* 2-D: $F_{ij}$ is the 2×2 (P-SV) or 1×1 (SH) block from §3.5, with the odd terms multiplied by
  ${\rm sgn}(x_i-x_j)$.

### 4.2 Impedance [S3, S18 Eq. 4.3, S12]

$$
X_{ff}(\omega)=F_{ff}(\omega)^{-1}.
$$

* The flexibility is full, symmetric and complex [S3: "An efficient in-place inversion subroutine is
  used"]. ACS writes it as $K+iD=(f+ig)^{-1}$ [S18 Eq. 4.3].
* Its size is $3n_f\times3n_f$ for $n_f$ interaction nodes [S12]. [R]: use a complex-symmetric
  $LDL^T$ factorization (LAPACK `zsytrf/zsytri`) or LU.

### 4.3 Excavated-soil dynamic stiffness [S3, S5, S12, S18]

* The excavated soil is an FE model (HOUSE) of the soil volume that the basement displaces. It uses the
  **free-field layer properties** of each embedment layer (the same complex moduli and densities as
  SITE) [S3: "mass and stiffness of the structure are reduced by the corresponding properties of the
  volume of soil excavated"; S18].
* Its dynamic stiffness is $C^e(\omega)=K^{e*}-\omega^2M^e$. The mass is ½ lumped + ½ consistent [S5,
  S18].
* [R] Element sizes in plan should match $R_0$ through the rule of §3.4. Element layers must coincide
  with the SITE interfaces.

### 4.4 The flexible-volume equation [S3 Eq. 11, S5 Eq. 2-5, S12 slide 9, S18 Eq. 2.1]

Use the DOF sets of the manual: $s$ = structure only, $i$ = nodes shared by structure and excavated soil,
$w$ = excavated-soil nodes not attached to the structure, and $f$ = interaction DOFs. Superscript $s$ =
structure, $e$ = excavated soil.

$$
\begin{bmatrix}
C^s_{ss}&C^s_{si}&0\\
C^s_{is}&C^s_{ii}-C^e_{ii}+X_{ii}&-C^e_{iw}+X_{iw}\\
0&-C^e_{wi}+X_{wi}&-C^e_{ww}+X_{ww}
\end{bmatrix}
\begin{Bmatrix}U_s\\U_i\\U_w\end{Bmatrix}=
\begin{Bmatrix}F_s\\X_{ii}U'_i+X_{iw}U'_w+F_i\\X_{wi}U'_i+X_{ww}U'_w\end{Bmatrix}.
$$

* $U$ are **total (absolute)** complex displacement amplitudes. $U'$ is the free-field motion from SITE
  (FILE1) at the interaction nodes, normalized to the unit control motion. $X_{ab}$ is the sub-block of
  $X_{ff}$, taken as zero if $a$ or $b$ is not in $f$. $F$ are external forces from FORCE (zero for
  seismic runs; $U'=0$ for pure vibration runs).
* When every node is an interaction node (FV), this reduces to [S3] Eq. 11:
  $\begin{bmatrix}C_{ss}&C_{si}\\C_{is}&C_{ii}-C_{ff}+X_{ff}\end{bmatrix}\begin{Bmatrix}U_s\\U_i\end{Bmatrix}=\begin{Bmatrix}0\\X_{ff}U'_f\end{Bmatrix}$.
  The term $C_{ii}-C_{ff}$ "simply indicates … that the stiffness and mass of the excavated soil are
  subtracted from the stiffness of the structure" [S3].
* **Built-in check [S18 Ch.2; V-plan].** If $C^s=C^e$ on the excavated volume, the equation gives
  $X(U-U')=0$, so $U=U'$ (no SSI). This is the first verification test (§9.2).
* **General assembly [R].** $\mathbf C=\mathcal A_s^TC^s\mathcal A_s-\mathcal A_e^TC^e\mathcal A_e+\mathcal A_f^TX_{ff}\mathcal A_f$
  and $\mathbf b=\mathcal A_f^TX_{ff}U'_f+\mathbf F$, with Boolean gather matrices $\mathcal A$.

### 4.5 Method variants (choice of the set $f$) [S12, S13–S16, S18]

| ACS name | Other names | Interaction set $f$ | Notes |
|---|---|---|---|
| **FV** | Direct, flexible volume | all excavated-soil nodes, $i\cup w$ | Reference solution. $X$ is full on all excavation nodes [S12 slide 9: "All Excavated Soil nodes are interaction nodes (include exact equations of motion)"]. |
| **FI-FSIN** | Subtraction (SM) | soil–foundation interface nodes only (outer perimeter of the basement) [S14, S15] | The $w$ rows keep only $-C^e_{ww}$, so the system is near-singular where $\det C^e_{ww}(\omega)\approx0$: at the **natural frequencies of the excavated volume with the interface nodes fixed**. This produces the documented spurious peaks "at and above the natural frequency of excavated soil volume" [S13, S14]. It is a method limitation, "not a discretization (λ/5) issue, not a programming error" [S13]. |
| **FI-EVBN** | Modified subtraction (MSM) | FSIN + **all excavation-boundary nodes, including the ground-surface face** | The surface constraint raises the excavated-volume frequencies [S16 Table 1]. It is much more accurate than SM, but anomalies can still appear above the new excavated-volume frequency [S13]. |
| **ESM** | Extended subtraction (SASSI2010) | perimeter + **additional horizontal planes** of interior nodes [S15] | Validated against ANSYS full harmonic analysis [S15]. |
| **FFV** | Fast flexible volume (Ghiocel 2013) | EVBN + **internal node layers at a regular vertical skip** ("FFV-Skip2/5") [S12] | Example [S12]: SMR model with 7 936 nodes for FV, 4 016 for Skip2 (20 % runtime) and 2 252 for MSM (6 %). Skip2 was "highly accurate". |

DOE screening frequencies for SM [S16]: soil-layer frequency $f_{SL}=V_s/4H$ ($H$ = embedment depth) and
the excavated-volume frequency $f_{EV}$ (FE model fixed at the interaction nodes). The SM is suspect when
input energy or structural frequencies lie at or above these.

### 4.6 Solving efficiently [R]

Order the unknowns as $n$ (non-interaction) and $f$. Factor the sparse $\mathbf C_{nn}$ once per
frequency and form the dense Schur complement
$S=\mathbf C_{ff}-\mathbf C_{fn}\mathbf C_{nn}^{-1}\mathbf C_{nf}$. Then solve
$(S+X_{ff})U_f=X_{ff}U'_f+\mathbf F_f-\mathbf C_{fn}\mathbf C_{nn}^{-1}\mathbf F_n$ and back-substitute.
Cost is dominated by $O(n_f^3)$ for $X_{ff}$; this matches the scaling noted in [S18].

### 4.7 Output [S18, S19]

For each frequency and load case, ANALYS stores the complex **transfer functions**
$H_k(\omega)=U_k/U_{cp}$ at all output DOFs in FILE8. "For typical problems, the transfer functions only
need to be solved for 60 to 80 frequencies" [S19]. The rest are obtained by interpolation (§5).

---

## 5. MOTION: frequency interpolation and convolution  — task item (5)

### 5.1 The SASSI (Tajirian) interpolation formula [S6, S7] [V]

The method was developed by Tajirian (1981) [S7] and is used in SASSI [S6]. It assumes that, locally,
every transfer function has the analytic form of a **two-degree-of-freedom system with hysteretic
damping**. That form is a ratio of two quartic polynomials in $\omega$ with only even powers and complex
coefficients ([S6] Eq. 1; [S7] Eq. 12):

$$
U(\omega)=\frac{C_1\omega^4+C_2\omega^2+C_3}{\omega^4+C_4\omega^2+C_5}.
$$

The five complex constants come from the five computed values $U_p=U(\omega_p)$, $p=1..5$, in a window
([S6] Eq. 2, written here with the sign derived from Eq. 1):

$$
\begin{bmatrix}\omega_p^4&\omega_p^2&1&-\omega_p^2U_p&-U_p\end{bmatrix}_{p=1..5}\begin{Bmatrix}C_1\\C_2\\C_3\\C_4\\C_5\end{Bmatrix}=\begin{Bmatrix}\omega_p^4U_p\end{Bmatrix}_{p=1..5}.
$$

**Note:** the MHI transcription [S6] prints $+U_p$ in the fifth column. Moving the denominator of Eq. 1
across gives $-U_p$. The minus sign is verified: the fitted form reproduces an exact 2-DOF hysteretic
transfer function to $10^{-14}$ over 0.5–12 Hz from five samples (`test5.py`).

* Interpolated values equal the computed values at the five frequencies [S6].
* Rule of thumb [S6, S7]: a 5-point window should contain no more than two resonant peaks; close peaks
  need denser frequencies. TLUSH: "normally 5 frequency points are needed in the vicinity of two single
  frequency peaks". Use 30–40 points per curve for dams [S7] and 60–80 for SSI [S19].
* **Rank deficiency [V][R].** If the data in a window come from fewer than two modes (e.g. an SDOF-like
  or nearly flat transfer function), the 5×5 system is singular (rank 4). The pole–zero pair is then
  not determined. In that case **every** solution of the family interpolates exactly. Solve by SVD or
  minimum-norm least squares (`lstsq`, rcond ≈ 1e-10). [V] (`test5b_interp_rank.py`): a hysteretic SDOF
  gives rank 4 and $6\times10^{-15}$ error; a viscous SDOF gives rank 5 and $6\times10^{-5}$ error, so
  it is not exact. Scale $\omega$ by the window's largest frequency for conditioning.
* **Correction to `05c` §B.6.4 and test C8.** `05c` assumed $(β_1+β_2x+β_3x^2)/(1+β_4x+β_5x^2)$ with a
  shifted linear variable $x$. The sourced form uses **$\omega^2$**, which has no odd powers. C8 should
  use a hysteretic SDOF/2DOF, $k^*/(k^*-m\omega^2)$, rather than a viscous SDOF.

### 5.2 Windowing (interpolation options) [S6, S18] + [R]

* Confirmed: the frequency range "is subdivided into smaller regions each of which contains the transfer
  function solution for 5 frequencies … For the last region, the solution from the previous region can be
  augmented, if necessary, to form the solution of 5 frequencies" [S6]. ACS options 0–5 use
  non-overlapping or overlapping, moving-average 5-point windows, with or without point shifts.
  **Option 0 = overlapping, moving average** [S6]. Option names and descriptions are in [S18] (see `05c`
  §B.5.5). Option 6 (spline) is an ACS addition.
* [R] Algorithms:
  - **Non-overlapping (option 1, "original SASSI 1982")**: windows $\{1..5\},\{5..9\},\dots$ share end
    points. The last window is the last 5 points.
  - **Moving average (option 0/2)**: for a target $\omega\in[\omega_j,\omega_{j+1}]$, evaluate every
    5-point window of consecutive computed points that contains $[\omega_j,\omega_{j+1}]$ (start index
    $m\in[\max(1,j-3),\min(j,N-4)]$). Average them, with weights decreasing with the distance of
    $\omega$ from the window centre for option 0 ("weighted"). Plain averaging gives option 2.
  - **Shifted variants (4, 5)**: non-overlapping windows whose starts are offset by 1 or 2 points;
    average with the unshifted solution.
  - **Guard**: if $\min|\omega^4+C_4\omega^2+C_5|$ over the window is $<10^{-3}\max|\cdot|$ (a spurious
    near-real pole inside the window), use a complex cubic through the 4 nearest points in that window.
* **Historical fallback [S7]:** LUSH and FLUSH (Lysmer et al. 1974, 1975) used "linear interpolation on
  the inverse of the amplification functions". This is a robust option when $N<5$.

### 5.3 Convolution [S18, S19, S7]

1. FFT the control acceleration $a_{cp}(t)$ ($N$ points, $\Delta t$), giving $A(\omega_k)$ with
   $\omega_k=2\pi k\Delta f$ and $\Delta f=1/(N\Delta t)$.
2. SASSI analysis frequencies are integer *frequency numbers* times $\Delta f$ [S18 SITE options].
3. Interpolate $H$ onto every $\omega_k\le\omega_{\rm cut}$ and set $H=0$ above the cut-off.
4. [R] Use $H(0)=1$ for the output component parallel to the input and $0$ otherwise: the system moves
   rigidly with the free field as $\omega\to0$. The manual's check that the frequency-1 ATF is about 1
   supports this [S18 4.1.2].
5. $a_k(t)={\rm IFFT}[H_k(\omega)A(\omega)]$, with a trailing quiet zone of zeros so that free vibration
   decays [S7, S18].
6. Response spectra follow.

STRESS uses the same interpolation on element strain and stress transfer functions [S18].

---

## 6. SOIL: SHAKE equivalent-linear site response  — task item (6)

The ACS SOIL module "uses the SHAKE methodology for simulating the vertically propagating wave within
horizontal soil layering" and the Seed–Idriss equivalent-linear model [S18]. The equations below are from
the SHAKE report [S8] and the public-domain SHAKE91 code [S9]. [V]: an independent implementation in
`test4.py` reproduces them.

### 6.1 Wave solution and recursion [S8 Eq. 8–9; S9 `CXSOIL`, `AMP`] (also BNL [S17] Eq. 1–11)

In layer $m$ (local depth $z_m$ downward from its top):

$$
u_m=E_me^{i(\omega t+k_m^*z_m)}+F_me^{i(\omega t-k_m^*z_m)},\qquad k_m^*=\omega/V_m^*,\quad V_m^*=\sqrt{G_m^*/\rho_m}.
$$

$E$ is the **upgoing (incident)** wave and $F$ the **downgoing (reflected)** wave. The stress is
$\tau=G^*\partial u/\partial z$. Continuity of $u$ and $\tau$ gives

$$
\begin{aligned}
E_{m+1}&=\tfrac12E_m(1+\alpha_m)e^{ik_m^*h_m}+\tfrac12F_m(1-\alpha_m)e^{-ik_m^*h_m},\\
F_{m+1}&=\tfrac12E_m(1-\alpha_m)e^{ik_m^*h_m}+\tfrac12F_m(1+\alpha_m)e^{-ik_m^*h_m},
\end{aligned}\qquad
\alpha_m=\frac{\rho_mV_m^*}{\rho_{m+1}V_{m+1}^*}=\sqrt{\frac{\rho_mG_m^*}{\rho_{m+1}G_{m+1}^*}} .
$$

* The free surface requires $E_1=F_1$; SHAKE starts with $E_1=F_1=1$ and recurses down to the halfspace.
* Within motion at the top of layer $m$ is $E_m+F_m$; outcrop motion is $2E_m$ (`AMP`: `AA = E+FF`, or
  `2*E` when outcropping).
* The transfer function between the object layer $n$ and layer $m$ is
  $(E_m+F_m)/(E_n+F_n)$, or $/(2E_n)$ for an outcrop object motion [S8, S9].
* The halfspace is an elastic or viscoelastic rock with its own $G^*$ and $\rho$, entering through
  $\alpha_N$.

### 6.2 Strain and effective strain [S9 `STRT`; S8]

* Strain is computed at **mid-height of each sublayer** from the acceleration wave amplitudes:
  $\gamma_m(\omega)=\big(E_m^ae^{ik^*h/2}-F_m^ae^{-ik^*h/2}\big)/(i\omega V_m^*)$. This is the code's
  `X = (E*EX − F/EX)·g/(i2πV)/f`. Then inverse FFT and take $\gamma_{\max}=\max_t|\gamma(t)|$.
* Effective strain is $\gamma_{\rm eff}=R_\gamma\,\gamma_{\max}$, with **$R_\gamma=(M-1)/10$**, where $M$
  is the earthquake magnitude: "for M = 5, the ratio would be 0.4, for M = 7.5, the ratio would be
  0.65" [S9 manual, Option 5]. The 1972 manual says "between 0.5 and 0.7 … 0.55 to 0.65 is usually
  adequate" [S8].

### 6.3 Property update and convergence [S9 `STRT`; S18] + [R]

* New $G=G_{\max}\cdot(G/G_{\max})(\gamma_{\rm eff})$ and $\beta=\beta(\gamma_{\rm eff})$. The curves are
  interpolated **linearly in $\log_{10}\gamma$** (code: `GN = AS*log10(SS)+BS`).
* Errors are reported as $\Delta G=(G_{\rm new}-G_{\rm used})/G_{\rm new}\times100\%$, and likewise for
  β. SHAKE91 runs a **fixed number of iterations** chosen by the user. Its tolerance test
  `IF (DGMAX.LT.ERR)` is commented out in the code, and the last iteration's errors are printed. The ACS
  manual recommends 8 iterations [S18 SOIL options].
* [R] Stop when $\max(|\Delta G|,|\Delta\beta|)<1$–$2\%$, or at the iteration limit.
* Vertical input (P waves) uses $V_p$ and no iteration [S18].
* **Hand-off to SSI [S18][R]:** the strain-compatible $V_s=\sqrt{G/\rho}$ and $\beta_s$ go to SITE and to
  the HOUSE excavated-soil layers (FILE88). For $V_p$ and $\beta_p$, either keep $\nu$ constant
  (unsaturated soil) or keep $V_p$ (saturated soil, $V_p\gtrsim1500$ m/s), with $\beta_p=\beta_s$
  capped by a user limit. This choice is not stated in the sources.
* Use one complex-modulus convention (§1) in both SOIL and SITE. Otherwise the SSI free field will not
  reproduce the SOIL motion at the control point.

---

## 7. Per-frequency SASSI data flow (summary) [S5, S18] + [R]

```
SOIL (optional):  SHAKE iteration -> strain-compatible Vs, beta  -> SITE/HOUSE layer properties
for each analysis frequency f (frequency number x df):
  SITE Mode 1 : build layers (+ 1.5*lambda_s halfspace sublayers + dashpots) ; solve Rayleigh/Love eigenproblems -> FILE2
  SITE Mode 2 : unit-normalized free field at all interfaces (wave mix, control point)       -> FILE1
  POINT       : for n = 1..L+1, mu = 0,1 : FE core (R0) + cylindrical TB -> alpha, axis values  -> FILE3
  ANALYS      : F_ff (3x3 blocks, cos/sin theta) -> symmetrize -> X_ff = F_ff^-1
                assemble C_s - C_e + X_ff ; RHS = X_ff U'_f (+F) ; solve -> transfer functions -> FILE8
MOTION/STRESS : 5-point 2DOF interpolation in omega^2 -> convolution with control-motion FFT -> IFFT
```

---

## 8. Confirmed vs reconstructed: summary table

| Item | Confirmed by source | Reconstructed here |
|---|---|---|
| TLM layer matrices A, B, G, M | **Exact** (Kausel Table 1 [S4]). Waas form with $iB$ [S5]. | Mixed-mass default for SITE layers |
| Eigenproblem, root choice, normalization | [S4] Eq. 14–22 and Appendix | — |
| Halfspace: 1.5λ, n ≤ 20 sublayers increasing with depth, dashpots ρV | [S18, S5] | Thickness series, cap λ/8, complex $V^*$ in the dashpot |
| Vertical body-wave free field; control point at top of layer | [S5, S18] | Dashpot load $2i\omega cE$, outcrop handling |
| Inclined body waves and surface-wave field | Form only [S5] | $k=\omega\sin\theta/V$, base load, mode-selection definitions |
| Point-load Green's functions (Hankel series) | **Exact** [S4] Eq. 88–89 | — |
| POINT3 = axisymmetric FE core + axisymmetric TB, cylindrical-to-Cartesian transform | [S3, S12, S18] | One radial element, full TB formula, storage of α ([V] against [S4]) |
| $R_0$ rules 0.9h / 0.85h / h | [S18] | Rationale |
| $X=F^{-1}$; FV, SM, MSM, ESM, FFV node sets | [S3, S5, S12–S16, S18] | Schur-complement solver, symmetrization |
| Interpolation formula (2DOF, ω², 5 points), windows | [S6, S7] | Window algorithms, rank-deficiency handling |
| Damping $G(1-2\beta^2+2i\beta\sqrt{1-\beta^2})$ | SHAKE91, TLUSH, Udaka–Lysmer 1973 [S7–S9]; SASSI inferred [S10+V] | Application to $V_p$/$\beta_p$ and HOUSE |
| SHAKE recursion, strain, $(M-1)/10$, iteration | [S8, S9] | Convergence tolerance, $V_p$ policy |

---

## 9. Verification

### 9.1 Checks run for this note (`docs/spec/R1_checks/`, `python run_all.py`)

| # | Script | Check | Result |
|---|---|---|---|
| V1 | `tlm.py` | Modal flexibility (Kausel Eq. 49–52) vs direct inverse of $Ak^2+Bk+C$, 10 layers, damped, real and complex k | rel. error $2\times10^{-14}$ |
| V2 | `test2.py` | Fundamental Rayleigh speed, uniform stratum, ν = 0.25, $h=\lambda/20$ | 0.9210 / 0.9194 / 0.9202 (consistent / lumped / mixed); exact 0.9194 |
| V3 | `test3.py` | Kausel point-load Green's functions vs Boussinesq/Cerruti (quasi-static, 115 graded layers, 400 m deep) | $u_\rho$, $u_\theta$, $u_z$ for both loads within 0.1–1.7 % at ρ = 1–5 m; signs correct |
| V4 | `test4*.py` | 1-D SH: TLM + variable-depth halfspace + dashpot vs exact SHAKE recursion | within TF ≤5.6 % at $h=\lambda/20$ (20 % at $\lambda/10$), independent of $n$; outcrop TF ≤0.7 % (uniform, n = 20), ≤1.5 % (n = 10); n = 5 up to 20 % |
| V5 | `test5.py`, `test5b` | 5-point ω² rational interpolation | exact ($10^{-14}$) for a 2-DOF hysteretic system; rank 4 but exact for a hysteretic SDOF; $6\times10^{-5}$ for a viscous SDOF |
| V6 | `test7.py` | Kausel line-load Green function (Eq. 47) vs numerical inverse Fourier transform | agreement to 4 digits |
| V7 | `point2.py` | 2-D FE strip + Waas TB vs Eq. 47 | ≤0.3 % for $|x|\ge2R_0$ |
| V8 | `point3.py` | 3-D FE core (one radial element) + cylindrical TB vs Kausel point load | 1–6 % at $2R_0$, ≤2.3 % at $3.3R_0$, ≤0.8 % at $6.7R_0$ on the dominant components, ≤0.25 % beyond $13R_0$ |
| V9 | `test6.py` | FE-core self-flexibility vs uniform-disk average | equivalent disk radius 0.25–0.29 $R_0$, so a tributary disk is not equivalent |
| V10 | `test8.py` | Reciprocity of assembled 3×3 blocks | $F_{ij}-F_{ji}^T$ ≈ 1e-4 to 8e-4 relative |
| V11 | `test9_damping.py` | SDOF 1/(2Hmax) vs Ostadan et al. 2004 SASSI table | matches the Udaka–Lysmer form (19.6 % vs 19.5 % at β = 20 %); $1+2i\beta$ gives 18.6 % |

### 9.2 Recommended acceptance tests for the implementation [R]

1. **No-SSI identity.** Use a structure with the same properties as the excavated soil (FV) and a
   vertical SV input. Expected: $U=U'$ at every node to round-off, and the surface ATF equals the SITE
   free-field TF.
2. **1-D column vs SOIL/SHAKE.** Use a free field without a structure. Expected: the TFs match SHAKE
   within V4 tolerances, and the control-point TF is 1.
3. **Rigid massless circular disk compliances** on a layer over a rigid base (Tajirian & Tabatabaie
   1985 Case 1 [S3]: $V_s=1$, $V_p=2$, $\rho=1$, $H=3$, β = 5 % and 15 %, compared with LUCON).
   Static limits use the Kausel/Gazetas layer formulas:
   $K_h\approx\frac{8GR}{2-\nu}(1+\frac R{2H})$ and $K_v\approx\frac{4GR}{1-\nu}(1+1.28\frac RH)$.
   Expected: vertical resonances at $f=(2n-1)V_p/4H$ and horizontal at $(2n-1)V_s/4H$ [S3].
4. **Halfspace disk** (Luco/Veletsos impedance tables), using 20 halfspace sublayers.
5. **FV vs SM vs MSM on an embedded box.** Reproduce the SM anomaly at $f_{EV}$ and its reduction with
   MSM and FFV [S13, S16].
6. **ANSYS full harmonic benchmark** (Anderson & Ostadan 2014 approach [S15]): a large soil island with
   a fixed base and rollers or absorbing sides, with damping mapped per §1.2. Compare ATFs and ISRS.
7. **Interpolation.** An exact 2-DOF hysteretic TF sampled at 5 points must be reproduced to round-off.
   ISRS from interpolated and dense TFs should agree within 1–2 %.

---

## 10. Open questions (not resolvable from the sources accessed)

1. The exact SASSI2000 sublayer-thickness law for the variable-depth halfspace. A reconstruction is
   given in §2.4; the original is in UCB/GT/81-02 and the SASSI2000 theory manual (NISEE, paywalled).
2. The internal discretization of Tajirian's "special cylindrical axisymmetric finite elements" (one
   radial element and linear interpolation are assumed), and whether their mass is mixed or consistent.
3. Whether SITE's surface-wave mode selection ("shortest wavelength" / "least decay") uses Re k or
   Re k with Im k as assumed, and how inclined body waves treat the halfspace (exact stiffness or the
   dashpot base).
4. The exact weighting in interpolation option 0 ("weighted averaging") and the point-shift rules of
   options 4 and 5.
5. Whether POINT stores mode amplitudes (assumed) or tabulated $u(r)$, and how ANALYS treats
   $0<r<R_0$.
6. $V_p$ and $\beta_p$ update policy after SOIL iterations.
7. Whether HOUSE structural elements use the same complex-modulus factor as soil. This is assumed; it is
   consistent with [S7] and [S18] item 13.

---

## 11. References (URLs as accessed on 2026-10-01)

* **[S1]** Lysmer, J., Tabatabaie-Raissi, M., Tajirian, F., Vahdani, S., Ostadan, F. (1981). *SASSI – A System for Analysis of Soil-Structure Interaction*, Rep. UCB/GT/81-02, UC Berkeley, April 1981. NISEE e-library, membership required: https://nisee.berkeley.edu/elibrary/Text/200701181 — *not accessed*.
* **[S2]** Lysmer, J., Ostadan, F., Chin, C.C. (1999). *SASSI2000 – Theoretical Manual / User's Manual*, UC Berkeley. NISEE: https://nisee.berkeley.edu/elibrary/Text/200703021 — *not accessed*.
* **[S3]** Tajirian, F.F., Tabatabaie, M. (1985). "Vibration analysis of foundations on layered media," ASCE Proc. *Vibrations for Soils and Foundations*, Detroit, Oct. 1985. http://mtrassoc.com/wp-content/uploads/2013/09/Vibration-Analysis-of-Foundations-on-Layered-Media.pdf
* **[S4]** Kausel, E. (1981). *An Explicit Solution for the Green Functions for Dynamic Loads in Layered Media*, MIT Research Report R81-13 (NSF/CEE-81060). https://nehrpsearch.nist.gov/static/files/NSF/PB82147893.pdf — Table 1, Eq. 14–22, 47–52, 61–72, 88–89, Table 2 and Appendix used. Journal version: Kausel, E., Peek, R. (1982), *BSSA* 72(5):1459–1481.
* **[S5]** Hsiung, S.-M., Chowdhury, A.H. (2011). *Seismic Effects on Soil-Structure Interactions*, CNWRA report for the US NRC (Contract NRC-02-07-006). It transcribes SASSI2000 theory (Eq. 2-1 to 2-15). https://www.nrc.gov/docs/ml1125/ml112580204.pdf
* **[S6]** Mitsubishi Heavy Industries (2013). Response to RAI 810-5874 Q03.07.02-91 (US-APWR), ACS SASSI interpolation equations (1)–(2). https://www.nrc.gov/docs/ML1305/ML13057A398.pdf
* **[S7]** Kagawa, T., Mejia, L.H., Seed, H.B., Lysmer, J. (1981). *TLUSH: A Computer Program for the Three-Dimensional Dynamic Analysis of Earth Dams*, UCB/EERC-81/14 — Eq. (3) complex modulus; §2.3 Tajirian 2DOF interpolation. https://nehrpsearch.nist.gov/static/files/NSF/PB82139940.pdf
* **[S8]** Schnabel, P.B., Lysmer, J., Seed, H.B. (1972). *SHAKE*, EERC 72-12, including Udaka & Lysmer revision of Sept. 1973. https://www.resolutionmineeis.us/sites/default/files/references/schnabel-lysmer-seed-1972.pdf
* **[S9]** Idriss, I.M., Sun, J.I. (1992). *User's Manual for SHAKE91*, and SHAKE91 Fortran source (public domain, "SHAKE16" repackaging). https://github.com/ocrickard/SHAKE16 (files `B1.f` CXSOIL/AMP, `C1.f` STRT/SHAKIT, `SHAKE91 User Manual.pdf`).
* **[S10]** Ostadan, F., Deng, N., Roesset, J.M. (2004). "Estimating total system damping for soil-structure interaction systems," 3rd UJNR SSI Workshop. https://www.pwri.go.jp/eng/ujnr/tc/a/ssi_w3/Contributions/Ostadan.pdf
* **[S11]** García, J.A., Kosbab, B., Tran, H., Richter, T., Rangelow, P. (2019). "Soil-pile-structure interaction analysis: implementation of flexible alternative theoretical approach in SASSI framework," SMiRT-25. https://scsolutions.com/wp-content/uploads/SMiRT25_Soil-Pile-Structure_Interaction_Analysis_Implementation_of_Flexible_Alternative_Theoretical_Approach_in_SASSI_Framework.pdf
* **[S12]** Ghiocel, D.M. (2019). *ACS SASSI Modeling, Part 1*, NRC training, June 25–27, 2019 (slides 9, 10, 14–16, 22–24). https://www.ghiocel-tech.com/storage/document/eb/50/82/a5aaf4c45eccd4453c4758c1e3c44a2d.nrc-acs-sassi-modeling-part-1-june-25-27-2019-.pdf
* **[S13]** Mertz, G.E., Costantino, M.C., Houston, T.W., Maham, A.S. (2011). "SASSI Subtraction Method Effects at Various DOE Projects," DOE NPH Workshop. https://www.energy.gov/sites/prod/files/SASSI%20Subtraction%20Method%20Effects%20at%20Various%20DOE%20Projects_1.pdf
* **[S14]** US DOE (2011). OE-3 2011-02, "SASSI Software Problem." https://www.energy.gov/sites/prod/files/2014/06/f16/OE-3_2011-02.pdf
* **[S15]** Anderson, L., Ostadan, F. (2014). "Validation of the SASSI2010 Subtraction Method Using Full Scale Independent Verification," DOE NPH Workshop. https://www.energy.gov/sites/prod/files/2014/12/f19/Anderson-%20Validation%20of%20SASSI%202010%20Subtraction%20Method.pdf
* **[S16]** STP 3&4 (2011). "Evaluation of D.O.E. SASSI Subtraction Method Analysis Recommendations" (f_SL, f_EV criteria). https://www.nrc.gov/docs/ML1121/ML112140192.pdf
* **[S17]** BNL (2008). Report N6112-051208 Rev.1, "Consistent Site Response – SSI…" (outcrop definitions). https://www.nrc.gov/docs/ML0919/ML091980384.pdf
* **[S18]** Ghiocel Predictive Technologies (2017). *ACS SASSI Version 3 User Manual* (local copy `reference/acs-sassi.txt`; see specs 01–05d).
* **[S19]** MTR & Associates, *MTR/SASSI Program Modules*. http://mtrassoc.com/mtrsassi/modules/
* Cited through the sources above, not accessed directly:
  - Waas, G. (1972), PhD thesis, UC Berkeley, *Linear two-dimensional analysis of soil dynamics problems in semi-infinite layered media*.
  - Lysmer, J., Waas, G. (1972), "Shear waves in plane infinite structures," *J. Eng. Mech. Div. ASCE* 98(EM1).
  - Lysmer, J., Kuhlemeyer, R.L. (1969), "Finite dynamic model for infinite media," *J. Eng. Mech. Div. ASCE* 95(EM4).
  - Kausel, E. (1974), *Forced Vibrations of Circular Foundations on Layered Media*, MIT R74-11.
  - Kausel, E., Roesset, J.M., Waas, G. (1975), *J. Eng. Mech. Div. ASCE* 101(EM5).
  - Kausel, E., Roesset, J.M. (1981), "Stiffness matrices for layered soils," *BSSA* 71(6).
  - Chen, J.-C., Lysmer, J., Seed, H.B. (1981), UCB/EERC-81/03.
  - Tabatabaie-Raissi, M. (1982), PhD thesis, UC Berkeley, *The flexible volume method for dynamic SSI analysis*.
  - Tajirian, F. (1981), PhD thesis, UC Berkeley, *Impedance matrices and interpolation techniques for 3-D interaction analysis by the flexible volume method*.
  - Ostadan, F., Deng, N. (2011), *SASSI2010*.
  - Ghiocel, D.M. (2013, 2015), FFV papers in SMiRT-22 and SMiRT-23. The NCSU repository blocks automated access, so these were not read.
