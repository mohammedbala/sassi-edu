---
id: 06-seismic-input
title: Seismic input: spectrum-compatible motion and strain-compatible soil
part: Design applications
order: 6
minutes: 35
example: ex04_site_response
summary: Turn a design response spectrum into an SSI input: a spectrum-compatible record checked against SRP 3.7.1, a strain-compatible soil profile from an equivalent-linear site response, and soil cases for property variation.
objectives: [Generate a spectrum-compatible motion with EQUAKE and read its acceptance table, Run an equivalent-linear (SHAKE) site response with SOIL and read the strain-compatible profile, Hand the strain-compatible soil to SITE with SITEX and FILE88, Run lower- and upper-bound soil cases and envelope their spectra]
prerequisites: [02-free-field, 04-surface-ssi]
---
In a fixed-base design you apply the design response spectrum at the base of the model, either
directly (response-spectrum analysis) or through a spectrum-compatible record (time-history analysis
for ISRS). An SSI analysis in SASSI needs two more things before the structure is even built:

1. a **time history**: SASSI works in the frequency domain and convolves its transfer functions with
   a record, so the design spectrum must become an acceleration history whose 5 % spectrum matches it;
2. a **soil profile at the strain level of that earthquake**: soil softens and dissipates more
   energy as it strains, and the strain depends on the motion. The low-strain shear-wave velocities of
   the geotechnical report are not the properties the structure sees.

This lesson runs example 4, the free-field chain **EQUAKE → SOIL → SITE**, on a 72 ft soil column
(30 ft of sand, 42 ft of clay, on rock) with no structure. EQUAKE builds the record, SOIL runs the
SHAKE equivalent-linear analysis and SITE takes over the strain-compatible properties for the SSI
analysis. The last two steps run the column again with softer and stiffer soil, the soil cases that
design practice envelopes. The lesson runs in about 20 s (EQUAKE takes most of it).

## Set the stage: the site, the time grid and the SSI frequencies

The model directory `ex04` holds every file of this lesson. The soil column is entered with its
**low-strain** properties, as given by the geotechnical investigation.

```sassi
MDL,ex04,ex04
TIT,Ex04 - free-field chain EQUAKE - SOIL - SITE (sand and clay on rock)
* L,<nm>,<thick>,<weight>,<Vp>,<Vs>,<pdamp>,<sdamp>  (low-strain properties)
L,1,3.0,0.120,1300,650,0.01,0.01
L,2,3.5,0.120,3000,1000,0.01,0.01
L,3,3.0,0.140,8000,4000,0.01,0.01
* 10 sublayers of sand 3 ft thick, 12 of clay 3.5 ft thick
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,2,2,2,2,2,2,2,2,2,2
TOPL,2,2
* frequency numbers n, f = n df with df = 1/(8192 x 0.005 s) = 0.0244 Hz
FREQ,1,4,20,40,60,80,100,120,140,160,180
FREQ,1,200,240,280,320,360,400,440,480,520,560
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
HOUSE,32.2,0,0,2,0,0,0,0,0
```

### What this does
* `L` defines three layer properties: sand ($V_s = 650\,\text{ft/s}$), clay ($V_s = 1{,}000\,\text{ft/s}$),
  both $0.120\,\text{kcf}$, and rock ($V_s = 4{,}000\,\text{ft/s}$, $0.140\,\text{kcf}$), all with 1 %
  damping, the small-strain values. `TOPL` stacks 22 sublayers, 10 of 3 ft in the sand and 12 of
  3.5 ft in the clay; L 3 is the half-space (SITE `<hs>` = 3), at 72 ft depth.
* `SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1` owns the time grid shared by every module:
  `<delt>` = 0.005 s, `<nft>` = 8192 Fourier points, so $\Delta f = 0.0244\,\text{Hz}$ and the Fourier
  period is 41 s, twice the 20 s record. The control point is the top of TOPL layer 1 (`<cl>` = 1),
  direction $x'$ (`<cm>` = 0); `WAVE,2,...` is a vertically propagating SV wave.
* `HOUSE,32.2,...` sets the gravity, 32.2 ft/s², and with it the British units (ft, kip, s) of the
  model; EQUAKE and the free-field modules take it from here.
* `FREQ` lists the SSI frequencies as integer numbers $n$ ($f = n\,\Delta f$): $n = 4$ (0.1 Hz), then
  multiples of 20, so that they fall on the 0.488 Hz grid of the SOIL amplification output used in
  step 6. The last one, $n = 560$ (13.7 Hz), is the cut-off.

### Why it matters
The time step and NFFT you choose here fix the frequency grid of every result, the length of the
quiet zone after the record (here 21 s) and the highest frequency the record carries (Nyquist
100 Hz). The cut-off is not arbitrary: 3 ft sublayers pass $V_s/(5h)$, and in step 5 the sand at
21-30 ft depth softens to 240-260 ft/s, which passes only 16-17 Hz (650 ft/s before softening:
43 Hz). A frequency set chosen with the low-strain velocities would exceed what the softened mesh
can carry.

### Technical basis
```math
\begin{aligned}
\Delta f &= \frac{1}{\Delta t\cdot\mathrm{NFFT}} = \frac{1}{0.005 \times 8192} = 0.0244\,\text{Hz}\\
f_\text{pass} &= \frac{V_s}{5h}
\end{aligned}
```
The second line is the one-fifth-wavelength rule (G-05) for a sublayer of thickness $h$; the SITE
listing prints $f_\text{pass}$ for every sublayer. See
[User Guide §9.1](docs/user/USER_GUIDE.md#91-time-step-nfft-and-the-frequency-step) and
[§7.2](docs/user/USER_GUIDE.md#72-layering-rules); the free-field theory is
[Theory §6](docs/theory/THEORY_MANUAL.md#6-free-field-motion-and-control-point-normalisation-site-mode-2)
(lesson 02).

### In ANSYS terms
`SITE <delt>` and `<nft>` are the `DELTIM` and the length of the transient, but chosen for an FFT: the
response is periodic with period $\mathrm{NFFT}\cdot\Delta t$, so the record needs a quiet tail in
which the response dies out, like the free-vibration tail you add to a transient run.

```action
plot-layers
explain: SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
```

## Generate a spectrum-compatible motion (EQUAKE)

The target is the RG 1.60 horizontal design spectrum anchored to 0.30 g at 5 % damping, given as 27
frequency / spectral-acceleration pairs in `data/rg160h_030g.rsi`. EQUAKE generates a 20 s record
whose 5 % spectrum matches it.

```sassi
* EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>
EQUAKE,0,27,11975,0.05,20,0,1
EQTIT,RG 1.60 horizontal, 5 %, anchored to 0.30 g
RSIN,1,../data/rg160h_030g.rsi
RSOUT,1,rg160h_eq.rso
ACCOUT,1,rg160h_eq.acc
AOPT,1,0,0,0,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNEQUAKE
```

### What this does
* `EQUAKE,0,27,11975,0.05,20,0,1`: random phases (`<accopt>` 0, no seed record), 27 target points
  (must equal the number of rows of the `.rsi` file), random seed 11975, target damping 5 %, duration
  20 s, no correlation between components, one seed trial.
* `RSIN` names the target, `RSOUT` the spectrum of the generated record, `ACCOUT` the record itself
  in g (EQUAKE also writes `.vel`, `.dis`, `.psd` and `.fft` with the same base name; with the
  gravity in ft/s², velocities are in in/s and displacements in inches).
* `AOPT,1,0,...` enables only EQUAKE for CHECK and AFWRITE; `RUNEQUAKE` runs the module (about 15 s).

### Why it matters
The record is the seismic input of everything that follows: the site response, the ISRS and the
member forces. A record whose spectrum dips below the target in a frequency band under-predicts the
response of any mode in that band, which is why the acceptance criteria limit how far, and over how
many adjacent frequencies, the spectrum may fall below the target.

### Technical basis
EQUAKE starts from a stationary Gaussian process with random phases and SIMQKE amplitudes under an
envelope, matches it in the frequency domain (Levy-Wilkinson: Fourier amplitudes multiplied by
$\mathrm{SA}_\text{target}/\mathrm{SA}_\text{computed}$, phases kept), refines it in the time domain
with the tapered-cosine wavelets of Al Atik and Abrahamson (2010) at the response-peak times, and
removes drift (high-pass filter and a constrained displacement polynomial with zero final velocity and
displacement). Response spectra are Nigam-Jennings exact recurrences at 100 points per decade.
[Theory §21](docs/theory/THEORY_MANUAL.md#21-spectrum-compatible-motions-equake); verified by
[VP-36](docs/verification/VERIFICATION_MANUAL.md#vp-36) (RG 1.60, SRP 3.7.1 criteria) and VP-30
(response-spectrum closed forms).

### Try this
Change the seed (`EQUAKE,0,27,12345,0.05,20,0,1`) and run again: a different record with the same
spectrum. Design practice uses several records or checks that one record is not unusually favourable;
the next step shows the quantities to compare.

```action
plot-spectrum: data/rg160h_030g.rsi, ex04/rg160h_eq.rso | log
plot-history: ex04/rg160h_eq.acc
open-listing: EQUAKE
explain: EQUAKE,0,27,11975,0.05,20,0,1
```

## Check the record against SRP 3.7.1

The EQUAKE listing ends with an acceptance table: the criteria of the ACS SASSI manual and those of
SRP 3.7.1 Rev. 4 (Option 1, Approach 2) reported separately. This step draws the spectral band of the
criteria so that you can see the match, not only read it.

```sassi
* line 1 = spectrum of the generated record, line 2 = target
READSPEC,rg160h_eq.rso,1,1
READSPEC,../data/rg160h_030g.rsi,1,2
* lines 3 and 4: the -10 % and +30 % limits around the target
LINECOMBIN,3,2,0.9
LINECOMBIN,4,2,1.3
WRITESPEC,srp_minus10.rs,3
WRITESPEC,srp_plus30.rs,4
```

### What this does
`READSPEC,<file>,<numLines>,<line>` loads spectrum columns into numbered line objects;
`LINECOMBIN,3,2,0.9` makes line 3, $y_3 = 0.9\,y_2$ ($y_k$ is the ordinate of line $k$); `WRITESPEC`
writes a line to a file (in the model directory). Plot them with the action below. The band is drawn
through the 27 target points, while EQUAKE checks against the log-log interpolated target at 100
points per decade, so the drawing is only a picture of the check.

### Why it matters
These are the numbers you put in a design report to show that the input is acceptable. For this
record the listing reports:

| Criterion (as checked by EQUAKE) | Required | This record |
|---|---|---|
| total duration | $\ge 20\,\text{s}$ | 20 s |
| time step / Nyquist | $\Delta t \le 0.005\,\text{s}$; Nyquist $\ge 50\,\text{Hz}$ | 0.005 s; 100 Hz |
| points per decade | $\ge 100$ | 100 |
| spectrum / target, 0.1-50 Hz (SRP) | none below 0.90, none above 1.30 | 0.957 to 1.166 |
| adjacent points below the target | $\le 9$ | 7 |
| strong-motion duration (Arias 5-75 %) | $\ge 6\,\text{s}$ | 9.03 s |

$\mathrm{PGA} = 0.324\,\text{g}$, $\mathrm{PGV} = 22.5\,\text{in/s}$, $\mathrm{PGD} = 16.0\,\text{in}$;
$V/A = 0.180\,\text{s}$ and $AD/V^2 = 3.95$, which the
listing says to compare with the controlling events. No target PSD was given, so the PSD criterion was
not checked (give one with `TPSD` and the eighth EQUAKE argument).

### Technical basis
The ACS SASSI manual's EQUAKE section lists the SRP criteria the module applies: duration
$\ge 20\,\text{s}$, $\Delta t \le 0.005\,\text{s}$, at least 100 points per decade, no point more than
10 % below or 30 % above the target, at most 9 adjacent points below. It names the strong-motion
duration (Arias 5-75 %), $V/A$ and $AD/V^2$ as computed parameters and says that a target PSD is
checked when one is given. The numerical limits for those, strong-motion duration $\ge 6\,\text{s}$,
cross-correlation $\lvert\rho\rvert \le 0.16$ between components and $\mathrm{PSD} \ge 80\,\%$ of the
target over 0.3-24 Hz, are the SRP 3.7.1 Rev. 4 checks as SASSI-EDU states them in
[Theory §21](docs/theory/THEORY_MANUAL.md#21-spectrum-compatible-motions-equake) item 6
(decision D-EQK-04: manual and SRP reported separately). The final velocity and displacement are zero
to round-off (-3e-13 in/s, -1.5e-12 in): no baseline drift.

### Check yourself
The spectrum of the record exceeds the target by up to 16.6 % (at 0.19 Hz). Is that a problem for the
SSI analysis of a stiff nuclear building whose SSI frequencies are 3-10 Hz?

Answer: Not for the acceptance criteria (the limit is +30 %), and it is at a frequency far below the
structure's modes. The band that drives the ISRS is 3-10 Hz, and there the record is closest to the
target: the ratio runs from its overall minimum, 0.957 at 4.67 Hz (a narrow dip), to about 1.12, with
a mean of 1.05 (from `rg160h_eq.rso`), inside the -10 % allowance. The dip, not the 0.19 Hz excess,
is what to look at for this building.

```action
plot-spectrum: ex04/rg160h_eq.rso, ex04/srp_minus10.rs, ex04/srp_plus30.rs | log
plot-history: ex04/rg160h_eq.vel
plot-history: ex04/rg160h_eq.dis
open-listing: EQUAKE
```

## Strain-dependent soil curves

SOIL needs, for each soil, how the shear modulus and the damping change with shear strain. Example 4
reads the SHAKE91 curves (Seed and Idriss 1970 for sand, Seed and Sun 1989 for clay, Idriss 1990
damping) from the library `sassi/data/dynp_library.pre`, copied into this lesson's workspace. `SPRO`
then assigns an L property and a curve to every sublayer of the SOIL column.

```sassi
* the DYNP curves Clay, Sand and Rock of SHAKE91
INP,../sassi/data/dynp_library.pre
* SPRO,<sublayer>,<L property>,<curve>: 1-10 sand, 11-22 clay, the last entry is the half-space
VAR,SANDL,1,2,3,4,5,6,7,8,9,10
VAR,CLAYL,11,12,13,14,15,16,17,18,19,20,21,22
FOREACH,SANDL,SPRO,@SANDL[#],1,Sand
FOREACH,CLAYL,SPRO,@CLAYL[#],2,Clay
SPRO,23,3
```

The library is a file of `DYNP` commands; the sand curve, for example, is

```sassi-show
* DYNP,<no>,<strain %>,<G/Gmax>,<strain %>,<damping %>,<label>
DYNP,1,0.0001,1.000,0.0001,0.24,Sand
DYNP,5,0.01,0.850,0.01,2.80,Sand
DYNP,7,0.1,0.370,0.1,9.80,Sand
DYNP,9,1.0,0.080,1.0,21.00,Sand
DYNP,11,10.0,0.035,10.0,28.00,Sand
```

### What this does
* `INP` executes the commands of a `.pre` file. Each `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>` of the
  library sets point `<no>` (of 11) of two curves of the property `<label>`: $G/G_{\max}$ at shear
  strain `<sg>` and damping (in percent) at shear strain `<sd>`, strains in percent.
* `FOREACH,SANDL,SPRO,@SANDL[#],1,Sand` runs `SPRO,k,1,Sand` for every sublayer $k$ of the list
  SANDL. `SPRO,23,3` without a label is the half-space: it stays linear.

### Why it matters
The curves carry most of the uncertainty of a site response: with the same record and the same
low-strain profile, a different curve gives a different strain-compatible profile. They come from
the geotechnical report, like the low-strain velocities, and design practice varies them together with
the profile.

### Technical basis
SOIL interpolates both curves linearly in $\log_{10}(\text{strain})$ and holds them constant outside
the tabulated range, as SHAKE91 does (D-SOL-03,
[Theory §12.2](docs/theory/THEORY_MANUAL.md#122-strains-and-the-equivalent-linear-iteration)).

### Check yourself
At 0.1 % shear strain, by how much has the sand lost stiffness, and what damping does it have?

Answer: $G/G_{\max} = 0.37$ (63 % of the low-strain shear modulus is lost, $V_s$ falls to
$\sqrt{0.37} = 61\,\%$ of its low-strain value) and the damping is 9.8 %, ten times the 1 % of the
low-strain profile.

```action
plot-soilprops: Sand
plot-soilprops: Clay
open-file: sassi/data/dynp_library.pre
explain: DYNP,7,0.1,0.370,0.1,9.80,Sand
```

## Equivalent-linear site response (SOIL)

SOIL applies the EQUAKE record as a **rock-outcrop** motion at the top of the half-space, computes the
strains in every sublayer, reads new properties from the curves and repeats.

```sassi
* SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>
SOIL,4000,32.2,1,1,1,8,0.65,1,0
* SOILX,<indir>,<mult>,<max>,<cl>: horizontal, factor 1, input at the top of sublayer 23
SOILX,0,1,0,23
THFILE,rg160h_eq.acc
THTIT,RG 1.60 spectrum-compatible rock outcrop motion (EQUAKE)
* outputs: surface history, outcrop maximum, spectra, strains, amplification
SACC,1,2,0
SACC,23,1,1
SRS,1,1,0
SRS,23,1,1
SSTR,5,1,0,1,1
SSTR,15,1,0,1,1
SSAF,1,1,0,0,23,0.48828125,surface / base (within)
DAMP,0.05
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSOIL
```

### What this does
* `SOIL,4000,32.2,1,1,1,8,0.65,1,0`: read 4000 values after 1 header line (the dt line of the EQUAKE
  file), gravity 32.2 ft/s² (SOIL takes its units from it: ft, kcf, ksf, ft/s), the input is an
  **outcrop** motion (`<outcrop>` 1), save the strain-compatible properties to FILE88
  (`<save>` 1), 8 iterations, effective strain
  $\gamma_\text{eff} = 0.65\,\gamma_{\max}$ ($\gamma_{\max}$ the maximum strain).
* `SOILX,0,1,0,23`: horizontal input ($V_s$ and shear damping), factor 1, applied at the top of
  sublayer 23, the rock.
* `SACC,1,2,0` saves the surface (within) acceleration history `ACC001.TH`; `SACC,23,1,1` reports
  the rock-outcrop maximum; `SRS` the 5 % spectra (`DAMP,0.05`) at the surface and of the outcrop;
  `SSTR` the strain histories at mid-height of sublayer 5 (sand) and 15 (clay); `SSAF` the surface /
  base Fourier amplification at 0.488 Hz steps.

### Why it matters
The listing `ex04_SOIL.out` gives the profile the structure will actually sit on. For this record:

* the iterations settle: the largest change of the shear modulus falls from 178 % (iteration 1) to
  2.0 % (iteration 8), of the damping from 90 % to 0.7 %;
* the sand softens with depth, from $V_s = 620\,\text{ft/s}$ at the surface ($G/G_{\max} = 0.91$,
  2.0 % damping) to **242 ft/s at 28.5 ft** ($G/G_{\max} = 0.14$, 17.8 % damping, effective strain
  0.49 %); the clay keeps $V_s = 843\text{–}903\,\text{ft/s}$ with 5.9-8.4 % damping;
* the column period lengthens from 0.34 s (low strain) to 0.43 s; the surface / base (within)
  amplification peaks at **4.96 at 1.97 Hz**, against 69 at 3.22 Hz with the low-strain properties
  and their 1 % damping;
* the surface motion peaks at **0.627 g** for a rock-outcrop input of 0.324 g; its 5 % spectrum peaks
  at 2.64 g at 2.04 Hz.

Applying the rock-outcrop design motion directly at the ground surface would miss this: the soil
moves the spectral peak from 2.5 Hz (target) to 2.04 Hz and nearly doubles the peak acceleration at
the ground surface ($0.627/0.324 = 1.94$), where a surface foundation would sit.

### Technical basis
In each sublayer a vertically propagating shear wave is the sum of an up-going and a down-going wave;
continuity of displacement and shear stress gives the recursion (12.2) from the free surface down. The
**within** motion at the top of sublayer $m$ is $E_m + F_m$, the **outcrop** motion $2E_m$, with $E_m$
and $F_m$ the amplitudes of the up-going and the down-going wave. Each iteration:

```math
\begin{aligned}
\gamma_\text{eff} &= R_\gamma\,\max_t \lvert\gamma(t)\rvert\\
G &= G_{\max}\,(G/G_{\max})(\gamma_\text{eff}), \qquad \beta = D(\gamma_\text{eff})
\end{aligned}
```

with $R_\gamma = 0.65$ here (SOIL `<ratio>`; SHAKE91 suggests $(M - 1)/10$ for magnitude $M$),
$(G/G_{\max})(\gamma)$ and $D(\gamma)$ the modulus-reduction and damping curves. SOIL runs exactly
`<iter>` iterations (SHAKE91 has no convergence test) and lists the changes per iteration.
[Theory §12](docs/theory/THEORY_MANUAL.md#12-the-shake-equivalent-linear-method-soil); verified by
[VP-04](docs/verification/VERIFICATION_MANUAL.md#vp-04) (the published SHAKE91 sample problem) and
VP-02a (closed-form layer transfer functions).

### In ANSYS terms
SOIL is a one-dimensional shear column whose secant modulus and damping are updated from the strain of
the previous pass, like a nonlinear material handled by repeated linear (harmonic) solutions with
secant properties instead of a Newton-Raphson time integration.

### Try this
Change the effective strain ratio to 0.5 (`SOIL,4000,32.2,1,1,1,8,0.5,1,0`), run AFWRITE and RUNSOIL
again, and compare the strain-compatible $V_s$ at 28.5 ft and the surface PGA in the listing. Then restore
0.65 and run AFWRITE and RUNSOIL once more before the next step: FILE88 and `ACC001.TH` always hold
the last SOIL run, and the next steps use them.

```action
open-listing: SOIL
plot-spectrum: ex04/RS023_01.RS, ex04/RS001_01.RS | log
plot-history: ex04/ACC001.TH
plot-history: ex04/SN005.TH, ex04/SN015.TH
explain: SOIL,4000,32.2,1,1,1,8,0.65,1,0
```

## Hand the strain-compatible soil to SITE (SITEX,1)

SOIL wrote the converged properties of the 22 sublayers to FILE88. `SITEX,1` (Non-Linear Soil) makes
SITE read them instead of the low-strain L table, so the free field of the SSI analysis is the one of
the earthquake.

```sassi
* SITEX,<soilmode>: 1 = SITE reads FILE88 written by SOIL
SITEX,1
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
```

### What this does
`SITEX,1` switches SITE to FILE88, position by position (the number of sublayers must match, EDU-07).
`AOPT` now enables SITE only. The SITE listing opens with "FILE88 strain-compatible properties read
for 22 layers" and prints the layer table with the iterated $V_s$, damping and the passing frequency
$f_\text{pass}$ of each sublayer.

### Why it matters
* **Consistency.** The SSI free field (FILE1) must reproduce the site response that defined the
  motion. Here it does: at 1.95 Hz the SITE surface / base amplification is $1/0.2021 = 4.949$, the
  SOIL value (`SAF001W_023W.TFU`) 4.954, 0.1 % apart; examples/README.md reports below 1 % up to
  7.8 Hz and 3.4 % at 13.7 Hz, the thin-layer discretisation error growing as $(kh)^2$.
* **Passing frequency.** The SITE listing shows $f_\text{pass} = 16.1\,\text{Hz}$ for the softened sand
  at 27-30 ft:
  the 13.7 Hz cut-off chosen in step 1 is inside it, 15 % below.
* **The control motion of the SSI run.** The SSI control point here is the top of layer 1, the free
  surface (at the free surface the within and the outcrop motions coincide). FILE1 is normalised to a
  unit motion there, so the SSI analysis must be driven by the SOIL surface motion `ACC001.TH`
  (`THFILE,ACC001.TH` in MOTION), not by the rock-outcrop record. Driving it with the rock record is
  a mistake, and here an unconservative one: 0.324 g instead of 0.627 g at the surface.

### Technical basis
FILE88 holds per sublayer the thickness, effective strain, $G$, $V_s$, $\beta_s$, $V_p$ and $\beta_p$.
For $V_p$ and $\beta_p$ the default policy keeps Poisson's ratio and sets $\beta_p = \beta_s$
(`EDUOPT,VPPOLICY`; use `VP` to keep $V_p$ for saturated soils). SOIL and SITE must use the same
complex-modulus form, otherwise the SSI free field would not reproduce the SOIL motion at the control
point. [Theory §12.3](docs/theory/THEORY_MANUAL.md#123-hand-off-to-ssi),
[User Guide §7.5](docs/user/USER_GUIDE.md#75-strain-compatible-properties-soil-and-sitex). For an
embedded structure the excavated soil must also get these properties (lesson 09 shows how).

### Check yourself
Your SSI model has its control point at the free surface (`SITE <cl>` = 1). Which record do you
give to MOTION: `rg160h_eq.acc` or `ACC001.TH`? And what if the design spectrum were defined at the
foundation level of a 30 ft deep basement?

Answer: `ACC001.TH`, the SOIL surface (within) motion, because FILE1 is normalised to a unit motion
at the control point. For a spectrum defined at foundation level, put the control point at the TOPL
interface at that depth. SITE's control motion is a *within* motion (lesson 02), so if the design
motion is defined there as an outcrop, run SOIL with it as an outcrop at that depth and drive the SSI
with the within motion SOIL computes there (`SACC` of that sublayer, within), not with the outcrop
record itself.

```action
open-file: ex04/FILE88
open-listing: SITE
plot-spectrum: ex04/SAF001W_023W.TFU
plot-history: ex04/ACC001.TH
```

## Soil property variation: lower- and upper-bound cases

The low-strain profile and the curves are uncertain. Design practice therefore runs the analysis for
several soil cases (best estimate, lower bound, upper bound) and envelopes the results. Each case is a
complete chain (SOIL → SITE → SSI); here we run SOIL for two more cases, in their own model
directories, with the same record.

```sassi
* copy the best-estimate model to models 1 and 2
CPMODEL,1
CPMODEL,2
* lower bound: G_max / 1.5, i.e. Vs and Vp / 1.225
ACTM,1
MDL,ex04lb,../ex04lb
L,1,3.0,0.120,1061,531,0.01,0.01
L,2,3.5,0.120,2449,816,0.01,0.01
THFILE,../ex04/rg160h_eq.acc
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
AFWRITE
RUNSOIL
* upper bound: G_max x 1.5, i.e. Vs and Vp x 1.225
ACTM,2
MDL,ex04ub,../ex04ub
L,1,3.0,0.120,1593,796,0.01,0.01
L,2,3.5,0.120,3675,1225,0.01,0.01
THFILE,../ex04/rg160h_eq.acc
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
AFWRITE
RUNSOIL
ACTM,0
```

### What this does
* `CPMODEL,1` copies the active model (with its name and path, hence the warning) to model 1;
  `ACTM,1` activates it and `MDL,ex04lb,../ex04lb` gives it its own name and directory, so its files
  do not overwrite the best-estimate ones.
* The two `L` commands scale the low-strain $V_s$ and $V_p$ of sand and clay (Poisson's ratio
  unchanged); the curves, the rock and the record are the same. The record lives in `ex04`, hence
  `THFILE,../ex04/rg160h_eq.acc`.
* The factor 1.5 on $G_{\max}$ is **illustrative only**: neither the ACS SASSI manual nor the SASSI-EDU
  documents prescribe a value. The bounds of a project come from its geotechnical basis and the
  regulatory guidance that applies to it.

### Why it matters
Compare the three listings (`ex04_SOIL.out`, `ex04lb/ex04lb_SOIL.out`, `ex04ub/ex04ub_SOIL.out`):

| Soil case | Surface PGA | Peak of the 5 % surface spectrum | Surface / base amplification peak |
|---|---|---|---|
| lower bound | 0.557 g | 2.02 g at 1.51 Hz | 4.28 at 1.46 Hz |
| best estimate | 0.627 g | 2.64 g at 2.04 Hz | 4.96 at 1.97 Hz |
| upper bound | 0.678 g | 3.31 g at 2.88 Hz | 6.96 at 2.83 Hz |

The site frequency moves with the stiffness, and the response does not scale with it: the upper
bound moves the site frequency (2.83 Hz) onto the high plateau of the RG 1.60 input spectrum, just
above its 2.5 Hz peak, and its surface spectrum peaks 25 % above the best estimate. Every structure
frequency near the site frequency sees a different amplification in each case. That is why the ISRS of the soil cases are
enveloped and broadened (lesson 07), and why the governing soil case is often a different one for
different floors and frequency bands.

Read the last line of each iteration table too: the upper-bound case still changes its moduli by
4.0 % at iteration 8 (best estimate 2.0 %, lower bound 0.6 %). Iterations are cheap; when the last
change is not small, increase `<iter>`.

### Technical basis
Each case is an independent equivalent-linear analysis, so its strain-compatible profile is not the
best-estimate profile multiplied by the factor. In the upper bound the sand at 28.5 ft strains 0.25 %
instead of 0.49 % and keeps $V_s = 363\,\text{ft/s}$, not $1.225 \times 242 = 296\,\text{ft/s}$. Bounds
are therefore applied to the **low-strain** properties (or to the curves) and the site response is
repeated. [User Guide §4.2](docs/user/USER_GUIDE.md#42-models-in-memory) explains the models in memory
(`CPMODEL`, `ACTM`, `MDL`).

```action
open-file: ex04lb/ex04lb_SOIL.out
open-file: ex04ub/ex04ub_SOIL.out
plot-spectrum: ex04lb/RS001_01.RS, ex04/RS001_01.RS, ex04ub/RS001_01.RS | log
```

## Envelope the soil cases

The line operations of the Plot menu combine spectra. `BROADEN` forms the envelope of several lines;
with a peak broadening of 0 % it is the plain envelope.

```sassi
READSPEC,../ex04lb/RS001_01.RS,1,21
READSPEC,RS001_01.RS,1,22
READSPEC,../ex04ub/RS001_01.RS,1,23
LINENAME,21,lower bound
LINENAME,22,best estimate
LINENAME,23,upper bound
* BROADEN,<dest>,<bridging %>,<broadening %>,<sources>: envelope only
BROADEN,24,0,0,21,22,23
WRITESPEC,surface_envelope.rs,24
```

### What this does
`READSPEC` loads the surface spectrum of each case into lines 21-23 (paths are relative to the
active model's directory, `ex04`), `LINENAME` labels them, `BROADEN,24,0,0,21,22,23` stores their
envelope in line 24 and `WRITESPEC` writes it to `ex04/surface_envelope.rs`.

### Why it matters
The envelope of the soil cases is the design spectrum at the ground surface for this column; in an SSI
analysis the same operation is applied to the ISRS of each soil case at each floor, followed by peak
broadening. Here the lower bound governs from 0.25 to 1.7 Hz, the best estimate from 1.7 to 2.3 Hz and
the upper bound above 2.3 Hz; the envelope peaks at 3.31 g at 2.88 Hz, from the upper bound.

### Technical basis
`BROADEN` computes, on the union of the frequencies of its sources, the envelope
$E(f) = \max_i y_i(f)$ ($y_i$ the ordinates of the source lines), then broadens its peaks by
$\pm\text{Smooth2}\,\%$ and bridges valleys between peaks within Smooth1 %
([User Guide §11.5](docs/user/USER_GUIDE.md#115-line-mathematics-and-spectrum-broadening)); verified
by [VP-48](docs/verification/VERIFICATION_MANUAL.md#vp-48). The ACS SASSI manual recommends, per SRP
3.7.1, spectra from 0.1 to 100 Hz with at least 301 frequencies and warns that fewer reduce the
accuracy of its broadening algorithm (the SOIL spectra here have 301).

### Check yourself
Why not simply run the best-estimate case and broaden its spectrum by ±15 % instead of running the
bounds?

Answer: Broadening covers the frequency shift of a peak, but not the change of its amplitude: the
cases differ in damping (a softer soil strains more and damps more) and in amplification, and the
soft case can also change which soil layer controls the response. Broadening is applied in addition to
the soil cases, to the enveloped ISRS, not instead of them.

```action
plot-spectrum: ex04lb/RS001_01.RS, ex04/RS001_01.RS, ex04ub/RS001_01.RS, ex04/surface_envelope.rs | log
open-doc: docs/user/USER_GUIDE.md#75-strain-compatible-properties-soil-and-sitex
```
