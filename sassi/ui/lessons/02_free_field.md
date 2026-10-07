---
id: 02-free-field
title: The site: layers, half-space and free-field motion
part: Fundamentals
order: 2
minutes: 25
example: ex01_surface_stick
summary: Build the layered site of example 1, run SITE, and check its free-field transfer functions against a one-dimensional SHAKE solution.
objectives: [Define soil layers and the layered column, Apply the one-fifth-wavelength rule to sublayer thickness, Set the wave field and the control point, Read the SITE listing (modes and half-space and free field), Tell within from outcrop motion and see why the half-space matters]
prerequisites: [01-why-ssi]
---
In a fixed-base analysis the design motion is applied directly at the base of the model. In SASSI
the earthquake reaches the structure through the soil: SITE computes the **free-field motion**,
the motion of the horizontally layered site *without* the structure, at every layer interface and
every analysis frequency, normalised so that the **control point** moves with unit amplitude. The
structure later receives it through the interaction nodes.

This is the same one-dimensional wave propagation you know from site-response analyses (SHAKE,
STRATA), with two additions SASSI needs for SSI: the layered soil is discretised into thin
sublayers so that the same model also gives the soil's response to point loads (the soil springs
of lesson 3), and the half-space below the layers is simulated with extra sublayers and dashpots.

In this lesson you build the site of example 1 (16 ft of sand and 39 ft of gravel on weathered
rock), run SITE, read its listing, and then check the free field against an exact one-dimensional
(SHAKE-type) solution of the same column computed by the SOIL module.

## Units and soil layer types

```sassi
MDL,ex01,ex01
TIT,Ex01 site - 16 ft of sand and 39 ft of gravel on weathered rock
* US customary units: ft, kip, s; weights are unit weights in kcf (kip/ft3), divided by gravity
GRAVITY,32.2
* L,<nm>,<thick>,<weight>,<Vp>,<Vs>,<pdamp>,<sdamp>
L,1,1.6,0.120,2000,1000,0.05,0.05
L,2,3.25,0.125,3300,1650,0.04,0.04
L,3,3.25,0.130,6600,3300,0.02,0.02
LLIST,1,3
```

```action
explain: L,1,1.6,0.120,2000,1000,0.05,0.05
```

### What this does
`MDL,ex01,ex01` names the model and creates its directory `ex01/`, which becomes the working
directory. `GRAVITY,32.2` sets the acceleration of gravity and with it the unit system: SASSI works
in any consistent unit set, takes **unit weights** and divides them by gravity to get densities (a
new model already starts with 32.2 ft/s², the ft-kip-s value; setting it explicitly documents the
units of the model).

`L` defines a soil layer **type**: `<nm>` its number, `<thick>` the thickness of each sublayer that
uses it, `<weight>` the unit weight, `<Vp>` and `<Vs>` the compression and shear-wave velocities,
`<pdamp>` and `<sdamp>` the P- and S-wave damping ratios (thickness in ft, unit weight in kcf =
kip/ft³, velocities in ft/s). Layer 1 is dense sand ($V_s = 1{,}000\,\text{ft/s}$, 5 %), layer 2
gravel ($V_s = 1{,}650\,\text{ft/s}$, 4 %), layer 3 weathered rock ($V_s = 3{,}300\,\text{ft/s}$, 2 %)
that will be the half-space; its thickness is not used. `LLIST` prints the derived properties:
$G = 3{,}727\,\text{ksf}$, 10,569 ksf and 43,966 ksf, mass densities 0.00373, 0.00388 and
0.00404 kip·s²/ft⁴, $\nu = 1/3$ for all three.

### Why it matters
$V_s$ and damping are the soil properties that control SSI. In a design analysis they are the
strain-compatible values from the site-response analysis (SOIL, lesson 6) for the lower-bound,
best-estimate and upper-bound soil cases, not the low-strain values from the geophysical survey.
Poisson's ratio follows from $V_p/V_s$; saturated soils with $\nu$ above 0.47 need care (CHECK
warns).

### Technical basis
$G = \rho V_s^2$, the constrained modulus $M = \lambda + 2G = \rho V_p^2$. Damping is hysteretic
(independent of frequency) and enters as a complex modulus on $G$ (S-wave damping) and on $M$
(P-wave damping):

```math
G^* = G\left(1 - 2\beta^2 + 2i\beta\sqrt{1-\beta^2}\right), \qquad \lvert G^*\rvert = G
```

with $\beta$ the damping ratio
([Theory §2](docs/theory/THEORY_MANUAL.md#2-complex-modulus-damping); VP-01, VP-02a).

### In ANSYS terms
A layer type is a material (`MP,EX`, `MP,PRXY`, `MP,DENS`) with a material damping ratio
(`MP,DMPR`). If you reproduce the SASSI complex modulus with ANSYS structural damping (stiffness
multiplied by $1 + ig$), use $E_{\text{ANSYS}} = E(1-2\beta^2)$ and
$g = 2\beta\sqrt{1-\beta^2}/(1-2\beta^2)$
([Theory §2](docs/theory/THEORY_MANUAL.md#2-complex-modulus-damping), section 2.3).

## The layered column and the one-fifth-wavelength rule

```sassi
* TOPL lists the layers from the ground surface down (it appends to the list)
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,2,2,2,2,2,2,2,2,2,2
TOPL,2,2
```

```action
plot-layers
explain: TOPL,1,1,1,1,1,1,1,1,1,1
```

### What this does
`TOPL` builds the column from the ground surface down: ten 1.6 ft sublayers of type 1 (the 16 ft
of sand) and twelve 3.25 ft sublayers of type 2 (the 39 ft of gravel), 22 sublayers in all. The
**interfaces** are the tops of the sublayers: interface 1 is the ground surface, interface 23 the top
of the half-space at 55 ft depth. The half-space is not part of `TOPL`; SITE adds it below
(next step). The layer plot shows the column and the property table.

### Why it matters
The sublayers are not the geological layers: they are the vertical finite-element mesh of the
soil. Their thickness limits the highest frequency the model transmits, exactly as the element
size does in an ANSYS soil model. Interaction nodes of an embedded structure must lie on
interfaces, so the sublayering also has to match the basement levels (lesson 5).

### Technical basis
Inside a sublayer the displacement varies linearly with depth (thin-layer method, a finite-element
discretisation in depth only). With the mixed mass matrix
($\tfrac12\,\text{lumped} + \tfrac12\,\text{consistent}$) used by SITE, the rule of the manual is
that a sublayer be no thicker than a fifth of the shortest wavelength:

```math
h \le \frac{V_s}{5\,f_{\text{cut}}}, \qquad f_{\text{pass}} = \frac{V_s}{5\,h}
```

with $h$ the sublayer thickness, $f_{\text{cut}}$ the cut-off frequency and $f_{\text{pass}}$ the
passing frequency of the sublayer. Here:

* sand: $f_{\text{pass}} = 1000/(5 \times 1.6) = 125\,\text{Hz}$;
* gravel: $f_{\text{pass}} = 1650/(5 \times 3.25) = 102\,\text{Hz}$.

Treat $f_{\text{pass}}$ as a limit, not a target: exactly at $h = \lambda/5$ the one-dimensional
amplification can still be off by up to about 9 % with the mixed mass (17 % with a consistent mass);
the error falls quickly with thinner sublayers, as the fourth power of $h/\lambda$ for the mixed
mass ([Theory §4.1](docs/theory/THEORY_MANUAL.md#4-the-thin-layer-method-site-mode-1)). Both layers
here pass the 20 Hz cut-off of this analysis by a factor of five or more. The manual also asks for
**more than 20 sublayers**, because too few degrade the Rayleigh and Love modes on which the
point-load solutions are built (22 here). CHECK warns when a layer does not pass the cut-off
(EDU-10) ([User Guide §7.2](docs/user/USER_GUIDE.md#72-layering-rules); VP-03 checks the dispersion
of linear elements with lumped, consistent and mixed mass, VP-05 the Rayleigh-wave velocity of a
deep stratum).

### Check yourself
A soft clay with a strain-compatible $V_s$ of 500 ft/s must pass a 30 Hz cut-off. What is the
largest sublayer thickness?

Answer: $h = 500/(5 \times 30) = 3.3\,\text{ft}$. Use the strain-compatible $V_s$, not the
low-strain value: if SOIL softens the clay to 330 ft/s, the limit drops to 2.2 ft.

## Frequencies, the wave field and the control point

```sassi
* frequency numbers n (f = n x df); 22 frequencies from 0.1 to 20 Hz
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,287,328,369,410,492,573,655
FREQ,1,737,819
LFREQ
* SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
* WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>: a vertically incident SV wave
WAVE,2,1,1,1,0
```

```action
explain: SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
open-dialog: ANALYSIS/SITE
```

### What this does
SASSI's frequencies are integer **frequency numbers** $n$ times the frequency step
$\Delta f = 1/(\Delta t \cdot \text{NFFT}) = 1/(0.005 \times 8192) = 0.0244\,\text{Hz}$; `LFREQ`
lists them in Hz (lesson 4 explains how to choose them). The `SITE` arguments:

| Argument | Value | Meaning |
|---|---|---|
| `<opmode>` | 0 | solution (1 = data check only) |
| `<mode1>`, `<mode2>` | 1, 1 | run Mode 1 (layer eigen-solutions, FILE2) and Mode 2 (free field, FILE1) |
| `<fstep>` | 0 | $\Delta f$ from $\Delta t$ and NFFT |
| `<nl>`, `<hs>` | 20, 3 | 20 generated half-space sublayers with the properties of layer type 3 (0 = rigid base) |
| `<wopt>` | 0 | in-plane waves (P, SV, Rayleigh); 1 = anti-plane (SH, Love) |
| `<freq1>`, `<freq2>` | 1, 4096 | frequency numbers at which the WAVE ratios are given |
| `<cl>`, `<cm>` | 1, 0 | control point at the top of TOPL layer 1 (the ground surface), direction $x'$ |
| `<delt>`, `<nft>`, `<freq>` | 0.005, 8192, 1 | time step, NFFT, frequency set 1 |

`WAVE,2,1,1,1,0` defines the wave field: type 2 = SV, `<opt>` 1 switches the wave on, ratio 1 at
both frequencies (the only wave), incidence angle 0 (vertical).

### Why it matters
Vertically propagating shear waves for the horizontal input and compression waves for the vertical
input are the standard assumption of design SSI analyses; the ACS SASSI manual recommends SV for X,
SH for Y and P for Z (lesson 7 runs all three). The control point is where your design motion is
defined. Here it is the ground surface, the usual location of a design (surface) response
spectrum.

### Technical basis
SITE normalises every wave field so that the control point moves with unit amplitude in the
control direction:

```math
\hat{U}(z) = \frac{u(z)}{u_{\text{cm}}(z_{\text{cp}})}
```

for one wave, with $u(z)$ the wave field at depth $z$ and $u_{\text{cm}}(z_{\text{cp}})$ its
component in the control direction at the control point (with several waves the ratios sum to 1), so
FILE1 holds transfer functions "free-field motion per unit control motion". The control motion is
the **within** motion at the top of layer `<cl>`: the actual motion of that point inside the column,
incident plus reflected waves
([Theory §6.2](docs/theory/THEORY_MANUAL.md#62-normalisation-to-the-control-point)).

### Check yourself
Why is the ground surface the most convenient control point?

Answer: at the free surface the up-going and down-going waves have equal amplitude, so the
within motion and the outcrop motion are the same. At depth they differ, and the control motion
must then be the within motion of that depth (last step of this lesson).

## Run SITE: modes, half-space and free field

```sassi
* only SITE: AOPT,<EQUAKE>,<SOIL>,-,<SITE>,<POINT>,<HOUSE>,-,<FORCE>,<ANALYS>,...
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
```

```action
open-listing: SITE
open-file: ex01/ex01.sit
```

### What this does
`AOPT` selects the modules that `CHECK` checks and `AFWRITE` writes decks for; here only SITE.
`AFWRITE` writes `ex01.sit`, and `RUNSITE` runs the module in the model directory (about a
second). It writes FILE2 (Mode 1) and FILE1 (Mode 2) and the listing `ex01_SITE.out`, which has
four parts:

1. **Soil layers**: the 22 sublayers and the half-space with $V_s$, damping, $\nu$ and the passing
   frequency $f_{\text{pass}}$ (125 and 102 Hz).
2. **Generated half-space sublayers**: below 55 ft SITE adds 20 sublayers of rock whose total
   depth is 1.5 shear wavelengths of the rock: 50,688 ft at 0.098 Hz (2,534 ft each), 4,945 ft at
   1.0 Hz, 248 ft at 20 Hz (12.4 ft each).
3. **Mode 1**: the propagating Rayleigh and Love modes of the column. At 20 Hz there are 7
   propagating Rayleigh modes (of 86) and 5 Love modes (of 43); the slowest propagating Rayleigh
   mode has a phase velocity of 1,352 ft/s and a wavelength of 67.6 ft.
4. **Mode 2**: the free-field amplitude at every interface for unit motion at the surface. At
   0.098 Hz it is 1.000 down to 55 ft. At 20 Hz it is 0.98 at 1.6 ft, falls to 0.09 at 12.8 ft,
   rises again to 0.68 at 29 ft and is 0.37 at 55 ft. The last line of each frequency is the ratio
   of the surface motion to the **outcrop** motion at the top of the half-space: 1.00 at 0.1 Hz,
   **2.13 at 7.0 Hz** and 2.11 at 8.0 Hz (the column frequency lies between the two), 1.48 at
   12 Hz, 2.18 at 18 Hz.

### Why it matters
FILE2 (the modes) feeds POINT, which computes the soil springs and dashpots (lesson 3). FILE1 is
the seismic input of ANALYS. The amplitudes with depth matter for embedded structures: at 20 Hz
the free field at about 13 ft depth is almost at rest (0.09) while the surface moves fully, so a
basement that spans these depths receives a strongly varying motion (kinematic interaction,
lesson 5).

### Technical basis
**Mode 1** solves the thin-layer eigenproblem $(A k^2 + B k + G - \omega^2 M)\,\phi = 0$ for the
wavenumbers $k$ of the free waves of the discrete column
([Theory §4](docs/theory/THEORY_MANUAL.md#4-the-thin-layer-method-site-mode-1)). **The half-space**
is simulated by two techniques together, as in the ACS SASSI manual: the variable-depth method
(sublayers with a total depth of $1.5\,V_{s,\text{hs}}/f$, where the fundamental Rayleigh mode has
died out) and viscous dashpots $\rho V_s$, $\rho V_p$ (with complex velocities) at their base
([Theory §5](docs/theory/THEORY_MANUAL.md#5-half-space-variable-depth-and-viscous-boundary)). The
manual recommends 10 to 20 generated sublayers. VP-02b compares SITE with the exact solution of a
layer on a half-space: 0.20 % with 20 uniform sublayers, 0.95 % with 10
([VP-02b](docs/verification/VERIFICATION_MANUAL.md#vp-02b)). **Mode 2** solves the column for
vertically incident waves ($k = 0$) and normalises to the control point
([Theory §6](docs/theory/THEORY_MANUAL.md#6-free-field-motion-and-control-point-normalisation-site-mode-2)).

### In ANSYS terms
Mode 1 is a modal analysis of the soil column, but in the wavenumber domain: for each frequency it
finds the wavelengths of the waves the layered soil can carry. There is no ANSYS equivalent of the
half-space buffer; in a direct (box-of-soil) ANSYS model you would need absorbing boundaries at the
base and sides.

## Free-field transfer functions: a one-dimensional SHAKE check

The SITE listing gives the free field only at the 22 analysis frequencies. The SOIL module (the
SHAKE methodology, lesson 6) gives the same column as a continuous transfer function, solved
exactly in each layer instead of with thin layers: an independent check of SITE.

```sassi
* SOIL needs a modulus/damping curve label; with 0 iterations it is never used
DYNP,1,0.0001,1.0,0.0001,5.0,Elastic
DYNP,2,10.0,1.0,10.0,5.0,Elastic
* SPRO,<sublayer>,<L>,<label>: the same 22 sublayers, then the half-space (last row)
VAR,SAND,1,2,3,4,5,6,7,8,9,10
VAR,GRAVEL,11,12,13,14,15,16,17,18,19,20,21,22
FOREACH,SAND,SPRO,@SAND[#],1,Elastic
FOREACH,GRAVEL,SPRO,@GRAVEL[#],2,Elastic
SPRO,23,3
* SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>
SOIL,4001,32.2,1,0,0,0,0.65,1,0
* SOILX,<indir>,<mult>,<max>,<cl>: horizontal, factor 1, motion given at the top of sublayer 1
SOILX,0,1,0,1
THFILE,../data/rg160h_030g.acc
DAMP,0.05
EDUOPT,SOILCUTOFF,25
* SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>
SSAF,1,1,0,1,23,0.0976562,surface / rock outcrop
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSOIL
```

```action
plot-spectrum: ex01/SAF001W_023O.TFU
open-listing: SOIL
open-listing: SITE
```

### What this does
* `DYNP` defines a strain-dependent curve with the label `Elastic` ($G/G_{\text{max}} = 1$, 5 %
  damping). SOIL requires a label (CHECK Error 95 otherwise), but with zero iterations the curves
  are not used: the analysis is linear with the `L` properties.
* `SPRO` assigns each SOIL sublayer its `L` type, the same 22 sublayers as `TOPL`, `FOREACH` doing
  the repetition; the last `SPRO` (23) is the half-space.
* `SOIL`: 4001 values after 1 header line (the time-step line of the record), within input
  (`<outcrop>` 0), no FILE88 (`<save>` 0), **0 iterations**. `SOILX` applies the motion at the top
  of sublayer 1, the surface, as SITE does.
* `SSAF` requests the Fourier amplification of the surface motion (within) over the **outcrop**
  motion at the top of sublayer 23 (the half-space), every 0.0977 Hz. `EDUOPT,SOILCUTOFF,25`
  stops SOIL at 25 Hz to keep this short.

The result (`SAF001W_023O.TFU`, computed in a few seconds): the amplification peaks at
**2.16 at 7.4 Hz** (the first mode of the column) and 2.19 at 17.8 Hz (the second); the listing's
quarter-wavelength estimate $4H/V$, with the thickness-weighted average velocity of 1,461 ft/s, is
0.15 s, i.e. 6.6 Hz. At the 22 SITE frequencies the SOIL curve and the "surface / outcrop" line of
the SITE listing agree within **0.5 %** (for example 2.136 against 2.126 at 7.0 Hz).

### Why it matters
SITE and your site-response analysis solve the same problem with the same normalisation, so the
free field of the SSI model can be checked against the site response, and the strain-compatible
properties from SOIL can be used in SITE consistently (lesson 6). The agreement within 0.5 %
confirms both the sublayering and the half-space simulation of SITE at the frequencies of
interest. Make this check whenever you change the layering: it costs seconds.

### Technical basis
SOIL uses the exact solution of the wave equation in each homogeneous layer and the transfer of
up- and down-going wave amplitudes across interfaces (the SHAKE recursion,
[Theory §12.1](docs/theory/THEORY_MANUAL.md#12-the-shake-equivalent-linear-method-soil)); SITE
uses linear thin layers on a simulated half-space. The difference has two sources. The thin-layer
discretisation error is small here (the sublayers are 70 to 90 times thinner than the shear
wavelength at 7 Hz): compared on the motion within the column, surface over 55 ft depth, which
does not depend on the half-space model, the two agree within 0.12 %. Most of the 0.5 % comes from
the half-space simulation, which enters the **outcrop** motion; VP-02b finds 0.20 % for the same
simulation on a uniform layer (VP-02a checks SOIL, VP-02b SITE against closed forms,
[VP-02a](docs/verification/VERIFICATION_MANUAL.md#vp-02a),
[Theory §5](docs/theory/THEORY_MANUAL.md#5-half-space-variable-depth-and-viscous-boundary)).
Example 4 makes the same check with strain-compatible properties
([examples/README.md](examples/README.md#example-4-the-free-field-chain-equake---soil---site)).

## Within or outcrop: why the half-space matters

Now compare the surface motion with the **within** motion at the same depth, the top of the
half-space.

```sassi
* the same SOIL run, ratio surface (within) / top of the half-space (within)
SSAF,1,1,0,0,23,0.0976562,surface / base of the soil (within)
AFWRITE
RUNSOIL
```

```action
plot-spectrum: ex01/SAF001W_023W.TFU, ex01/SAF001W_023O.TFU
open-listing: SOIL
```

### What this does
`SSAF` for layer 1 is redefined with `<outcrop2>` = 0; `RUNSOIL` writes `SAF001W_023W.TFU` next to
the outcrop curve of the previous step. The surface-to-within ratio peaks at **17.1 at 7.3 Hz**,
almost eight times the outcrop amplification of 2.2 (SOIL listing: "Maximum amplification
surface / base WITHIN").

### Why it matters
Two lessons for design practice:

1. **Know where and how your design motion is defined.** The outcrop motion at the top of rock is
   twice the incident wave; the within motion at the same depth is incident plus reflected wave,
   and at the column frequency the two nearly cancel. A rock-outcrop motion applied in SASSI as a
   within motion at 55 ft would produce a free field many times too strong near 7 Hz. SASSI's
   control motion is always the within motion at the top of layer `<cl>`; an outcrop design motion
   is first converted with SOIL (input as outcrop, then the computed surface motion is used as the
   control motion at the surface, as example 4 shows).
2. **A half-space carries energy away; a rigid base does not.** The within ratio is exactly the
   transfer function of the 55 ft column on a rigid base: with no energy leaving through the base,
   only the 4-5 % material damping limits the resonance (17.1). With the half-space, waves radiate
   into the rock and the column resonance is 2.2 times the outcrop motion. SITE therefore
   simulates the half-space (`<nl>` = 10-20 generated sublayers plus dashpots); `<nl>` = 0 is only
   for bedrock much stiffer than the soil. The same radiation into the half-space gives the
   foundation its radiation damping (lesson 3).

### Technical basis
For vertically propagating waves the within motion at the top of the half-space is $E + F$ (incident
plus reflected amplitude) and the outcrop motion is $2E$; the column above a given depth responds to
the within motion there, whatever lies below, so $\text{surface}/\text{within}(55\,\text{ft})$ is
the rigid-base transfer function of the 55 ft column
([Theory §5](docs/theory/THEORY_MANUAL.md#5-half-space-variable-depth-and-viscous-boundary), "Within
and outcrop motions"; [User Guide §7.3](docs/user/USER_GUIDE.md#73-half-space-or-rigid-base) and
[§7.4](docs/user/USER_GUIDE.md#74-wave-field-and-control-point-site-options)).

```figure
soil-column f=7.3
The column of this lesson solved exactly, per unit motion at the surface. Near 7.3 Hz the surface
moves 2.2 times the rock outcrop but 17 times the within motion at 55 ft, the numbers of the two
SOIL runs. Raise the rock velocity towards a rigid base and the outcrop curve climbs towards the
within curve: less and less energy leaves through the base.
```

### Try this
Run SITE with a rigid base at 55 ft (`<nl>` = 0) in a copy of the model and compare its Mode 2
table with the half-space run: the amplitudes at the interfaces are identical (at 7.0 Hz the motion
at 55 ft is 0.0836 of the surface motion in both). Only the last line of each frequency changes: on
a rigid base the outcrop and the within motion at 55 ft are the same, so it reads 12.0 at 7.0 Hz
instead of 2.13. With the control motion given at the surface, the free field in the layers above
the base does not depend on what lies below. The base model matters when the motion is defined at
depth or as an outcrop motion, and for the soil springs and dashpots of lesson 3.

```sassi-show
CPMODEL,2
ACTM,2
MDL,ex01rb,../ex01_rigidbase
SITE,0,1,0,0,3,1,0,1,4096,1,0,0.005,8192,1
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
AFWRITE
RUNSITE
```

### Check yourself
The site report gives the design motion as a 0.30 g outcrop motion at the top of the rock at 55 ft
depth. What do you give SASSI as control motion, and where?

Answer: not the outcrop record at 55 ft as a within motion. Run SOIL with the record as an outcrop
motion at the top of the half-space (`<outcrop>` 1, control at sublayer 23, with the
strain-dependent curves), take the computed surface motion, and use it as the control motion at the
ground surface (`<cl>` = 1) with the strain-compatible properties in SITE. Lesson 6 does exactly
this.
