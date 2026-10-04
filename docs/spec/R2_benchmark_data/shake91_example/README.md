# SHAKE91 sample problem (public domain)

Source: https://github.com/ocrickard/SHAKE16 (public-domain modernisation of SHAKE91; `test-data/`).
The SHAKE16 test suite reports byte-identical output to the original SHAKE91 executable.

- `INP.DAT`  : SHAKE91 sample input (Idriss & Sun 1992, User's Manual Table B-1), English units
               (ft, kcf, fps; g = 32.2 ft/s^2 for mass, 981 cm/s^2 for spectra).
- `DIAM.ACC` : Loma Prieta 1989, Diamond Heights (H1_90), dt = 0.02 s, 2000 pts (1900 read), g units,
               PGA 0.1129 g, scaled to 0.10 g by SHAKE91; applied as OUTCROP motion at layer 17 (half-space).
- `output1_SHAKE16.txt` : output file #1 (strain-compatible properties, PGA profile, spectra, amplification).

Published reference values (SHAKE91 manual, Table B-2) are tabulated in `../../R2_benchmarks.md`, section H.
