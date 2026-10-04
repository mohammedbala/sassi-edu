# SHAKE91 sample problem input (public domain)

Copied from `docs/spec/R2_benchmark_data/shake91_example/` (origin: github.com/ocrickard/SHAKE16,
public-domain modernisation of SHAKE91, `test-data/`) so that verification problem VP-04 runs
from an installed package.

- `INP.DAT`  : SHAKE91 sample input (Idriss & Sun 1992, Table B-1), ft, kcf, ft/s, g = 32.2 ft/s^2.
- `DIAM.ACC` : Loma Prieta 1989, Diamond Heights (H1_90), dt = 0.02 s, 3 header lines, 8 values
               per line, accelerations in g; 1900 values are used, scaled to 0.10 g.

VP-04 (`sassi/verify/problems/vp_soil.py`) parses `INP.DAT`, builds the SOIL deck and compares the
results with Table B-2 of the SHAKE91 manual (see docs/spec/R2_benchmarks.md section H.1).
