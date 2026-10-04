# RG 1.60 targets for EQUAKE

Generated with `sassi.core.equake_lib` (`rg160_spectrum`, `rg160_target_psd`); see
docs/spec/R2_benchmarks.md sections G.2 and G.3.

- `RG160H_5pct_1g.rsi`, `RG160V_5pct_1g.rsi`: RSIN target spectra (Hz, SA in g), 5 % damping,
  anchored to 1.0 g, 27 frequencies from 0.1 to 100 Hz including every RG 1.60 control point
  (0.25, 2.5, 3.5, 9, 33 Hz), so EQUAKE's log-log interpolation reproduces the RG 1.60 lines
  exactly. Use `EQUAKE,0,27,...` (Number of Frequencies = 27). Scale linearly for other PGA.
- `RG160H_AppA_1g_cm2s3.tpsd` / `RG160H_AppA_1g_in2s3.tpsd`: SRP 3.7.1 Rev. 4 Appendix A minimum
  PSD for the RG 1.60 horizontal spectrum at 1 g (TPSD files for SI and British gravity units).
  Scale by PGA^2. No vertical target PSD is defined in SRP 3.7.1 Appendix A.
