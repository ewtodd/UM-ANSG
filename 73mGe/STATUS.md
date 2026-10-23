# 73mGe — status and plan (2026-10-23)
<!---->
The measurement is the energy of the bare 73mGe gamma ray seen in the pixelated CZT detector; ENSDF gives 68.752 ± 0.007 keV.
Ge-detector values (68.79–68.82 keV) include a recoil deposit this measurement does not see, so a CZT result there would mean a bias.
The previously reported 68.7502 ± 0.0145 keV (c260617) came from a pol2 that gave the Pb Kβ1,3 blend Kβ1's energy; that alone pulled Ge down by about 27 eV (Calibration, below).
<!---->
## Current result (tau20, 01/14 dropped)
<!---->
| Run | Method | Ge [keV] | fit [eV] | cal [eV] | χ²/ndf |
|---|---|---|---|---|---|
| Cu shield 10%, 01/13 | in-situ | 68.7788 | 20.7 | 6.3 | 1.72 |
| Cd shield 10%, 01/13 | in-situ | 68.7537 | 30.7 | 12.9 | 1.08 |
| Cd shield 25%, 01/13 | in-situ | 68.7584 | 33.0 | 12.1 | 1.48 |
| No shield 5%, 01/15 | rate-sub | 68.7858 | 35.1 | 6.3 | 1.55 |
| Graphite castle 10%, 01/16 | rate-sub | 68.7830 | 56.2 | 9.3 | 0.91 |
<!---->
- Method 1 (BLUE, three in-situ runs): 68.7695 ± 0.0171 keV, χ²/ndf 0.26.
- Method 2 (01/15–16): 68.7850 ± 0.0302 keV.
- Quoted: **68.7695 ± 0.0231 keV**, the 15.5 eV method gap as a systematic: +17.5 eV, 0.7σ from ENSDF.
- Paper-style: 68.7723 ± 0.0148 keV.
<!---->
## Pipeline
<!---->
`macros/CalibrateAndFit.cpp` with `macros/GainTransferFit.hpp`, on Analysis-Utilities 5d7032b.
<!---->
- **Method 1, 01/13, in-situ.** Each background run (Pb Kα2/Kα1) is fitted together with its signal run (the same lines plus Ge) in raw energy: one Pb shape and one low tail, each run's Pb centroids free.
  The Cu 10% pair is calibrated with a pol1 through Am-241 and its background-run Pb Kα; the two Cd pairs inherit it through their own Pb Kα.
- **Method 2, 01/15–16, rate subtraction.** A pol1 through Am-241 and Ba-133 80.998 keV from the source runs; Ge from the live-time-normalised signal-minus-background residual (62–71.5 keV, 255 eV bins, linear residual background, tail interpolated between Am and Ba).
  The Ba doublet is fitted with one width, one tail and ENSDF's 2.65 / 32.9 yield ratio.
- **01/14 is not used.** The detector was probably unbiased overnight and re-biased that morning, and there was no Am run, so no transfer from 01/13 applies. It read +40 to +99 eV in every configuration. Its three fit states are still tracked, unused.
- **Lineshape.** Gaussian plus a low exponential tail, τ/σ ≤ 20; testbed closure of the in-situ fit −1.7 ± 3.0 eV.
- **Combiners.** `CombineGeResult(tag)`: BLUE of Method 1, |Method 1 − Method 2| as the method systematic. `CombineGeResultAlt(tag)`: paper-style mean with outlier removal.
<!---->
Switches, read when the macro loads; each combination writes its own `root_files/calibrated_<tag>`, `plots/calibration_low_<tag>` and `results/ge_<tag>.result`:
<!---->
| Variable | Tag part | Effect |
|---|---|---|
| (none) | `tau20` | production |
| `LINESHAPE=legacy` | `tau8` | lineshape before 2026-10-22 |
| `FIT_HARD=1` | `_hard` | Minuit strategy 2, final HESSE, 200k calls |
| `PB_ANCHOR=signal` | `_sigpb` | each signal run calibrated on its own Pb Kα (the pair calibration plus an offset) |
| `TRANSFER=offset` | `_oxfer` | the Cd pairs inherit the reference calibration shifted, not scaled |
| `FIDUCIAL=interior` | `_int` | every run keeps only each crystal's 9 × 9 interior pixels |
| `BKG=linear` | `_linbkg` | a sloped continuum in the in-situ fit |
<!---->
A new tag needs fit states in `plots/calibration_low_<tag>/fits` copied from `tau20`'s, or the fitter opens the editor; for `BKG=linear`, free the copied `BkgSlope` (`0 0.0001 0`).
Only `tau20`'s result and states are tracked.
<!---->
## Configuration matrix (01/14 dropped)
<!---->
Ge − 68.752 keV in eV per in-situ row; σ counts the ENSDF error. The rate-sub rows are +34 / +31 eV with all pixels and +76 / −10 eV with interior pixels.
<!---->
| Tag | Cu 01/13 | Cd 10% | Cd 25% | Method 1 (χ²/ndf) | Gap | Quoted | Paper-style |
|---|---|---|---|---|---|---|---|
| `tau20` | +27 | +2 | +6 | 68.7695 ± 0.0171 (0.26) | 15 | 68.7695 ± 0.0231 (0.7σ) | 68.7723 ± 0.0148 |
| `tau20_oxfer` | +27 | +13 | +18 | 68.7742 ± 0.0172 (0.07) | 11 | 68.7742 ± 0.0203 (1.0σ) | 68.7762 ± 0.0145 |
| `tau20_sigpb` | +22 | −32 | +10 | 68.7582 ± 0.0184 (1.09) | 27 | 68.7582 ± 0.0325 (0.2σ) | 68.7639 ± 0.0182 |
| `tau20_sigpb_oxfer` | +22 | −22 | +22 | 68.7625 ± 0.0177 (0.77) | 22 | 68.7625 ± 0.0286 (0.4σ) | 68.7676 ± 0.0172 |
| `tau20_int` | +38 | −16 | +22 | 68.7762 ± 0.0185 (0.89) | 25 | 68.7762 ± 0.0312 (0.8σ) | 68.7805 ± 0.0181 |
| `tau20_oxfer_int` | +38 | −8 | +33 | 68.7802 ± 0.0186 (0.67) | 21 | 68.7802 ± 0.0281 (1.0σ) | 68.7839 ± 0.0172 |
| `tau20_sigpb_int` | +36 | −44 | +4 | 68.7632 ± 0.0284 (2.16) | 38 | 68.7632 ± 0.0475 (0.2σ) | 68.7705 ± 0.0226 |
| `tau20_sigpb_oxfer_int` | +36 | −33 | +17 | 68.7675 ± 0.0244 (1.59) | 34 | 68.7675 ± 0.0417 (0.4σ) | 68.7744 ± 0.0213 |
| `tau20_linbkg` | +29 | +5 | −14 | 68.7684 ± 0.0170 (0.60) | 17 | 68.7684 ± 0.0238 (0.7σ) | 68.7709 ± 0.0154 |
| `tau20_sigpb_oxfer_int_linbkg` | +44 | −33 | +24 | 68.7736 ± 0.0270 (1.80) | 28 | 68.7736 ± 0.0387 (0.6σ) | 68.7793 ± 0.0220 |
<!---->
Every quoted value is within 1σ of ENSDF (+6 to +28 eV). Only the `tau20`, `tau20_oxfer`, `tau20_sigpb` and `tau20_sigpb_oxfer` plot directories are kept; the others are rebuilt when needed.
<!---->
## Open choices
<!---->
- **Pb anchor.** `PB_ANCHOR=signal` moves essentially only Cd 10% (−34 eV).
  The a priori case: in the same crystal, pixel ring and depth slice, that signal run's Pb reads +36 ± 11 eV above its background run's, so the shift is run-wide (where the Pb lands contributes +0.2 eV), and Ge in the same run should share it. Cu 01/13 and Cd 25% show no such shift (+5 ± 12, +11 ± 25 eV).
  The cause is unknown; the rate law predicts the opposite sign.
  It cannot be chosen because it lands closest to ENSDF.
- **Transfer.** Gain against offset cannot be decided inside the calibrated range: Pb Kβ lies 10 keV above it, and the Cd background runs' in-run Kα1 → Kβ1 slopes split between the two. It moves Method 1 by 4–5 eV and the Cd rows by 8–13 eV; a systematic.
- **Interior pixels.** Found while chasing 01/14; with 01/14 gone it mostly widens the method gap.
- **Method 2's binning and residual background.** The summed spectra carry a 264 eV comb (raw energy, every day, aligned between signal and background runs to 1 eV), and production's 255 eV bins beat against it: shifting the bins moves Ge by 13–16 eV rms, 8–10 eV with bins of the comb period, which also lower 01/15 by 15 eV. A quadratic residual background instead of the linear one moves Ge −19 (01/15) and −24 eV (01/16) with no χ² preference. Neither effect is in the error budget.
<!---->
## Findings
<!---->
### Calibration
<!---->
- **The old pol2.** With the Kβ1,3 blend given 84.936 keV (its centroid is about 84.770 with Kβ3 at 0.52 of Kβ1), a pol2 through Am, Pb Kα and Kβ puts Cu 01/13's Ge 27.4 eV below the pol1; with the blend's centroid, 5.4 eV above.
- **Curvature, 01/13.** The in-run Kα1 → Kβ1 slope over 75–85 keV (eight runs, 0.98817 ± 0.00052) matches Am → Kα over 59.5–75 keV (0.98822): a quadratic of that size moves Ge by about 1 eV. Both lines lie outside the anchors, so this bounds a curvature across the range, not structure between Am and Pb.
- **Curvature, 01/15–16.** Ba-133 53.16 keV reads −5.6 ± 16.7 eV on the 01/15 calibration; a pol2 through Ba 53 + Am + Ba 81 moves Ge −3.6 (01/15) and +17 eV (01/16, Ba 53 at ±52 eV). The 81 keV line's Te-Kα escape (53.53 / 53.80 keV) can pull Ba 53. The response is strongly nonlinear above 81 keV: the low pol1 misses Ba 276–384 keV by 1.8–3.9 keV, a pol2 cannot span 53–384 keV, and a pol3 through everything moves Ge −10 and −13 eV; those lines do not usefully constrain curvature near Ge.
- **No Pb in the source runs**, so 01/15–16 have no anchor between Am and Ba 81.
- **Depth.** Lines read higher deeper in the crystal, and not as a pure offset: across four depth slices Am moves 220 eV and Pb Kα1 130 eV, so the slope varies 0.984–0.992. A per-depth calibration is statistically possible for Cu and Cd 10% (marginal for Cd 25%, whose signal-run Pb is 6–11k per slice) and, in throwaway TF1 fits, leaves Method 1 near 68.763 keV. Ge sits in slightly different depth cells than Pb; Pb measured on Ge's cells reads 10.6 ± 5.0 eV lower.
- **Edge pixels.** In each crystal's outer ring and its ring along the gap, Pb Kα1 reads 7–155 eV off its all-pixel value.
- **The ~43 keV "peak"** is the threshold: split by each pixel's turn-on (40–45 keV), it moves with it (42.14 / 42.94 / 43.75 keV for turn-ons at 40.13 / 41.12 / 42.19). Not a line, not an anchor.
<!---->
### Lineshape and fit model
<!---->
- **Pb Kα2/Kα1.** ENSDF gives 0.595; the background runs fit 0.595 ± 0.015 (Cu), 0.567 ± 0.013 (Cd 10%) and 0.577 ± 0.013 (Cd 25%), as Pb and Cd absorb 72.8 keV more than 75.0 keV. Fixing it moves both lines together by ≤ 6 eV and Method 1 by about −1.5 eV; not adopted.
- **High-side tail (pileup).** Background runs never take one; shared by a pair it is zero in three of four; per run, the signal runs take a 3–5% tail at the τ/σ = 1 bound (a high-side widening), Cd 25% rails at τ/σ = 20 (a pedestal) and Cd 10% fails. Where sensible it moves Ge up 5–8 eV. Not adopted.
- **Signal-run-only in-situ fit** (a one-channel `FitSimultaneous`, no library change): the signal runs' own tails come out shorter, moving Ge −3 (Cu) and −16 eV (Cd 10%); Cd 25%'s signal run is too weak to fix its shape alone.
- **Ba doublet.** With free yields the weak line takes 0.0866 (01/15) and 0.1066 (01/16) of the strong one against 0.0805; fixing it leaves χ²/ndf unchanged.
<!---->
### Runs and rates
<!---->
- 01/13 start times: Test 5% 07:04, Cd bkg 25% 07:32, Test 90% 08:09, Cd sig 25% 08:32, Cd sig 10% 08:52, Cd bkg 10% 10:04, Cu sig 10% 11:41, Cu bkg 10% 13:45, Am-241 15:43.
- From the two 01/13 test runs Pb Kα reads −7.5 eV per 1000 counts/s per crystal; corrected for that it is flat to 10:04 and 130–137 eV higher from 11:41, at the Cd → Cu shield change.
- The rate law does not carry over to signal runs: at 1.7–2.5 times their background runs' rates (filtered 1.0–12k counts/s per crystal) they read 12–57 eV above it.
- Rates against the calibration runs: 01/13 Am 33 counts/s per crystal; 01/15 Am 86, Ba 1393, beam 5.8–7.5k; 01/16 Am + Ba 254, castle 26.5–33.5k.
- In template fits Ge − Pb within a signal run is the same in the Cd and Cu runs (±45 eV) while Pb moves up to 213 eV between them.
<!---->
### Method 2 checks
<!---->
- **Fits without subtraction** read +91 (01/15) and +134 eV (01/16) on comb-period bins, and a template with a free scale degenerates (it fits 0.3–1.3 against live-time ratios of 3.6 and 1.8): the subtraction removes real structure near Ge. Pb is not visible in these beam runs, and the 95.9 keV line, whose escapes fall under Ge, is 1.2 times stronger in the signal runs than live time predicts.
- **The isomer's 53.440(9) keV line** is not cleanly visible in the residuals (a shelf at 51–54 keV); no use as an in-run anchor.
- **Unbinned fits are not viable.** Every pixel's energies sit on their own discrete 1–8 eV grid (the 264 eV comb is the coarse imprint), and the signal and background runs populate those grids differently. Signal − α × background on the continuum below Ge is overdispersed 12–16 (01/15) and 7–12 times (01/16) at 0.5–10 eV bins, about 4 at 50 eV and 1.4–1.7 at 255 eV, with an additive or a multiplicative smooth residual alike. A parametric two-run fit and an on/off likelihood with the background profiled per bin both fail below ~20 eV; from there up the on/off likelihood (per-crystal live time and residual) agrees with rate-sub: 68.76–68.81 keV on 01/15 and 68.77–68.80 keV on 01/16 at 255 and 262 eV bins.
- **Rate-sub errors** are not scaled for the residual's overdispersion (χ²/ndf 1.4–1.9 at 255 eV): they are about 20–35% too small.
<!---->
### Ruled out
<!---->
- **95.9 keV escapes under Ge in the in-situ runs**: at the measured 0.5% (Te) and 2.5% (Cd) of the line, Ge moves ≤ 7 eV, and unmodelled escapes pull it low.
- **The 95.9 keV line as a probe**: in Cu runs it reads 31–67 eV low against Pb, and −270 eV in the 01/16 castle runs.
<!---->
### Analysis-Utilities
<!---->
- `FitResult::cov_qual`'s comment in `include/FittingUtils.hpp` has Minuit's scale inverted (since 27d7022); the code reads it correctly.
- The retry recovers the 01/15 Ba refit (status 0, χ²/ndf 1.488).
<!---->
## Plan
<!---->
1. **Configuration**, on grounds that do not use the expected energy: Pb anchor (the run-wide Cd 10% shift argues for the signal run), transfer model (a systematic), interior cut. Make the choice the default and track its states and result.
2. **The Cd pairs** inherit the Cu reference across the 130 eV step at the Cd → Cu change; if the bias was cycled then, they have 01/14's problem.
3. **Method 2**: settle the binning (comb period or the on/off fit) and carry the residual-background shape as a systematic; decide on 01/16.
4. **Per-depth calibration** with production-quality fits, if wanted.
5. **Analysis-Utilities**: report the `cov_qual` comment; review the five checklists.
6. **Injection sweep**: `InjectionTest.cpp` still requires the `gaintransfer` tag and lists the 01/14 pair.
7. **Write-up**: the move from 68.7502 keV step by step — the Kβ blend in the old pol2 (−27 eV), the RooFit fix (≤ 4 eV), fit convergence (up to 36 eV), the lineshape, the Pb centroid link, dropping 01/14, and the configuration choices.
