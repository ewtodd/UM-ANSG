# 73mGe — status and plan (2026-10-23)
<!---->
The measurement is the energy of the bare 73mGe gamma ray seen in the pixelated CZT detector (`czt/`), and in the PSI CdTe (`cdte/`, PSI CdTe below); ENSDF gives 68.752 ± 0.007 keV.
Ge-detector values (68.79–68.82 keV) include a recoil deposit this measurement does not see, so a CZT result there would mean a bias.
The previously reported 68.7502 ± 0.0145 keV (c260617) came from a pol2 that gave the Pb Kβ1,3 blend Kβ1's energy; that alone pulled Ge down by about 27 eV (Calibration, below).
<!---->
## Current result (tau20_sigpb_esc)
<!---->
| Run | Method | Ge [keV] | fit [eV] | cal [eV] | sys [eV] | transfer model [eV] | escapes [eV] | χ²/ndf | BLUE weight | Ge − mean |
|---|---|---|---|---|---|---|---|---|---|---|
| Cu shield 10%, 01/13 | in-situ | 68.7627 | 20.6 | 7.3 | — | 0 | −11.7 | 1.36 | 0.455 | +9.4 eV, +0.38σ |
| Cd shield 10%, 01/13 | in-situ | 68.7194 | 30.8 | 7.3 | — | +10.0 | −0.7 | 1.14 | 0.247 | −34.0 eV, −1.02σ |
| Cd shield 25%, 01/13 | in-situ | 68.7611 | 33.4 | 20.1 | — | +11.9 | −0.8 | 1.47 | 0.156 | +7.7 eV, +0.19σ |
| No shield 5%, 01/15 | rate-sub | 68.7662 | 42.7 | 6.3 | 30.6 | 0 | −9.7 | 1.43 | 0.085 | +12.9 eV, +0.24σ |
| Graphite castle 10%, 01/16 | rate-sub | 68.7853 | 55.9 | 9.3 | 31.3 | 0 | −9.2 | 0.95 | 0.057 | +32.0 eV, +0.49σ |
<!---->
- Quoted: **68.7533 ± 0.0206 keV**: +1.3 eV, 0.06σ from ENSDF counting its error.
  - The BLUE of all five runs, 68.7533 ± 0.0174 keV (χ²/ndf 0.39): statistics 14.0, calibration offset 4.4, calibration multiplicative 3.6, rate-sub systematics 3.2, transfer model 4.3, escapes 6.7 eV.
  - Outlier removal, as in the published budget: 11.0 eV (the BLUE without Cd 10%, the lowest run).
- Method 1 (01/13) 68.7503 ± 0.0185 keV and Method 2 (01/15–16) 68.7739 ± 0.0420 keV agree to 24 ± 46 eV.
- The uncorrelated weighted mean of the same runs is 68.7536 ± 0.0164 keV.
- The escapes column is each row's correction for the 95.902 keV line's K X-ray escapes (95.902 keV escapes, below), charged in full as a systematic: fully correlated among the runs whose escapes come from captures inside the crystal (Cu, 01/15, 01/16) and among the Cd-shield runs, independent between the two sets. Without the escapes (`tau20_sigpb`) the quoted value was 68.7630 ± 0.0191 keV; with them as a correction only, 68.7553 ± 0.0183 keV.
- The last two columns are each run's BLUE weight and its distance from the mean in its own total error. Cd 10% is the one run below the mean; its offset is the size of the signal-run anchor's shift, which Pb Kβ supports (Choices made, below). Without 01/16 the quoted value is 68.7515 ± 0.0210 keV.
<!---->
## Pipeline
<!---->
`czt/CalibrateAndFit.cpp` with `czt/GainTransferFit.hpp`; `czt/CombineGeResult.cpp` combines.
The escapes need `RooFitUtils::AddFixedLines`, on Analysis-Utilities branch `review/simultaneous-fixed-lines` (47ac109; 6ceba20 also draws the pulls of a residual's negative bins; on 5d7032b), which is not pushed: the lock still pins 5d7032b, so until it is pushed and the lock bumped, enter the shell with `nix develop .#73mGe --override-input utils 'git+file:///home/e-work/Analysis-Utilities?ref=review/simultaneous-fixed-lines'`.
<!---->
- **Method 1, 01/13, in-situ.** Each background run (Pb Kα2/Kα1) is fitted together with its signal run (the same lines plus Ge) in raw energy: one Pb shape and one low tail, each run's Pb centroids free.
  The Cu 10% pair is calibrated with a pol1 through Am-241 and its background-run Pb Kα; the Cd pairs inherit it through their background-run Pb Kα as a gain (`TRANSFER=offset`: as an offset).
  Each signal run is then shifted onto its own Pb Kα (the signal-run anchor), so only the pair's slope is inherited.
- **Method 2, 01/15–16, rate subtraction.** A pol1 through Am-241 and Ba-133 80.998 keV from the source runs (the Ba doublet with one width, one tail and ENSDF's 2.65 / 32.9 yield ratio); Ge from the live-time-normalised signal-minus-background residual, 62–71.5 keV, with the tail interpolated between Am and Ba.
  The residual is binned at the 264 eV comb period (times the day's slope) at five starting points and the fits averaged.
  The fit error is scaled by √(χ²/ndf); the systematic is the spread over starting points, the shift to a quadratic residual background, and 20 eV for carrying the source-run calibration to the beam runs, in quadrature.
- **95.902 keV escapes.** Both methods carry the line's Cd and Te K X-ray escapes as fixed lines, at the run's own 95.902 keV centroid minus the X-ray energy, with that line's Gaussian yield times the escape fractions (95.902 keV escapes, below): fractions for captures inside the crystal in the no-shield and Cu-shield runs, for photons from the Cd shield in the Cd-shield runs.
  In the in-situ fits each channel gets its own run's lines, drawn with the Pb Kα1 shape; the line's yield and centroid come from its fits in comb-period bins at five starting points over a fixed 91.5–100.5 window (in the events' units), with the in-situ fit's converged shared tail, for which the raw fit first runs once without the escapes.
  In the rate-sub fit they are a fixed term with Ge's width and tail, scaled by the residual's own 95.902 keV line, which is 1.2 times stronger in the signal runs than live time predicts.
- **01/14 is not used.** The detector was probably unbiased overnight and re-biased that morning, and there was no Am run, so no transfer from 01/13 applies. Its fit states and its testbed pair are removed.
- **Lineshape.** Gaussian plus a low exponential tail, τ/σ ≤ 20; testbed closure of the in-situ fit −1.7 ± 3.0 eV.
- **Combination**, `CombineGeResult(tag = "tau20_sigpb_esc", alt_tag = "tau20_sigpb_esc_oxfer", escape_tag = "tau20_sigpb")` (the defaults): a BLUE of every run of both methods.
  The correlated blocks are the Am reference (shared by the 01/13 rows), the gain transfer, the transfer model (each run's Ge under `alt_tag` minus under `tag`), and the escape correction (each run's Ge under `tag` minus under `escape_tag`) in two blocks, the Cd-shield runs and the rest.
  `escape_tag` must be the same configuration without the escapes; pass `""` when combining any other tag.
  Fit errors, the signal-run anchor errors and the rate-sub systematics are independent, on the diagonal.
  On top, in quadrature: the published budget's outlier removal (the BLUE without its highest run, its lowest, and both; the largest shift).
  Method 1 alone and the uncorrelated mean are printed for comparison.
<!---->
Result columns: `method ge_mu fit_err cal_err bkg_err chi2 gain_err cal_off sys_err anchor_err label`; `bkg_err` is a diagnostic, and files without the later columns read them as 0.
<!---->
Switches, read when the macro loads; each combination writes its own `root_files/calibrated_<tag>`, `plots/calibration_low_<tag>` and `results/ge_<tag>.result`:
<!---->
| Variable | Tag | Effect |
|---|---|---|
| (none) | `tau20_sigpb_esc` | production |
| `ESCAPES=off` | drops `_esc` | no 95.902 keV escapes; `tau20_sigpb` is production without them |
| `PB_ANCHOR=background` | drops `_sigpb` | signal runs keep their background run's calibration |
| `TRANSFER=offset` | `…_oxfer` | the Cd pairs inherit the reference shifted, not scaled; `tau20_sigpb_esc_oxfer` is production's transfer-model partner |
| `FIT_HARD=1` | `…_hard` | Minuit strategy 2, final HESSE, 200k calls |
| `FIDUCIAL=interior` | `…_int` | every run keeps only each crystal's 9 × 9 interior pixels |
| `BKG=linear` | `…_linbkg` | a sloped continuum in the in-situ fit |
| `RESULT_TAG_OVERRIDE=<name>` | `<name>` | plots and results land under `<name>`; the calibrated trees are still read from the switches' tag |
<!---->
Tracked: the `tau20_sigpb_esc`, `tau20_sigpb_esc_oxfer`, `tau20_sigpb`, `tau20_sigpb_oxfer` and `tau20` results, and the fit states of `tau20_sigpb_esc`, `tau20_sigpb` and `tau20` (identical files; the escapes add no fitted parameter).
A new tag needs fit states in `plots/calibration_low_<tag>/fits` copied from `tau20_sigpb_esc`'s, or the fitter opens the editor; for `BKG=linear`, free the copied `BkgSlope` (`0 0.0001 0`).
<!---->
## Configurations
<!---->
Ge − 68.752 keV per in-situ row in eV; the rate-sub rows are +23.9 (01/15) and +42.6 eV (01/16) without the escapes, +14.2 and +33.3 eV with them. Each combined alone, without a transfer-model partner.
<!---->
| Tag | Cu 01/13 | Cd 10% | Cd 25% | Method 1 | Both methods, BLUE | Outlier | Quoted |
|---|---|---|---|---|---|---|---|
| `tau20` | +26.8 | +1.7 | +6.4 | 68.7695 ± 0.0171 | 68.7715 ± 0.0158 | 3.6 | 68.7715 ± 0.0162 |
| `tau20_oxfer` | +26.8 | +12.8 | +17.8 | 68.7742 ± 0.0172 | 68.7755 ± 0.0158 | 2.1 | 68.7755 ± 0.0160 |
| `tau20_sigpb` | +22.4 | −31.9 | +9.9 | 68.7580 ± 0.0168 | 68.7615 ± 0.0153 | 12.0 | 68.7615 ± 0.0195 |
| `tau20_sigpb_oxfer` | +22.4 | −22.0 | +21.5 | 68.7627 ± 0.0165 | 68.7656 ± 0.0153 | 10.2 | 68.7656 ± 0.0184 |
| `tau20_sigpb_esc` | +10.7 | −32.6 | +9.1 | 68.7510 ± 0.0166 | 68.7542 ± 0.0154 | 10.2 | 68.7542 ± 0.0185 |
| `tau20_sigpb_esc_oxfer` | +10.7 | −22.6 | +20.9 | 68.7556 ± 0.0166 | 68.7582 ± 0.0154 | 8.5 | 68.7582 ± 0.0176 |
<!---->
The interior-pixel and sloped-continuum variants (under the earlier combination) are in the matrix of d3e02d7; with 01/14 gone the interior cut mostly widened the method gap, and the sloped continuum moved Method 1 by −1 to +6 eV.
<!---->
## PSI CdTe (cdte/)
<!---->
The same line in the PSI CdTe (Amptek, DT5780 channel 1, BOA beamline, October 2025), from our own trapezoid energies, calibrated in each run by the lines around Ge.
<!---->
- Result: **68.7448 ± 0.0242 keV** (−7.2 eV from ENSDF), the rise-time slices up to 320 ns. The collaborator's CdTe value is 68.7484 ± 0.0189 keV.
<!---->
| Term | eV | |
|---|---|---|
| statistics | 13.9 | weighted mean of five slices, 6.6 eV, × √(χ²/ndf) 2.11 of the slices ≤ 300 ns (13.4/3) |
| outlier removal | 2.1 | largest shift without the highest, the lowest, or both slices |
| lineshape | 5.1 | Ge's tail fitted on its own |
| escapes | 4.5 | Te Kα/Cd Kβ ± its error and per slice; Cd Kβ/Kα 0.324 (no-Ge runs) against Geant4's 0.275 |
| Pb Kβ blend | 1.3 | Kβ3/Kβ1 (and K-N2/K-N3) ± 10%; Kβ5, K-N4,5, K-O2,3 ± 20% |
| extrapolation | 13.9 | the scale's curvature free (Ge lies below the lines that fix the scale) |
| slice range | 1.4 | slices ≤ 300 ns instead (statistical spread 1.1 eV) |
| slice combination | 6.1 | a random-effects mean of the slices (between-slice spread 27.8 eV) |
<!---->
- **Extrapolation** is the full shift with the curvature free, +13.9 eV, not the 7.1 eV the curvature adds to the error in quadrature, because the shift is twice that. Most of it comes from the 0–270 ns slice (+165 eV, its Cd Kβ escape yield dropping from 0.19 to 0.05 of Ge, its error shrinking) and the 290–300 ns slice (+43 eV): the curvature is nearly degenerate with the escapes there, so this is conservative; charging the quadrature cost instead gives ± 0.0210 keV.
- **Pipeline** (each macro's header has the details):
  - `Trapezoid.cpp`: the raw CoMPASS waveforms in `/labdata/PSI/CdTe` (read only) to one tree per run in `/labdata/ANSG/PSI/processed/CdTe/trap`: trigger, baseline, pole-zero deconvolution (τ = 325 samples), trapezoid energies for four shapings, 10–90% rise time.
  - `CdTeAnalysis.hpp`: runs, rise slices, the EMG lineshape, the 95.902 keV line's fit, and the in-run fit: every line at its true energy on a free scale x = a + b E (Pb Kα2/Kα1, the Pb Kβ1′ and Kβ2′ groups, the 95.902 keV line's Cd and Te escapes at its measured position minus the X-ray energies, Ge free), one EMG with linear width, linear continuum, 62–89 keV, on the five Ge beam runs.
  - `GeBudget.cpp`: the slices, the budget above, `results/cdte_ge.result`, the plots in `plots/studies/psi_cdte/ge_rise_slices`.
  - `NoGeEscapes.cpp`: the escape tie from the beam runs without Ge (Te Kα/Cd Kβ 0.333 ± 0.042 at rise < 300 ns, per slice, and Cd Kβ/Kα 0.324 with Pb Kα2/Kα1 fixed), which `GeBudget.cpp` takes as constants.
  - `RiseTime.cpp` (charge loss against rise time) and `DraftComparison.cpp` (our spectra in the draft's format) draw the supporting figures.
- The import reproduces the /tmp study's numbers to the printed digits (68.7448, 13.9, 2.1, 5.1, 4.5, 10.3 eV; the escape tie). The reprocessed trapezoid trees agree with the study's to float rounding, about 1 in 10⁹.
- Line energies and the escape tables are `include/ReferenceLines.hpp`, shared with the CZT code (a production rerun on it is byte-identical).
<!---->
### CZT + CdTe
<!---->
`CombineGeResult.cpp`'s `COMBINE_WITH_CDTE` (off until the CdTe analysis is final) adds the error-weighted mean of the quoted CZT value and `results/cdte_ge.result`, taken as independent (different detectors, beam times, calibrations and escape inputs; they share only the reference line energies and Geant4's Te Kβ/Kα ratio, below an eV here). Now: **68.7498 ± 0.0157 keV** (CZT weight 0.58), the two differing by 8.5 ± 31.8 eV.
<!---->
## Findings
<!---->
### Choices made
<!---->
- **The signal-run Pb anchor.** Ge follows its own run's Pb (template fits: Ge − Pb within a signal run is the same in the Cd and Cu runs to ±45 eV while Pb moves up to 213 eV), and the Cd 10% signal run reads +36 ± 11 eV above its background run within the same crystal, pixel ring and depth slice (where the Pb lands contributes +0.2 eV): a shift of the whole run. Cu 01/13 and Cd 25% show none (+5 ± 12, +11 ± 25 eV). The cause is unknown.
  Pb Kβ, 10 keV above the window and not used by the anchor, follows: on the anchored calibration the Cd 10% signal run's Kβ reads +46 ± 30 eV against its background run's, where a shift of the Kα lines alone would leave it at −31 eV (2.6σ away) and a shift of the whole run at 0 (1.5σ). Cu 01/13 and Cd 25% read −28 ± 37 and +30 ± 77 eV, Kα itself 0 to +7 eV in all three (comb-period bins at 20 starting points; `/tmp/73mGe-escape/PbKbetaCheck.C`, plots in `plots/studies/pb_kbeta_check`). So Cd 10%'s low Ge, −34 eV from the mean and −1.0σ, is not an artefact of the anchor.
- **Gain against offset transfer is a systematic.** Pb Kβ lies 10 keV above the calibrated range, and the Cd background runs' in-run Kα1 → Kβ1 slopes split between the two models; no line inside the range resolves the 0.15% slope difference. The two move the Cd rows by +9.9 and +11.7 eV.
- **Both methods in the result.** The gap between them was not significant (0.3–0.5σ), and charging it as a systematic used Method 2 only through that noise.
<!---->
### Calibration
<!---->
- **The old pol2.** With the Kβ1,3 blend given 84.936 keV (its centroid is about 84.770 with Kβ3 at 0.52 of Kβ1), a pol2 through Am, Pb Kα and Kβ puts Cu 01/13's Ge 27.4 eV below the pol1; with the blend's centroid, 5.4 eV above.
- **Curvature, 01/13.** The in-run Kα1 → Kβ1 slope over 75–85 keV (eight runs, 0.98817 ± 0.00052) matches Am → Kα over 59.5–75 keV (0.98822): a quadratic of that size moves Ge by about 1 eV. Both lines lie outside the anchors, so this bounds a curvature across the range, not structure between Am and Pb.
- **Curvature, 01/15–16.** Ba-133 53.16 keV reads −5.6 ± 16.7 eV on the 01/15 calibration; a pol2 through Ba 53 + Am + Ba 81 moves Ge −3.6 (01/15) and +17 eV (01/16, Ba 53 at ±52 eV). The response is strongly nonlinear above 81 keV (the low pol1 misses Ba 276–384 keV by 1.8–3.9 keV), so those lines cannot constrain curvature near Ge.
- **Carrying the source-run calibration to the beam runs (01/15–16)** cannot be checked in those runs: no line appears in both (the source runs show only their own lines and Cs-137, the beam runs' other lines are prompt captures), and Pb Kα in the beam runs is a few-percent bump under the 95.9 keV line's Cd-Kα escape (72.73 / 72.92 keV) that does not fit stably. On 01/13 the same step, Am source run to beam run, holds to about 20 eV (the Kβ slope above); that is the systematic used.
- **Run-to-run shifts.** The scale moves between runs by tens of eV, up to 130 eV at the 01/13 Cd → Cu shield change. From the two 01/13 test runs Pb Kα read −7.5 eV per 1000 counts/s per crystal, but the signal runs read 12–57 eV above that line and 01/15 and 01/16 agree to −3 ± 66 eV where it predicts 200 eV, so it is not used as a model. Count rates (filtered, per crystal): 01/13 pairs 1.0–12k; 01/15 beam 5.8–7.5k against Am 86 and Ba 1393; 01/16 castle 26–34k against Am + Ba 254.
- **Depth.** Lines read higher deeper in the crystal, and not as a pure offset: across four depth slices Am moves 220 eV and Pb Kα1 130 eV (slope 0.984–0.992). A per-depth calibration is statistically possible for Cu and Cd 10% (marginal for Cd 25%) and, in throwaway fits, leaves Method 1 near 68.763 keV. Pb measured on Ge's own depth and pixel cells reads 10.6 ± 5.0 eV lower than on Pb's.
- **Edge pixels.** In each crystal's outer ring and its ring along the gap, Pb Kα1 reads 7–155 eV off its all-pixel value.
- **The ~43 keV "peak"** is the threshold: split by each pixel's turn-on (40–45 keV), it moves with it (42.14 / 42.94 / 43.75 keV for turn-ons at 40.13 / 41.12 / 42.19). Not a line.
<!---->
### Energy grids and binning
<!---->
- Every pixel's energies sit on their own discrete 1–8 eV grid; summed, they show a comb with a 264 eV period (raw, on every day, aligned between a day's runs).
- Bins near the comb period beat against it: 255 eV bins move rate-sub Ge by 13–16 eV rms with where they start, comb-period bins by 8–10 eV.
- The signal and background runs fill the grids differently, so signal − α × background is overdispersed: 7–16 times Poisson at 0.5–10 eV bins, about 4 at 50 eV, 1.4–1.9 at 255 eV, with an additive or a multiplicative smooth residual alike. Unbinned or fine-binned two-run fits therefore fail (a parametric two-run fit and an on/off likelihood with the background profiled per bin both did), and from ~20 eV bins up the on/off likelihood agrees with the subtraction.
- The unbinned in-situ fits are not affected: dithering every energy by ±132 eV (a boxcar one comb period wide) moves the 01/13 rows by −8 to +4 eV against 2.5–3 eV of dither noise, Method 1 by about −3 eV, and Am by < 1 eV. It does lower the binned χ²/ndf the fits print (Cu 01/13 signal 1.73 → 0.96–1.56): part of it is the comb in the 255 eV display bins.
<!---->
### Lineshape and fit model
<!---->
- **Pb Kα2/Kα1.** ENSDF gives 0.595; the background runs fit 0.595 ± 0.015 (Cu), 0.567 ± 0.013 (Cd 10%) and 0.577 ± 0.013 (Cd 25%), as Pb and Cd absorb 72.8 keV more than 75.0 keV. Fixing it moves both lines together by ≤ 6 eV and Method 1 by about −1.5 eV; not adopted.
- **High-side tail (pileup).** Background runs never take one; shared by a pair it is zero in three of four; per run, the signal runs take a 3–5% tail at the τ/σ = 1 bound (a high-side widening), Cd 25% rails at τ/σ = 20 (a pedestal), Cd 10% fails. Where sensible it moves Ge up 5–8 eV. Not adopted.
- **Signal-run-only in-situ fit** (a one-channel `FitSimultaneous`): the signal runs' own tails come out shorter, moving Ge −3 (Cu) and −16 eV (Cd 10%); Cd 25%'s signal run is too weak to fix its shape alone.
- **Ba doublet.** With free yields the weak line takes 0.0866 (01/15) and 0.1066 (01/16) of the strong one against 0.0805; fixing it leaves χ²/ndf unchanged.
<!---->
### Method 2 checks
<!---->
- **Without subtraction** the fits read +91 (01/15) and +134 eV (01/16), and a template with a free scale degenerates (it fits 0.3–1.3 against live-time ratios of 3.6 and 1.8): the subtraction removes real structure near Ge. The 95.9 keV line, whose escapes fall around Ge, is 1.2 times stronger in the signal runs than live time predicts; the escapes this leaves in the residual are now fitted.
- **The isomer's 53.440(9) keV line** is not cleanly visible in the residuals (a shelf at 51–54 keV); no use as an in-run anchor.
- **The 01/16 Ge-on-CZT run** (0.5% power, sample on the detector) has a Ge line about 20 times weaker than the castle run's; rate-sub fits give ±240–310 eV. No use.
<!---->
### 95.902 keV escapes
<!---->
- **The line** is 113Cd(n,γ)'s E2 transition in 114Cd, converted in the K shell 47% of the time and emitted as a photon 37% (BrIcc: αK 1.288, αT 1.735). Its escapes come from the Cd K X-rays of conversions at captures near a surface, and from Cd and Te K X-rays when its photon is absorbed: Cd Kα at 72.73 / 72.92 keV under Pb Kα2, Cd Kβ at 69.25–69.81 keV beside Ge, Te Kα at 68.43 / 68.70 keV under Ge, Te Kβ at 64.2–64.9 keV.
- **Captures inside the crystal** (no shield, Cu shield), as fractions of the line's Gaussian yield: Cd Kβ 0.0119, the mean of the 01/15 and 01/16 beam background runs (0.0122, 0.0115, each ± ~0.0017, the one escape clear of other lines); Cd Kα 0.0450 from it through the DDEP intensities times escape probabilities; Te Kα 0.0024 from Geant4's 0.053 of Cd Kα (the runs give 0.0006 and 0.0021, ± 0.002–0.004); Te Kβ 0.0007. Geant4 alone (CZT, captures at 0.25 mm mean depth) gives 0.061, 0.017, 0.0033 and 0.0010.
- **Photons from the Cd shield**, which stops the thermal neutrons: Geant4 with the line's photons onto the face, 0.0115, 0.0033, 0.0057 and 0.0017. The 01/13 Cd background runs show no Cd Kβ escape (0 ± 0.004 at 10%, ± 0.006 at 25%), where captures inside would give 0.012; the Cu run cannot tell (± 0.008–0.026).
- **Where they sit**: at the run's own 95.902 keV centroid minus the X-ray energy. In the 01/16 background run the line and its escapes read low together (−244 eV).
- **Effect**: Cu 01/13 −11.7 eV, the Cd pairs −0.7 and −0.8 eV, the rate-sub rows −9.7 and −9.2 eV; the quoted value −7.7 eV. A TF1 replica of the in-situ fit with the same lines gives −12.7, +0.1 and −1.2 eV.
- **Sensitivity**: with Te Kα (and Kβ) doubled, as PSI's CdTe runs without Ge suggest (Te Kα / Cd Kα about twice Geant4's for captures inside), the rows move +2.4 (Cu) and +1.4 / +3.1 eV (rate-sub), well inside the escape systematic, which is the whole correction.
- **The fixed lines** (Analysis-Utilities `AddFixedLines`) reproduce a testbed's doublet without bias (−0.1 ± 2.1 and −0.3 ± 1.4 eV, three known lines among them; +12.4 and +7.6 eV without them), and with the escapes switched off production reproduces `tau20_sigpb` exactly.
<!---->
### Ruled out
<!---->
- **The 95.9 keV line as a probe**: in Cu runs it reads 31–67 eV low against Pb, and −270 eV in the 01/16 castle runs.
<!---->
### Removed (2026-10-23)
<!---->
Code that no longer fed the result, cut so it cannot be mistaken for part of it:
<!---->
- **`LINESHAPE=legacy`** (`tau8`, τ/σ ≤ 8): its tail bound railed (2026-10-22 lineshape tests), so that lineshape was set by the bound, not fitted; τ/σ ≤ 20 replaced it.
- **The earlier combination scheme** in `CombineGeResult` (Method 1 with the gap to Method 2 as a systematic, and its split of the error without calibration): superseded by the BLUE of both methods (Choices made).
- **The calibration-peak-selection term** (`peak_sel_eV`): the quoted result took the default, 0.
- **`RateSubAlpha.cpp`**, a subtraction-scale scan and Gaussian-against-tail test on Gaussian-fitted 255 eV residuals, and **`AdditionalLevels.cpp`**, a likelihood-ratio search for other sample lines (73Ge, 77Ge, 77As) in the beam runs: one-off studies on the earlier rate-sub fit and calibration, superseded by the shared-lineshape, comb-period fit.
- **01/14**: its fit states in `calibration_low_tau20` and its pair in `InjectionTest.cpp`.
- **Diagnostics in `CalibrateAndFit`**: the Ge/Pb rate-ratio printout, the Am-plus-Kα variant of the validation (the Pb-only validation and the linearity check stay), the unit check and the pre-cal Pb comparison table. In `GainTransferFit.hpp`: the in-situ path with one Pb line and the rate-sub fit's free-tail model, unreachable since the two-line model and the shared lineshape; the single-histogram rate binning gave way to the five-phase binning, which `CalibrateAndFit` and `InjectionTest` now share with the phase-averaged rate-sub fit.
- **`Filter.cpp`'s `GetTrimmedPixelCenters`**: never called.
<!---->
### Analysis-Utilities
<!---->
- `FitResult::cov_qual`'s comment in `include/FittingUtils.hpp` has Minuit's scale inverted (since 27d7022); the code reads it correctly.
- The retry recovers the 01/15 Ba refit (status 0, χ²/ndf 1.488).
<!---->
## Plan
<!---->
1. **Method 1's remaining model choices** are not in the budget: the in-situ background shape (flat against sloped, −1 to +6 eV on Method 1), the lineshape (shared tail, τ/σ ≤ 20) and the pixel selection. Decide which become systematics, as the transfer model did.
2. **The Cd pairs** inherit the Cu reference across the 130 eV step at the Cd → Cu change; if the bias was cycled then, they have 01/14's problem.
3. **01/16**: kept, with its larger errors (BLUE weight 0.057; without it the quoted value moves −1.8 eV to 68.7515 ± 0.0210 keV). Its residual runs from −5 counts/s at 62 keV to 0 at 71 keV, the signal run below the live-time-scaled background, and the signal/background rate ratio is below 1 over 55–70 keV; it is the highest-rate run (26–34k counts/s per crystal) and has no line in both its runs to check the source-run calibration. Decide whether it stays.
4. **Analysis-Utilities**: report the `cov_qual` comment; review the five checklists and `review/simultaneous-fixed-lines`, then push it and bump the lock (`nix flake lock --update-input utils`), after which the 73mGe shell needs no override.
5. **Injection sweep**: `InjectionTest.cpp` still requires the `gaintransfer` tag and lists the 01/14 pair.
6. **Write-up**: the move from 68.7502 keV step by step — the Kβ blend in the old pol2 (−27 eV), the RooFit fix (≤ 4 eV), fit convergence (up to 36 eV), the lineshape, the Pb centroid link, dropping 01/14, the signal-run anchor, the transfer model, and both methods in the combination.
