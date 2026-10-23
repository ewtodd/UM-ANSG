# 73mGe — status and plan (2026-10-23)
<!---->
The measurement is the energy of the bare 73mGe isomeric gamma ray seen in the pixelated CZT detector.
ENSDF gives 68.752 ± 0.007 keV for the gamma ray.
The Ge-detector values (68.79–68.82 keV) include the recoil deposit that a CZT measurement of an external source does not see, so a CZT result there would mean a bias.
The previously reported 68.7502 ± 0.0145 keV (c260617) came from a pol2 calibration that assigned 84.936 keV to the Pb Kβ1,3 blend; the change from it has to be explained step by step, not as one bug.
<!---->
## Production pipeline
<!---->
`macros/CalibrateAndFit.cpp` with `macros/GainTransferFit.hpp` (replaces `CalibrationLow.cpp`):
<!---->
- **Method 1, 01/13–14, in-situ.**
  Simultaneous fit of each background run (Pb Kα2/Kα1) and its signal run (the same Pb lines plus Ge) in raw energy.
  The reference pair (Cu 10%, 01/13) is calibrated with a pol1 through Am-241 (source run) and the background-run Pb Kα lines; the other pairs inherit it through a Pb gain factor.
  Ge is read from the post-calibration fit of the signal run.
- **Method 2, 01/15–16, rate subtraction.**
  Pol1 through Am-241 and Ba-133 80.998 keV from the post-reactor source runs; Ge from the live-time-normalised signal minus background residual, with the tail interpolated between Am and Ba.
  The 01/14 rate-subtracted row has been dropped: 01/14 is either fixed in the in-situ fit or dropped entirely.
- **Lineshape.** Gaussian plus a low-side exponential tail, τ/σ ≤ 20, the tail amplitude a full-range yield ratio (Analysis-Utilities 7dc3ae8).
  Testbed closure of the in-situ fit with physical truth: −1.7 ± 3.0 eV over 50 testbeds.
- **Combiners.** `CombineGeResult(tag)` quotes Method 1 by BLUE and takes |Method 1 − Method 2| as the method systematic; `CombineGeResultAlt(tag)` is the paper-style mean with outlier removal.
<!---->
Switches, fixed when the macro loads; each combination writes its own `root_files/calibrated_<tag>`, `plots/calibration_low_<tag>` and `results/ge_<tag>.result`:
<!---->
| Variable | Tag part | Effect |
|---|---|---|
| (none) | `tau20` | production |
| `LINESHAPE=legacy` | `tau8` | lineshape before 2026-10-22 (τ/σ ≤ 8, Gaussian rate-sub) |
| `FIT_HARD=1` | `_hard` | Minuit strategy 2, final HESSE, 200k calls |
| `PB_ANCHOR=signal` | `_sigpb` | each signal run calibrated on its own Pb Kα (see below) |
<!---->
A new tag needs saved fit states in `plots/calibration_low_<tag>/fits`, copied from `tau20`'s; without them the fitter opens the interactive editor.
Only the `tau20` result and its fit states are tracked, besides `results/ge_amxfer.result` from the removed amxfer configuration.
<!---->
## Current result (tau20)
<!---->
| Run | Method | Ge [keV] | fit [eV] | cal [eV] | χ²/ndf |
|---|---|---|---|---|---|
| Cu shield 10%, 01/13 | in-situ | 68.7760 | 20.1 | 4.1 | 1.58 |
| Cd shield 10%, 01/13 | in-situ | 68.7291 | 30.4 | 8.1 | 1.14 |
| Cd shield 25%, 01/13 | in-situ | 68.7596 | 33.3 | 8.5 | 1.49 |
| Cu shield 10%, 01/14 | in-situ | 68.8375 | 19.7 | 7.4 | 2.34 |
| No shield 5%, 01/15 | rate-sub | 68.7759 | 34.4 | 6.3 | 1.57 |
| Graphite castle 10%, 01/16 | rate-sub | 68.7877 | 57.6 | 11.3 | 0.85 |
<!---->
- Method 1 (BLUE, the four in-situ runs): 68.7887 ± 0.0256 keV, χ²/ndf 3.71.
- Method 2 (01/15–16): 68.7790 ± 0.0301 keV.
- Quoted: **68.7887 ± 0.0274 keV**, with the 9.7 eV method gap as a systematic.
- Without 01/14 (`results/ge_tau20_no0114.result`, untracked): Method 1 68.7621 ± 0.0158 keV, quoted 68.7621 ± 0.0232 keV with a 16.9 eV method gap.
<!---->
None of these is final: 01/14 and 01/16 are in doubt (below), and the in-situ fit has had the signal and background Pb centroids tied.
<!---->
## Findings, 2026-10-22/23
<!---->
### The detector reading moves with count rate
<!---->
- Within one data-taking session, Pb Kα reads 8–11 eV lower per 1000 counts/s per crystal: −77 eV from 10% to 90% power on 01/14 (2.4k → 9.4k /s), −140 eV from the 5% to the 90% test run on 01/13 (1.2k → 19.8k /s).
- Pb Kα (75 keV) and Pb Kβ (85 keV) move together, so over 75–85 keV the shift is an offset, not a gain change; an offset reaches the Ge line at full size.
- The 01/13 Cd-shield runs, taken before the Cu-shield runs, sit 110–140 eV below them at the same rate: a drift over time, the shield change, or both.
- At fixed conditions a run is flat to ±20 eV (01/14 Cu 10%, 3.9 h); the 90% run settles by about 20 eV after its first quarter.
- The Am-241 and Ba-133 source runs ran at 33–1440 counts/s per crystal, the beam runs at 1k–33k.
<!---->
### The 95.902 keV 113Cd(n,γ) line
<!---->
- ENSDF 95.902 ± 0.003 keV, E2, α(K) = 1.22, 1.25 per 100 captures; present in every beam run with up to 1.6 million counts.
- Read through each day's production calibration: 01/15 −11 / −4 eV (signal / background); 01/16 −80 eV in the 0.5% runs and −270 eV in the 10% graphite-castle runs, on the same calibration.
- It is not a clean probe of the detector.
  With the Cd shield (captures in the shield, gammas entering from outside) it tracks Pb; without it (captures in the CZT) it reads 31–67 eV low relative to Pb, more at higher rate, and drifts −70 eV through the 01/14 90% run while Pb drifts −20 eV.
- Its Cd-Kα escape at 72.79 keV is too weak to anchor anything (±100–150 eV) and sits under Pb Kα2.
- Every line reads higher deeper in the crystal: +140 to +380 eV at Z = 8.0–8.5 mm against the top 0.25 mm, for Am, Ba and the 96 keV line alike.
  The 96 keV events are not uniform in reconstructed depth; on 01/15 45% sit in the top 0.25 mm, as for Am.
<!---->
### 01/15 and 01/16
<!---->
- Pb is too weak in these beam runs to fit, so there is no clean in-situ line.
- 01/16: the castle runs (27–33k /s) are 100–130 times the rate of the 254 /s calibration run; on the Pb slope that is a shift of a few hundred eV, the size of the 96 keV line's −270 eV.
- 01/15: the 96 keV line reads within about 10 eV of its energy, but only that line checks it.
- A template fit (signal = Ge + smooth factor × background run, background-run noise in the χ²) gives 01/15 68.7533 ± 0.0358 and 01/16 68.7978 ± 0.0532 keV with production's tails; with one tail shared by both Ba lines, 68.763 and 68.807 keV.
<!---->
### Signal and background runs were forced to one Pb centroid
<!---->
- `RooFitUtils::LinkPeakShape` links `Mu` with the shape parameters, contrary to its documentation, since AU bb1721d (2026-05).
  `RunBkgGeSim` uses it to share the Pb shape between the channels, so every in-situ fit has fitted one compromise Pb centroid to both runs; the printed signal and background centroids agree to 0.1 eV in all four pairs.
- With the centroid freed (`PB_ANCHOR=signal`, which links the other eight shape parameters one by one) the signal runs read +4 ± 6 (Cu 01/13), +31 ± 6 (Cd 10%), −4 ± 20 (Cd 25%) and +16 ± 4 eV (Cu 01/14) above their background runs, as independent fits found.
- Ge moves by only −9 to +2 eV, because the background run's Pb, no longer pulled toward the signal run, moves the pair calibration the other way; Method 1 goes from 68.7887 to 68.7891 keV.
- Open: the 01/14 in-situ χ²/ndf rose from 2.34 to 2.62, which a freed parameter should not do.
  Either the fit starts the freed centroid from a default about 0.35 keV away (the saved states come from the tied fit) and lands in a worse minimum, or the signal run's Pb Kα has a component that does not follow the detector (for example fluorescence excited by the activated sample, from elsewhere and at other depths).
<!---->
### The 01/14 outlier
<!---->
- In raw energy, 01/14's Ge sits 99 ± 42 eV higher relative to the Pb of its own run than on 01/13: the outlier, seen without any calibration.
  It is not a gain change: the Pb Kα–Kβ spacing agrees with 01/13 to about 15 eV over 10 keV, at most 10 eV at Ge.
- It is confined to the cathode surface.
  01/14 − 01/13 in (Ge − Pb): +24 ± 70 eV at Z = 8–9 mm, −25 ± 74 at 9–9.5 mm, **+240 ± 73 at 9.5–10 mm**.
  On 01/13, Ge reads 200–250 eV lower at the surface than deeper, as Pb does; on 01/14 it does not (+13 ± 77 eV), while Pb does.
- The Pb depth profile is the same on both days, the top-slice residuals show no extra line with the Ge width fixed to Pb's, and the Ge illumination per crystal is unchanged.
  The excess appears in crystals 0–2 and not in crystal 3 (2σ).
<!---->
## Plan
<!---->
1. **Analysis-Utilities (in progress).**
   `LinkPeakShape` stops tying `Mu`, and old saved states start the freed centroid from the source's; `FitDoublePeak(..., link_sigma, link_tails)` (written, under review); the `PlottingUtils::LegendWidthPx` regression under ROOT 6.41.
   Then `nix flake lock --update-input utils`.
2. **After the pin.**
   Remove the explicit links in `RunBkgGeSim`; rerun `tau20`, `tau20_hard` and `tau20_sigpb`, and check that the 01/14 χ² recovers under the hard fit.
   Decide whether the signal run is calibrated on its own Pb.
3. **Ba-133 anchors.**
   Fit the 79.6/81.0 keV doublet with one shared σ and tail (yield ratio fixed at 0.0805 if needed); fit Am together with Ba 53.16 keV in the 01/16 combined source run; keep 53.16 keV as a validation line.
4. **01/14.**
   Find what moves the surface Ge events (per pixel, Ge − Pb against depth, the detector bias log and the run order); fix the in-situ fit or drop 01/14.
5. **01/15–16.**
   Decide: drop 01/16 for the rate gap, and keep 01/15 with a calibration-transfer systematic of about ±50 eV, unless an in-situ handle turns up.
   The 96 keV line's Te-Kα and Cd-Kβ escapes (68.4–68.7 and 69.3–69.8 keV) sit under and next to Ge and are not yet in the template fit.
6. **Injection sweep.**
   `InjectionTest.cpp` still requires the `gaintransfer` tag, production does not save fit endpoints for it to draw truth from, and trial plots are kept by default.
7. **Write-up.**
   Account for the move from 68.7502 keV step by step: the RooFit fix (≤ 4 eV), fit convergence (Offset, up to 36 eV), the lineshape, the calibration anchors, and the Pb centroid link.
