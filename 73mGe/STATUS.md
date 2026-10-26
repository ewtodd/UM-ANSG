# 73mGe — status and plan (2026-10-26)
<!---->
The energy of the bare 73mGe γ ray: in the pixelated CZT at the OSU reactor (`czt/`, January 2026), and in the Amptek CdTe and Si at PSI's BOA beamline (`cdte/`, `si/`, October 2025). ENSDF gives 68.752 ± 0.007 keV. Ge-detector values (68.79–68.82 keV) include a recoil deposit these detectors do not see, so a result there would mean a bias.
<!---->
Why each correction is applied, how large it is, and how the CZT result moved from the 68.7502 ± 0.0145 keV reported earlier: `CORRECTIONS.md`.
<!---->
## Results
<!---->
| Detector | Ge [keV] | ± [keV] | From |
|---|---|---|---|
| CZT, quoted | **68.7511** | 0.0211 | `czt/CombineGeResult.cpp` on `tau20_sigpb_esc` |
| PSI CdTe | 68.7501 | 0.0331 | `cdte/GeBudget.cpp` → `results/cdte_ge.result` |
| PSI Si | 68.7648 | 0.0205 | `si/SiGeBudget.cpp` → `results/si_ge.result` |
| CZT + CdTe | — | — | `COMBINE_WITH_CDTE`, off until the CdTe is final |
<!---->
The collaborator's draft gives CdTe 68.7484 ± 0.0189, Si-1 68.7543 ± 0.0115 and CZT 68.7606 ± 0.0146 keV.
<!---->
## Running
<!---->
- **Shell.** `nix develop .#73mGe` from the repo root; the lock pins Analysis-Utilities 26.10.26 (3465b43), which has everything production needs.
- **CZT.** `czt/CalibrateAndFit.cpp` (with `czt/GainTransferFit.hpp`) per configuration, then `czt/CombineGeResult.cpp`. Production is the default configuration, `tau20_sigpb_esc`; its partners for the systematics are `TRANSFER=offset` (`tau20_sigpb_esc_oxfer`) and `ESCAPES=off` (`tau20_sigpb`). Each writes `root_files/calibrated_<tag>`, `plots/czt/calibration_low_<tag>` and `results/ge_<tag>.result` (columns `method ge_mu fit_err cal_err bkg_err chi2 gain_err cal_off sys_err anchor_err label`). A new tag needs `plots/czt/calibration_low_<tag>/fits` copied from production's, or the fitter opens the editor.
- **CdTe.** `cdte/ProcessWaveforms.cpp` (raw waveforms to `features` trees), `cdte/NoGeEscapes.cpp` (the escape tie, entered as constants in `GeBudget.cpp`), then `cdte/GeBudget.cpp`.
- **Si.** `si/ProcessWaveforms.cpp`, then `si/SiGeBudget.cpp`.
- **Data.** The raw PSI waveforms in `/labdata/PSI/{CdTe,Si}` are read only; the processed trees are in `/labdata/ANSG/PSI/processed/{CdTe,Si}/features`. Line energies and escape tables shared by the three analyses are in `include/ReferenceLines.hpp`.
<!---->
## CZT (`czt/`)
<!---->
| Run | Method | Ge [keV] | fit | cal | sys | transfer model | escapes | χ²/ndf | BLUE weight | Ge − mean |
|---|---|---|---|---|---|---|---|---|---|---|
| Cu shield 10%, 01/13 | in-situ | 68.7621 | 20.6 | 7.2 | — | 0 | −12.3 | 1.39 | 0.478 | +11.0 eV, +0.44σ |
| Cd shield 10%, 01/13 | in-situ | 68.7193 | 30.8 | 7.3 | — | +10.0 | −0.7 | 1.14 | 0.261 | −31.8 eV, −0.96σ |
| Cd shield 25%, 01/13 | in-situ | 68.7611 | 33.4 | 20.1 | — | +11.9 | −0.8 | 1.46 | 0.165 | +9.9 eV, +0.24σ |
| No shield 5%, 01/15 | rate-sub | 68.7659 | 42.7 | 6.3 | 27.4 | 0 | −10.0 | 1.43 | 0.096 | +14.8 eV, +0.28σ |
<!---->
Errors and corrections in eV; Ge − mean in each run's own total error.
<!---->
- **Quoted** 68.7511 ± 0.0211 keV, −0.9 eV from ENSDF: the BLUE of the four runs, 68.7511 ± 0.0178 keV (χ²/ndf 0.43; statistics 14.4, calibration offset 4.6, multiplicative 3.8, rate-sub systematics 2.6, transfer model 4.6, escapes 6.8 eV), with 11.2 eV of outlier removal (the BLUE without Cd 10%, the lowest run) in quadrature. Method 1 (68.7498 ± 0.0186) and Method 2 (68.7659 ± 0.0521 keV) agree to 16 ± 55 eV.
- **Method 1, 01/13, in-situ.** Each background run (Pb Kα2/Kα1) is fitted together with its signal run (the same lines plus Ge) in raw energy, with one Pb shape and one low tail (a Gaussian with an exponential tail, τ/σ ≤ 20). Cu 10% is calibrated with a pol1 through Am-241 and its background run's Pb Kα; the Cd pairs inherit it as a gain through their own Pb Kα, and each signal run is shifted onto its own Pb Kα.
- **Method 2, 01/15, rate subtraction.** A pol1 through Am-241 and Ba-133 80.998 keV from the source runs; Ge from the live-time-normalised signal − background residual over 62–71.5 keV, in comb-period bins at five starting points, the fit error scaled by √(χ²/ndf). Systematic: the spread over starting points, a quadratic residual background, and 20 eV for carrying the source-run calibration to the beam runs.
- **Both methods** carry the 95.902 keV line's Cd and Te K X-ray escapes as fixed lines.
- **Combination** (`CombineGeResult`'s defaults: `tag = "tau20_sigpb_esc"`, `alt_tag = "tau20_sigpb_esc_oxfer"`, `escape_tag = "tau20_sigpb"`): a BLUE with correlated blocks for the Am reference (the 01/13 rows), the gain transfer, the transfer model (`alt_tag` − `tag`) and the escape correction (`tag` − `escape_tag`, the Cd-shield runs and the rest as two blocks); fits, anchors and rate-sub systematics on the diagonal; the outlier removal on top.
- **01/14 and 01/16 are not used** (`CORRECTIONS.md`): 01/14's calibration cannot be carried over from 01/13, and the 01/16 graphite castle differs from the other runs in too many ways.
<!---->
## PSI CdTe (`cdte/`)
<!---->
The CdTe (DT5780 channel 1) on `PulseHeightProcessor` energies, in rise-time slices because slow pulses lose charge to hole trapping. Each slice is fitted over 62–89 keV with every line at its true energy on a free scale (Pb Kα2/Kα1, the Pb Kβ groups, the 95.902 keV line's Cd and Te escapes, Ge free), one EMG with a linear width and a linear continuum. The result is the weighted mean of the five slices up to 320 ns, on the five Ge beam runs (10/09–10/11).
<!---->
| Term | eV | |
|---|---|---|
| statistics | 16.4 | weighted mean of five slices, 6.6 eV, × √(χ²/ndf) 2.47 of the slices ≤ 300 ns (18.4/3) |
| outlier removal | 10.7 | largest shift without the highest (0–270 ns), the lowest (290–300 ns), or both slices |
| lineshape | 5.5 | Ge's tail fitted on its own |
| escapes | 4.5 | Te Kα/Cd Kβ ± its error and per slice; Cd Kβ/Kα 0.324 against Geant4's 0.275 |
| Pb Kβ blend | 1.4 | Kβ3/Kβ1 (and K-N2/K-N3) ± 10%; Kβ5, K-N4,5, K-O2,3 ± 20% |
| extrapolation | 22.8 | the scale's curvature free (Ge lies below the lines that fix the scale) |
| slice range | 0.8 | slices ≤ 300 ns instead (statistical spread 1.1 eV) |
| slice combination | 4.7 | a random-effects mean of the slices (between-slice spread 32.7 eV) |
<!---->
- **Extrapolation** is the full shift with the curvature free. It is a noisy statistic (20.6 ± 9.5 eV over replicas of the same data), mostly the 0–270 and 300–320 ns slices, where the curvature is nearly degenerate with the escapes; the curvature's quadrature cost instead gives ± 0.0249 keV (Plan 5).
- **10/10 carries most of the slices' disagreement**: without its two Ge-target files (a scratch rerun) the slices ≤ 300 ns agree at χ² 9.0/3 instead of 18.4/3 and Ge reads 68.7537 ± 0.0283 keV. 10/10 stays in.
- **The move to `PulseHeightProcessor`** (+5.3 eV, from 68.7448) is reprocessing scatter, not a correction (`cdte/ReprocessingComparison.cpp`, `CORRECTIONS.md`). The old `trap/` trees are kept until that comparison is no longer needed.
<!---->
## PSI Si (`si/`)
<!---->
The Si (DT5780 channel 0) on the trapezoid of the processed waveforms (rise 4 µs, flat top 1 µs, baseline up to 75 samples before the trigger; `si/SiShapingScan.cpp` chose it by the Am and Pb Kα widths), the collaborator's way: one Gaussian on a linear background over 66.5–71 keV in the summed main Ge runs (10/10–10/12; BeamTurnOn_4 of 10/09 is outside the result), calibrated with Am-241 from the source runs and the beam runs' own Pb Kα (cal A). Errors: fit, calibration offset and slope, background (|linear − flat|), window (RMS over 15 windows).
<!---->
- **Result** 68.7648 ± 0.0205 keV (fit 17.2, calibration 2.2, background 10.6, window 2.5 eV).
- **The digitiser's on-board energy is not used.** On the same records it read 68.7698 ± 0.0196 keV, −4.9 ± 15.1 eV from the trapezoid by bootstrap, so the choice is not a systematic. The trapezoid has the steadier Pb Kα1 run to run (χ² 6.4/7 against 18.1/7), per-run Ge χ² 9.8/6 against 13.9/6, physical per-run widths and no channel comb.
- **Checks** (eV from cal A): cal B (Am + Ba-133) +5.0; per-day mean −5.9; without GeTarget_1 −5.5; with the W Kβ group −10.4; quadratic cal B −6.5 against linear cal B; isolated records +2.3; the 53.44 keV line 53.4201 ± 0.0319 keV.
- **GeTarget_1** (the first 1.6 h of the 10/10 acquisition) reads 2.9σ high, but the seven runs give χ² 9.8/6 and nothing besides its value argues for dropping it.
- **The collaborator's Si-1** averages six day/night runs from 10/09 to 10/12, so it includes 10/09.
- **Records whose baseline sits at the bottom of the ADC range** are dropped (`CORRECTIONS.md`).
<!---->
## Findings behind open questions
<!---->
- **Calibration curvature.** 01/13: the in-run Pb Kα1 → Kβ1 slope over 75–85 keV (0.98817 ± 0.00052) matches Am → Kα (0.98822), so a quadratic of that size moves Ge by about 1 eV. 01/15: Ba-133 53.16 keV reads −5.6 ± 16.7 eV on the linear calibration; above 81 keV the response is strongly nonlinear.
- **Run-to-run shifts.** The scale moves between runs by tens of eV, up to 130 eV at the 01/13 Cd → Cu shield change. No count-rate model holds across days.
- **Depth and edges.** Lines read higher deeper in the crystal (Am 220 eV, Pb Kα1 130 eV across four depth slices); in each crystal's outer ring and along the gap Pb Kα1 reads 7–155 eV off. A per-depth calibration leaves Method 1 near 68.763 keV in throwaway fits.
- **The energy comb** does not affect the unbinned in-situ fits: dithering every energy by ±132 eV moves the 01/13 rows by −8 to +4 eV and Method 1 by about −3 eV.
- **Testbed** (`czt/InjectionTest.cpp`, `czt/SummarizeInjection.cpp`). A line injected into a synthetic replica of each pair is recovered: Cu +1.8 ± 2.7 eV (50 trials), Cd 10% −11.4 ± 8.3, Cd 25% −4.3 ± 9.2 eV. On halves of the real background runs, rate-sub closes (−2.1 ± 2.7 eV together), while in-situ reads +12 to +33 eV; a shared fluctuation of the one background run explains that (χ² 3.9/3) but does not exclude a real-background effect (Plan 4).
<!---->
## Plan
<!---->
1. **Method 1's remaining model choices** are not in the budget: the in-situ background shape (flat against sloped, −1 to +6 eV on Method 1), the lineshape (shared tail, τ/σ ≤ 20) and the pixel selection. Decide which become systematics, as the transfer model did.
2. **The Cd pairs** inherit the Cu reference across the 130 eV step at the Cd → Cu change; if the bias was cycled then, they have 01/14's problem.
3. **The Cd-shield rows' escapes.** Both capture lines say the Cd runs see captures in the CZT's cathode layer as the other runs do (`CORRECTIONS.md`): the 95.9 keV line sits there in every configuration, and against the 558 keV line it is as inside-type in the Cd runs (0.88–0.92 of Cu's) as in the no-shield run (0.85–0.92), which carries the inside escapes. So the Cd rows likely need the inside escapes too (roughly −4 eV on the quoted value, an estimate), against the missing Cd Kβ escape in their background runs (0 ± 0.004 where inside escapes give 0.012).
4. **The in-situ fit on real background**: running the halves procedure on several synthetic background runs, each held fixed through its trials, would measure the spread the shared fluctuation causes; a real-background effect would show beyond it.
5. **The CdTe extrapolation term**: the full curvature-free shift or the curvature's quadrature cost, and whether to move the slice edges off the rise estimator's sawtooth peaks (280, 300, 320 ns).
6. **The Si in the combination**: whether it joins, and with which correlations (the Si and CdTe share the PSI beam time and the Pb Kα and Am reference energies, not the calibration fits).
7. **Write-up**: the step-by-step move from 68.7502 keV (`CORRECTIONS.md`), with each step's size measured.
