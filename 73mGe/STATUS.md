# 73mGe — status and plan (2026-11-19)
<!---->
The energy of the bare 73mGe γ ray: in the pixelated CZT at the OSU reactor (`czt/`, January 2026), and in the Amptek CdTe and Si at PSI's BOA beamline (`cdte/`, `si/`, October 2025). ENSDF gives 68.752 ± 0.007 keV. Ge-detector values (68.79–68.82 keV) include a recoil deposit these detectors do not see, so a result there would mean a bias.
<!---->
Why each correction is applied, how large it is, and how the CZT result moved from the 68.7502 ± 0.0145 keV reported earlier: `CORRECTIONS.md`.
<!---->
## Results
<!---->
| Detector | Ge [keV] | ± [keV] | From |
|---|---|---|---|
| CZT, final | **68.7567** | 0.0191 | `czt/CombineGeResult.cpp`; sent to the collaborator 2026-10-27 |
| PSI CdTe | 68.7501 | 0.0331 | `cdte/GeBudget.cpp` → `results/cdte_ge.result` |
| PSI Si | 68.7579 | 0.0174 | `si/SiGeBudget.cpp` → `results/si_ge.result` |
<!---->
The collaborator's draft gives CdTe 68.7484 ± 0.0189, Si-1 68.7543 ± 0.0115 and CZT 68.7606 ± 0.0146 keV.
<!---->
## Running
<!---->
- **Shell.** `nix develop .#73mGe` from the repo root; the lock pins Analysis-Utilities 26.11.19 (a046ecd), which has everything production needs.
- **CZT.** `czt/Filter.cpp` (event selection) and `czt/PixelCalibration.cpp` (pixel gains) first. Then `czt/CalibrateAndFit.cpp` (with `czt/GainTransferFit.hpp` and `czt/PixelAlignment.hpp`) four times: production, and `TRANSFER=offset`, `ESCAPES=off` and `CD_ESCAPES=inside` for its systematics. Then `czt/CombineGeResult.cpp`. For the dissertation's table of what the pixel gain matching does (not a systematic), two more runs, `PIXEL_MAP=off` and `PIXEL_GAINS=off` (tags `tau20_sigpb_esc_pxoff` and `_nogain`, their saved fit states copies of production's), then `czt/GainMatchingTable.cpp`, which writes the table into the dissertation. Each run takes 2–3 minutes and writes `root_files/calibrated_<tag>`, `plots/czt/calibration_low_<tag>` (with the saved fit states its fits start from) and `results/ge_<tag>.result` (columns `method ge_mu fit_err cal_err bkg_err chi2 gain_err cal_off sys_err anchor_err label`). The tags are `tau20_sigpb_esc_pxday`, `tau20_sigpb_esc_oxfer_pxday`, `tau20_sigpb_pxday` and `tau20_sigpb_esc_cdin_pxday`. `czt/FitFigure.cpp` redraws the production fits of the 01/13 signal runs for the paper (`plots/czt/figures`).
- **CdTe.** `cdte/ProcessWaveforms.cpp` (raw waveforms to `features` trees), `cdte/NoGeEscapes.cpp` (the escape tie, entered as constants in `GeBudget.cpp`), then `cdte/GeBudget.cpp`.
- **Si.** `si/ProcessWaveforms.cpp`, then `si/SiGeBudget.cpp`.
- **Recommended value.** `evaluation/RecommendedValue.cpp` prints the averages, draws the comparison of measurements and the Jones and Kraner recalibration (`plots/evaluation`), and writes `results/recommended_value.result` (with the ionization of the 2024 summed peak). Then `evaluation/RecoilEnergy.cpp`: the ⁷³Ge recoil on ENSDF 2019 data against Jones and Kraner's 254.1 eV, the share of the 68.75 keV feeding the two cascades explain, the quenching factor, and the Lindhard comparison up to the kinematic maximum (`results/recoil_energy.result`).
- **Data.** The raw PSI waveforms in `/labdata/PSI/{CdTe,Si}` are read only; the processed trees are in `/labdata/ANSG/PSI/processed/{CdTe,Si}/features`. Line energies and escape tables shared by the three analyses are in `include/ReferenceLines.hpp`.
<!---->
## Reference lines (`include/ReferenceLines.hpp`)
<!---->
One source per quantity, for every analysis:
<!---->
- **Decay γ rays: DDEP.** Am-241 59.5409(1) keV ([table](http://www.lnhb.fr/nuclides/Am-241_tables.pdf)); Ba-133 79.6142(19) and 80.9979(11) keV ([table](http://www.lnhb.fr/nuclides/Ba-133_tables.pdf)).
- **Capture γ ray: CapGam.** ¹¹³Cd(n,γ) 95.902(3) keV ([dataset](https://www.nndc.bnl.gov/capgam/by-target.html?dsid=113CD%28N%2CG%29%20E%3DTHERMAL)).
- **X-ray energies: NIST SRD 128** (Deslattes et al. 2003), the measured value, else the theoretical one. Pb Kα2 72.80542(24) and Kα1 74.97011(17) keV ([Pb](https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Pb&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV), [Cd](https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Cd&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV), [Te](https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Te&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV)).
- **X-ray intensities: Table of Isotopes**, 8th ed. (Firestone 1996), Table 7a ([table](https://grundpraktikum.physik.uni-saarland.de/gpalt/Anleitungen/Ergaenzungen/Versuche/R+oe+ntgenfluoreszenz/X_ray_Energies_Intensities.pdf)). The Kβ groups sit at their intensity-weighted energies; K-O2,3, which SRD 128 does not list, takes Table 7a's energy.
<!---->
## CZT (`czt/`): done
<!---->
Final (2026-10-27): this value was sent to the collaborator and the analysis is closed.
<!---->
| Run | Method | Ge [keV] | fit | cal | sys | transfer model | escapes | escape source | χ²/ndf | BLUE weight | Ge − mean |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cu shield 10%, 01/13 | in-situ | 68.7660 | 20.4 | 7.0 | — | 0 | −12.1 | 0 | 1.33 | 0.510 | +9.3 eV, +0.38σ |
| Cd shield 10%, 01/13 | in-situ | 68.7302 | 31.1 | 7.1 | — | +9.6 | −0.5 | −12.3 | 1.31 | 0.216 | −26.5 eV, −0.74σ |
| Cd shield 25%, 01/13 | in-situ | 68.7435 | 33.2 | 19.8 | — | +11.8 | −1.1 | −13.5 | 1.65 | 0.138 | −13.2 eV, −0.31σ |
| No shield 5%, 01/15 | rate-sub | 68.7771 | 42.4 | 6.2 | 20.4 | 0 | −2.5 | 0 | 1.22 | 0.136 | +20.4 eV, +0.43σ |
<!---->
Errors and corrections in eV; Ge − mean in each run's own total error.
<!---->
- **Quoted** 68.7567 ± 0.0191 keV, +4.7 eV from ENSDF: the BLUE of the four runs, 68.7567 ± 0.0180 keV (χ²/ndf 0.30; statistics 14.4, calibration offset 4.2 and multiplicative 3.6, systematics: rate-sub 2.8, transfer model 3.7, escapes 6.5 and escape source 4.5 eV), with 6.4 eV of outlier removal (the BLUE without Cd 10%, the lowest run) in quadrature. Method 1 (68.7536 ± 0.0193) and Method 2 (68.7771 ± 0.0475 keV) agree to 24 ± 51 eV, the shared escape fractions included.
- **Selection and pixel gains.** Single-interaction events within 2 mm of the cathode (z ≥ 8 mm) on a pixel centre, without pile-up (under 3 µs of live time and in the same millisecond as the previous event), each energy times its pixel's gain: 558.456 keV over its 558 keV centroid in the 01/16 graphite-castle run.
- **Day pixel map.** Near 70–96 keV the 558 keV gains leave the crystals 80–300 eV apart and the pixels 110 eV rms (`CORRECTIONS.md`). Every event is aligned on its day's map (`czt/PixelAlignment.hpp`): each pixel's line position against its run's summed spectrum, averaged over the day's beam runs, on Pb Kα on 01/13 and on the 95.902 keV line on 01/15, where Pb Kα barely stands above the continuum. A pixel measured to 100 eV or better at 68.75 keV gets the factor μ_run/μ_pixel; a day's source runs take its map.
- **Method 1, 01/13, in-situ.** Each background run (Pb Kα2/Kα1) is fitted together with its signal run (the same lines plus Ge) over 65–80 keV in an unbinned likelihood, first in raw energy, with one Pb shape and one low tail (a Gaussian with an exponential tail, τ/σ ≤ 20) on a flat continuum. Cu 10% is calibrated with a pol1 through Am-241 and its background run's Pb Kα; the Cd pairs inherit it as a gain through their own Pb Kα, and each signal run is shifted onto its own Pb Kα.
- **Method 2, 01/15, rate subtraction.** A pol1 through Am-241 and Ba-133 80.998 keV from the source runs; Ge from the signal − background residual of rates (each crystal's counts over its own live time, the live intervals of all its events) over 62–71.5 keV, in comb-period bins at five starting points, the fit error scaled by √(χ²/ndf). Systematic: the spread over starting points (8.2 eV) and the shift to a quadratic residual background (18.7 eV).
- **Both methods** carry the 95.902 keV line's Cd and Te K X-ray escapes as fixed lines. Whether the Cd-shield rows' escapes come from the shield or from captures inside the crystal is open, so their difference is a systematic.
- **Combination** (`CombineGeResult`): a BLUE with correlated blocks for the Am reference (the 01/13 rows), the gain transfer, the transfer model (offset − gain transfer), the escape correction (with − without the escapes, the Cd-shield runs and the rest as two blocks) and the Cd-shield rows' escape source (inside − shield); fits, anchors and rate-sub systematics on the diagonal; the outlier removal on top.
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
The Si (DT5780 channel 0) on the trapezoid of the processed waveforms (rise 4 µs, flat top 1 µs, baseline up to 75 samples before the trigger; `si/SiShapingScan.cpp` chose it by the Am and Pb Kα widths), built like the collaborator's Si-1: the Ge runs of 10/09–10/12 in the same six day/night groups, each fitted with one Gaussian on a linear background over 66.5–71 keV and calibrated with Am-241 from the source runs and the group's own Pb Kα, then the groups' error-weighted mean. Errors: fit, calibration (the Am position shared by the groups, the Pb positions each group's own), background (|linear − flat|), window (RMS over 15 windows of the summed spectrum).
<!---->
- **Result** 68.7579 ± 0.0174 keV (fit 14.8, calibration 2.2, background 8.6, window 2.5 eV). The six groups give χ² 3.2/5.
- **Against the draft's groups** (ours − draft): Thu 10/09 −39.4, Fri night +62.8, Sat day −9.9, Sat night +8.5, Sun day +27.6, Sun night +2.3 eV, each group with a 25–74 eV error. The means differ by +3.6 eV (Si-1 68.7543 ± 0.0115). The draft does not say how each group was fitted. Plot: `plots/si/ge_budget/psi_si_ge_groups.png`.
- **The summed spectrum** of the 10/10–10/12 runs reads 68.7653 ± 0.0205 keV, +7.4 eV from the groups.
- **The digitiser's on-board energy is not used.** On the summed spectrum the two energies differ by 4.9 ± 15.1 eV (bootstrap on the same records), so the choice is not a systematic. The trapezoid has the steadier Pb Kα1 run to run (χ² 6.4/7 against 18.1/7), per-run Ge χ² 9.8/6 against 13.9/6, physical per-run widths and no channel comb.
- **Checks** on the summed spectrum (eV from its Am + Pb Kα value): Am + Ba-133 calibration +4.5; per-day mean −5.9; quadratic against linear Am + Ba −6.5; isolated records +2.3; the 53.44 keV line 53.4201 ± 0.0319 keV.
- **Every fit is checked for convergence** (MIGRAD and HESSE status 0, full covariance). Only the Am fit of one low-count source run fails, and the line-stability check skips that run.
- **Records whose baseline sits at the bottom of the ADC range** are dropped (`CORRECTIONS.md`).
<!---->
## Findings behind open questions
<!---->
- **Calibration curvature.** 01/13: the in-run Pb Kα1 → Kβ1 slope over 75–85 keV (0.98817 ± 0.00052) matches Am → Kα (0.98822), so a quadratic of that size moves Ge by about 1 eV. 01/15: Ba-133 53.16 keV reads −5.6 ± 16.7 eV on the linear calibration; above 81 keV the response is strongly nonlinear.
- **Run-to-run shifts.** The scale moves between runs by tens of eV, up to 130 eV at the 01/13 Cd → Cu shield change. No count-rate model holds across days.
- **Depth.** Lines read higher deeper in the crystal (Am 220 eV, Pb Kα1 130 eV across four depth slices); a per-depth calibration left Method 1 near 68.763 keV in throwaway fits, before the pixel map.
- **The energy comb** does not affect the unbinned in-situ fits: dithering every energy by ±132 eV moves the 01/13 rows by −8 to +4 eV and Method 1 by about −3 eV.
- **Testbed** (`czt/InjectionTest.cpp`, `czt/SummarizeInjection.cpp`). A line injected into a synthetic replica of each pair is recovered: Cu +1.8 ± 2.7 eV (50 trials), Cd 10% −11.4 ± 8.3, Cd 25% −4.3 ± 9.2 eV. On halves of the real background runs, rate-sub closes (−2.1 ± 2.7 eV together), while in-situ reads +12 to +33 eV. Halves of eight synthetic background runs, each held fixed through its trials, scatter by 16.0 (Cu), 21.1 (Cd 10%) and 9.1 eV (Cd 25%) from one background run to the next, so the real runs' offsets are that run's fluctuation, not a real-background effect. These numbers are from 10/23; a 3-trial recheck on the current code (2026-10-27, `plots/testbed/injection_recheck_*`) agrees for the synthetic replicas and the rate-sub halves.
<!---->
## Plan
<!---->
1. **The CdTe extrapolation term** (dissertation only: the letter cites the collaborator's preprint, arXiv:2609.20848, for the PSI values, and our CdTe and Si reprocessing goes only in the dissertation): the full curvature-free shift or the curvature's quadrature cost, and whether to move the slice edges off the rise estimator's sawtooth peaks (280, 300, 320 ns).
2. **Recommended value** for the paper: an error-weighted average of the CZT (this work), the collaborator's CdTe and Si, and Jones and Kraner (1971, 1975), as in the collaborator's paper (adopted 68.7536 ± 0.0034 keV there); our CdTe and Si (the collaborator's six Si groups as they are) are cross-checks of the same data. Only measurements without a Ge detector enter: neutrons cannot be kept out of a detector (the CZT sees ¹¹³Cd(n,γ) from its own crystal in every configuration, behind the Cd foil too), and captures on a Ge detector's own ⁷²Ge add the recoil to the 68.75 keV line. That leaves out Collar et al. (with it the value moves −0.6 eV) and any Jones and Kraner value taken with a Ge(Li) detector. Their 1971 paper (Phys. Rev. C 4, 125) lists seven energies, (p,p′γ) ×4, (n,γ), (α,α′γ) ×2, from "both Ge(Li) and Si(Li) detectors" without saying which; only a (p,p′) spectrum is shown as Si(Li). The 1975 value (Phys. Rev. A 11, 1347), 68.755 ± 0.005 keV, is ⁷²Ge(n,γ) on a GeO₂ sample seen with a Si(Li). Both papers calibrate on Am 59.537 keV and Bearden's Pb Kα; on DDEP and SRD 128 their energies rise by +1.99 to +2.17 eV for any linear or quadratic calibration through those three lines (their lines on a fake ADC, linear on their scale, recalibrated every way), so they take +2.08 ± 0.09 eV (`evaluation/RecommendedValue.cpp`, `plots/evaluation`). With the 1971 (p,p′γ) group (its spectrum shown taken with a Si(Li)) and the 1975 value: 68.7547 ± 0.0037 keV; the PSI values, taken to be calibrated like the collaborator's 2024 paper (Phys. Rev. D 110, 083043: DDEP Am, Bearden Pb Kα; the 2026 preprint does not say), take +0.64 ± 0.21 eV (lines with Am weighted at least as each Pb line: 0.42–0.85 eV), giving **68.7548 ± 0.0037 keV** (χ²/ndf 0.8/5; 68.7532 unadjusted). The 2024 summed peak, 68.8158 keV on the same energies, takes the same +0.64 eV: ionization 61.6 ± 5.9 eV, quenching factor 24.2 ± 2.3 % (254.1 eV recoil). The Si inputs are the draft's results table (Si-1 68.7529 ± 0.0114, Si-2 68.7555 ± 0.0128 keV; its per-group values average to 68.7543); nothing more is asked of the collaborator (2026-11-19, deadlines: dissertation 12/01, papers 12/15). Open: correlations among the Jones and Kraner values (J&K 1975 carries 56% of the weight), treated as independent for now.
3. **Write-up**: a PRC letter (the measurement, the CZT fit, the non-Ge recommended value), drafted 2026-11-19 in `~/Papers/Writing/WIP-73mGe-PRC` (`experiment.tex`, `gamma.tex`, figures from `czt/FitFigure.cpp` and `evaluation/RecommendedValue.cpp`); `CORRECTIONS.md` stays out of it.
