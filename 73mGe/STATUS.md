# 73mGe — status and plan (2026-10-27)
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
| PSI Si | 68.7648 | 0.0205 | `si/SiGeBudget.cpp` → `results/si_ge.result` |
<!---->
The collaborator's draft gives CdTe 68.7484 ± 0.0189, Si-1 68.7543 ± 0.0115 and CZT 68.7606 ± 0.0146 keV.
<!---->
## Running
<!---->
- **Shell.** `nix develop .#73mGe` from the repo root; the lock pins Analysis-Utilities 26.10.26 (3465b43), which has everything production needs.
- **CZT.** `czt/Filter.cpp` (event selection) and `czt/PixelCalibration.cpp` (pixel gains) first. Then `czt/CalibrateAndFit.cpp` (with `czt/GainTransferFit.hpp` and `czt/PixelAlignment.hpp`) four times: production, and `TRANSFER=offset`, `ESCAPES=off` and `CD_ESCAPES=inside` for its systematics. Then `czt/CombineGeResult.cpp`. Each run takes 2–3 minutes and writes `root_files/calibrated_<tag>`, `plots/czt/calibration_low_<tag>` (with the saved fit states its fits start from) and `results/ge_<tag>.result` (columns `method ge_mu fit_err cal_err bkg_err chi2 gain_err cal_off sys_err anchor_err label`). The tags are `tau20_sigpb_esc_pxday`, `tau20_sigpb_esc_oxfer_pxday`, `tau20_sigpb_pxday` and `tau20_sigpb_esc_cdin_pxday`.
- **CdTe.** `cdte/ProcessWaveforms.cpp` (raw waveforms to `features` trees), `cdte/NoGeEscapes.cpp` (the escape tie, entered as constants in `GeBudget.cpp`), then `cdte/GeBudget.cpp`.
- **Si.** `si/ProcessWaveforms.cpp`, then `si/SiGeBudget.cpp`.
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
| Cu shield 10%, 01/13 | in-situ | 68.7660 | 20.4 | 7.0 | — | 0 | −12.2 | 0 | 1.05 | 0.509 | +9.3 eV, +0.38σ |
| Cd shield 10%, 01/13 | in-situ | 68.7302 | 31.1 | 7.2 | — | +9.7 | −0.6 | −12.3 | 1.48 | 0.216 | −26.5 eV, −0.74σ |
| Cd shield 25%, 01/13 | in-situ | 68.7435 | 33.2 | 19.8 | — | +11.8 | −1.1 | −13.5 | 1.54 | 0.138 | −13.1 eV, −0.31σ |
| No shield 5%, 01/15 | rate-sub | 68.7771 | 42.4 | 6.2 | 20.4 | 0 | −2.6 | 0 | 1.22 | 0.136 | +20.4 eV, +0.43σ |
<!---->
Errors and corrections in eV; Ge − mean in each run's own total error.
<!---->
- **Quoted** 68.7567 ± 0.0191 keV, +4.7 eV from ENSDF: the BLUE of the four runs, 68.7567 ± 0.0180 keV (χ²/ndf 0.30; statistics 14.4, calibration offset 4.2 and multiplicative 3.6, systematics: rate-sub 2.8, transfer model 3.7, escapes 6.6 and escape source 4.5 eV), with 6.4 eV of outlier removal (the BLUE without Cd 10%, the lowest run) in quadrature. Method 1 (68.7535 ± 0.0193) and Method 2 (68.7771 ± 0.0475 keV) agree to 24 ± 51 eV, the shared escape fractions included.
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
- **Depth.** Lines read higher deeper in the crystal (Am 220 eV, Pb Kα1 130 eV across four depth slices); a per-depth calibration left Method 1 near 68.763 keV in throwaway fits, before the pixel map.
- **The energy comb** does not affect the unbinned in-situ fits: dithering every energy by ±132 eV moves the 01/13 rows by −8 to +4 eV and Method 1 by about −3 eV.
- **Testbed** (`czt/InjectionTest.cpp`, `czt/SummarizeInjection.cpp`). A line injected into a synthetic replica of each pair is recovered: Cu +1.8 ± 2.7 eV (50 trials), Cd 10% −11.4 ± 8.3, Cd 25% −4.3 ± 9.2 eV. On halves of the real background runs, rate-sub closes (−2.1 ± 2.7 eV together), while in-situ reads +12 to +33 eV. Halves of eight synthetic background runs, each held fixed through its trials, scatter by 16.0 (Cu), 21.1 (Cd 10%) and 9.1 eV (Cd 25%) from one background run to the next, so the real runs' offsets are that run's fluctuation, not a real-background effect. These numbers are from 10/23; a 3-trial recheck on the current code (2026-10-27, `plots/testbed/injection_recheck_*`) agrees for the synthetic replicas and the rate-sub halves.
<!---->
## Plan
<!---->
1. **The CdTe extrapolation term**: the full curvature-free shift or the curvature's quadrature cost, and whether to move the slice edges off the rise estimator's sawtooth peaks (280, 300, 320 ns).
2. **The Si in the combination**: whether it joins, and with which correlations (the Si and CdTe share the PSI beam time and the Pb Kα and Am reference energies, not the calibration fits).
3. **Write-up**: the step-by-step move from 68.7502 keV (`CORRECTIONS.md`), with each step's size measured.
