# 73mGe — the corrections and why (2026-10-24)
<!---->
What each correction in the result does, how large it is, and the evidence for it; the current numbers and the plan are in `STATUS.md`.
A shift between two analysis choices is charged as a systematic only when more than one of the choices can be defended; fixing a defect, or changing an implementation detail that should give the same answer, is not a systematic (the CdTe reprocessing, below).
<!---->
## CZT: from 68.7502 to 68.7533 keV
<!---->
The value reported earlier, 68.7502 ± 0.0145 keV (c260617), moved through these steps (sizes where measured; the full table is for the write-up, Plan 7 in `STATUS.md`):
<!---->
1. **The Pb Kβ1,3 blend in the old pol2** was given Kβ1's energy, 84.936 keV, where the blend's centroid is about 84.770 (Kβ3 at 0.52 of Kβ1). That pulled Cu 01/13's Ge 27.4 eV below the pol1; with the centroid it sits 5.4 eV above. Production uses a pol1 through Am-241 and Pb Kα only.
2. **The RooFit fix**: ≤ 4 eV.
3. **Fit convergence**: up to 36 eV.
4. **The lineshape**: a Gaussian with a low exponential tail, τ/σ ≤ 20. The earlier τ/σ ≤ 8 railed at its bound (2026-10-22 tests), so that lineshape was set by the bound, not fitted.
5. **The Pb centroid link** between the background and signal runs.
6. **01/14 dropped**: the detector was probably unbiased overnight and re-biased that morning, and there was no Am run, so no transfer from 01/13 applies.
7. **The signal-run Pb anchor** (below).
8. **The gain transfer model** (below).
9. **Both methods in the combination**: the gap between them was not significant (0.3–0.5σ), and charging it as a systematic used Method 2 only through that noise.
10. **The 95.902 keV escapes** (below): −7.7 eV on the quoted value.
<!---->
How the anchor (`_sigpb`), the transfer model (`_oxfer`) and the escapes (`_esc`) move each run: Ge − 68.752 keV per in-situ row in eV; the rate-sub rows are +23.9 (01/15) and +42.6 eV (01/16) without the escapes, +14.2 and +33.3 eV with them. Each tag combined alone.
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
The quoted 68.7533 keV is `tau20_sigpb_esc` combined with its transfer-model and escape partners as correlated terms (`STATUS.md`).
<!---->
## CZT corrections
<!---->
### Signal-run Pb anchor
<!---->
Each 01/13 signal run is shifted onto its own Pb Kα, so it inherits only its pair's slope.
<!---->
- Ge follows its own run's Pb: in template fits Ge − Pb within a signal run is the same in the Cd and Cu runs to ±45 eV while Pb moves up to 213 eV.
- The Cd 10% signal run reads +36 ± 11 eV above its background run within the same crystal, pixel ring and depth slice (where the Pb lands contributes +0.2 eV): a shift of the whole run. Cu 01/13 and Cd 25% show none (+5 ± 12, +11 ± 25 eV). The cause is unknown.
- Pb Kβ, 10 keV above the window and not used by the anchor, follows: on the anchored calibration the Cd 10% signal run's Kβ reads +46 ± 30 eV against its background run's, where a shift of the Kα lines alone would leave it at −31 eV (2.6σ away) and a shift of the whole run at 0 (1.5σ). Cu 01/13 and Cd 25% read −28 ± 37 and +30 ± 77 eV, Kα itself 0 to +7 eV in all three (comb-period bins at 20 starting points; `/tmp/73mGe-escape/PbKbetaCheck.C`, plots in `plots/czt/studies/pb_kbeta_check`).
- So Cd 10%'s low Ge, −34 eV from the mean and −1.0σ, is not an artefact of the anchor.
<!---->
### Gain transfer to the Cd pairs
<!---->
The Cd pairs inherit the Cu 10% calibration through their background runs' Pb Kα, as a gain (production) or as an offset (`TRANSFER=offset`). Pb Kβ lies 10 keV above the calibrated range and the Cd background runs' in-run Kα1 → Kβ1 slopes split between the two; no line inside the range resolves the 0.15% slope difference. Both can be defended, so the difference (+9.9 and +11.7 eV on the Cd rows) is a correlated systematic in the combination.
<!---->
### 95.902 keV escapes
<!---->
- **The line** is 113Cd(n,γ)'s E2 transition in 114Cd, converted in the K shell 47% of the time and emitted as a photon 37% (BrIcc: αK 1.288, αT 1.735). Its escapes come from the Cd K X-rays of conversions at captures near a surface, and from Cd and Te K X-rays when its photon is absorbed: Cd Kα at 72.73 / 72.92 keV under Pb Kα2, Cd Kβ at 69.25–69.81 keV beside Ge, Te Kα at 68.43 / 68.70 keV under Ge, Te Kβ at 64.2–64.9 keV.
- **Captures inside the crystal** (no shield, Cu shield), as fractions of the line's Gaussian yield: Cd Kβ 0.0119, the mean of the 01/15 and 01/16 beam background runs (0.0122, 0.0115, each ± ~0.0017, the one escape clear of other lines); Cd Kα 0.0450 from it through the DDEP intensities times escape probabilities; Te Kα 0.0024 from Geant4's 0.053 of Cd Kα (the runs give 0.0006 and 0.0021, ± 0.002–0.004); Te Kβ 0.0007. Geant4 alone (CZT, captures at 0.25 mm mean depth) gives 0.061, 0.017, 0.0033 and 0.0010.
- **Photons from the Cd shield**, which stops the thermal neutrons: Geant4 with the line's photons onto the face, 0.0115, 0.0033, 0.0057 and 0.0017. The 01/13 Cd background runs show no Cd Kβ escape (0 ± 0.004 at 10%, ± 0.006 at 25%), where captures inside would give 0.012; the Cu run cannot tell (± 0.008–0.026).
- **In the fits** the escapes are fixed lines (`RooFitUtils::AddFixedLines`) at the run's own 95.902 keV centroid minus the X-ray energy (in the 01/16 background run the line and its escapes read low together, −244 eV), with that line's Gaussian yield times the fractions. In the in-situ fits each channel gets its own run's lines, drawn with the Pb Kα1 shape; the line's yield and centroid come from its fits in comb-period bins at five starting points over 91.5–100.5 keV, with the in-situ fit's converged tail (the raw fit first runs once without the escapes). In the rate-sub fit they are a fixed term with Ge's width and tail, scaled by the residual's own 95.902 keV line, which is 1.2 times stronger in the signal runs than live time predicts.
- **Effect**: Cu 01/13 −11.7 eV, the Cd pairs −0.7 and −0.8 eV, the rate-sub rows −9.7 and −9.2 eV; the quoted value −7.7 eV (68.7630 ± 0.0191 keV without, 68.7553 ± 0.0183 with them as a correction only). A TF1 replica of the in-situ fit with the same lines gives −12.7, +0.1 and −1.2 eV.
- **Charged in full** as a systematic, correlated among the runs whose escapes come from captures inside the crystal (Cu, 01/15, 01/16) and among the Cd-shield runs, independent between the two sets.
- **Sensitivity**: with Te Kα (and Kβ) doubled, as PSI's CdTe runs without Ge suggest, the rows move +2.4 (Cu) and +1.4 / +3.1 eV (rate-sub), well inside that systematic.
- **Validation.** The fixed lines reproduce a testbed's doublet without bias (−0.1 ± 2.1 and −0.3 ± 1.4 eV; +12.4 and +7.6 eV without them), and with the escapes switched off production reproduces `tau20_sigpb` exactly. In production's own fit, on a synthetic replica of Cu 01/13 with the escapes in the data, the fit with the escape lines recovers the injected 68.750 keV (+1.8 ± 2.7 eV, 50 trials); the same events fitted without them read +11.1 ± 0.1 eV higher in every trial, the size of production's Cu correction.
<!---->
### Rate subtraction (Method 2)
<!---->
- **Without the subtraction** the fits read +91 (01/15) and +134 eV (01/16), and a template with a free scale degenerates (it fits 0.3–1.3 against live-time ratios of 3.6 and 1.8): the subtraction removes real structure near Ge.
- **The energy comb.** Every pixel's energies sit on their own 1–8 eV grid; summed, they show a comb with a 264 eV period. Bins near that period beat against it: 255 eV bins move rate-sub Ge by 13–16 eV rms with where they start, comb-period bins by 8–10 eV, hence comb-period bins at five starting points, averaged, with their spread in the systematic.
- **Overdispersion.** The signal and background runs fill the grids differently, so signal − α × background is 7–16 times Poisson at 0.5–10 eV bins, about 4 at 50 eV, 1.4–1.9 at 255 eV; unbinned or fine-binned two-run fits fail, and from ~20 eV bins up an on/off likelihood agrees with the subtraction. The fit error is scaled by √(χ²/ndf).
- **Source runs to beam runs.** 01/15–16 are calibrated on source runs; no line appears in both, so the step cannot be checked there. On 01/13 the same step holds to about 20 eV (the in-run Kβ slope), which is the systematic charged.
<!---->
## PSI CdTe
<!---->
- **Per-slice calibration.** Slow pulses lose charge to hole trapping, so Ge is fitted in 10–90% rise-time slices, each calibrated on its own lines; the result uses the slices up to 320 ns.
- **Escape tie.** The 95.902 keV line's Cd and Te escapes are in the fit, their ratio tied by the beam runs without Ge: Te Kα/Cd Kβ 0.353 ± 0.043 (rise < 300 ns; per-slice values as a variation) and Cd Kβ/Kα 0.324 against Geant4's 0.275 (`cdte/NoGeEscapes.cpp`).
- **The Pb Kβ groups** are in the in-run fit at their true energies; their sub-line intensities ± 10–20% move Ge 1.4 eV.
- **The move to `PulseHeightProcessor`** (2026-10-24) moved Ge +5.3 eV (from 68.7448), single slices by up to 26 eV. It is not a correction: the old pole-zero was the continuous approximation (a 1/(2τ) gain error and a tilt on the step), and the trigger, baseline window and peaking point read the same flat top. Event by event (`cdte/ReprocessingComparison.cpp`, five Ge runs) the new heights are the old times 0.9993–0.9995 per slice plus a random 76 eV, and the rise times differ by −0.3 ± 0.9 ns, which moves 7–15% of each slice's records into a neighbour (the 10–90% rise has a 2-sample sawtooth with peaks on the 280, 300 and 320 ns slice edges). The +5.4 eV splits into escape tie +0.4, old-only records +0.1, the same records reprocessed +3.0, new-only records +1.9; 104 pseudo-reprocessings of the old records give +3.2 ± 3.1 eV. The outlier term's growth (2.1 → 10.6 eV) is a rank flip between the two lowest slices, and the extrapolation term's (13.9 → 22.7 eV) follows the 280–290 ns slice's curvature-free fit, whose replicas scatter by ± 17 eV.
<!---->
## PSI Si
<!---->
- **Records at the bottom of the ADC range** are dropped (baseline RMS below 8 ADC): 10–11% of the records in the 10/11 Ge and 10/13 beam runs (under 1% elsewhere) have the raw baseline near 16 ADC instead of ~590, followed by a non-decaying ~650 ADC step at digitiser energy ~273 channels (~11.4 keV), 19–40 ms apart. They are not detector pulses.
- **The trapezoid against the digitiser's energy** is a choice, not a correction: on the same records the two differ by −4.9 ± 15.1 eV, less than their scatter.
