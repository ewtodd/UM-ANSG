# 73mGe — the corrections and why (2026-11-24)
<!---->
What each correction in the result does, how large it is, and the evidence for it; the current numbers and the plan are in `STATUS.md`.
A shift between two analysis choices is charged as a systematic only when more than one of the choices can be defended; fixing a defect, or changing an implementation detail that should give the same answer, is not a systematic (the CdTe reprocessing, below).
<!---->
## CZT: from 68.7502 to 68.7506 keV
<!---->
The value reported earlier, 68.7502 ± 0.0145 keV (c260617), moved through these steps (sizes where measured; the full table is for the write-up, Plan 7 in `STATUS.md`):
<!---->
1. **The Pb Kβ1,3 blend in the old pol2** was given Kβ1's energy, 84.936 keV, where the blend's centroid is about 84.770 (Kβ3 at 0.52 of Kβ1). That pulled Cu 01/13's Ge 27.4 eV below the pol1; with the centroid it sits 5.4 eV above. Production uses a pol1 through Am-241 and Pb Kα only.
2. **The RooFit fix**: ≤ 4 eV.
3. **Fit convergence**: up to 36 eV.
4. **The lineshape**: a Gaussian with a low exponential tail, τ/σ ≤ 20. The earlier τ/σ ≤ 8 railed at its bound (2026-10-22 tests), so that lineshape was set by the bound, not fitted.
5. **The Pb centroid link** between the background and signal runs.
6. **01/14 dropped**: the detector was unbiased on 01/13 and re-biased on 01/14, and there was no Am run that day, so no transfer from 01/13 applies.
7. **The signal-run Pb anchor** (below).
8. **The gain transfer model** (below).
9. **Both methods in the combination**: the gap between them was not significant (0.3–0.5σ), and charging it as a systematic used Method 2 only through that noise.
10. **The 95.902 keV escapes** (below): −7.9 eV on the quoted value.
11. **01/16 dropped** (2026-10-26): the graphite-castle run differs from the others in too many ways to be understood. Its residual runs from −5 counts/s at 62 keV to 0 at 71 keV, it is the highest-rate run (26–34k counts/s per crystal), no line appears in both its runs to check the source-run calibration, and its photons reach deeper into the crystal than in any other configuration. Its background run no longer enters the inside escape fraction either, which is now 01/15's alone. 68.7533 → 68.7511 keV: −1.8 eV for the run, −0.4 eV for the fraction.
12. **The day pixel map** (2026-10-27, below): 68.7511 → 68.7543 keV, +3.2 eV.
13. **The audit** (2026-10-27, below): the live time, the reference lines and smaller fixes, 68.7543 → 68.7551 keV, +0.8 eV.
14. **The Cd-shield escape source** as a systematic (2026-10-27, below): 68.7551 ± 0.0189 → 68.7567 ± 0.0191 keV, +1.6 eV through the BLUE weights, which move from the Cd rows to Cu.
15. **The escape fractions measured, every run in-crystal** (2026-11-24, below): the in-crystal fractions from an unbinned fit with Geant4's series ratios (Cd Kβ 0.0122 → 0.0140), the Cd-shield runs on them too, the escape-source systematic dropped, and the escape systematic from the fractions' own uncertainty instead of the full correction. 68.7567 ± 0.0191 → 68.7506 ± 0.0183 keV, −6.1 eV.
<!---->
How the anchor (`_sigpb`), the transfer model (`_oxfer`), the escapes (`_esc`) and the day pixel map (`_pxday`) move each run: Ge − 68.752 keV per in-situ row in eV; the rate-sub row (01/15) is +23.9 eV without the escapes and +13.9 eV with them, +29.0 and +20.3 eV with the pixel map, +27.6 and +25.1 eV after the audit. Each tag combined alone, without 01/16; the `_pxday` rows are after the audit.
<!---->
| Tag | Cu 01/13 | Cd 10% | Cd 25% | Method 1 | Both methods, BLUE | Outlier | Quoted |
|---|---|---|---|---|---|---|---|
| `tau20` | +26.8 | +1.7 | +6.4 | 68.7695 ± 0.0171 | 68.7701 ± 0.0163 | 10.3 | 68.7701 ± 0.0192 |
| `tau20_oxfer` | +26.8 | +12.8 | +17.8 | 68.7742 ± 0.0172 | 68.7744 ± 0.0163 | 5.4 | 68.7744 ± 0.0172 |
| `tau20_sigpb` | +22.4 | −31.9 | +9.9 | 68.7580 ± 0.0168 | 68.7596 ± 0.0158 | 12.4 | 68.7596 ± 0.0200 |
| `tau20_sigpb_oxfer` | +22.4 | −22.0 | +21.5 | 68.7627 ± 0.0165 | 68.7639 ± 0.0158 | 10.5 | 68.7639 ± 0.0189 |
| `tau20_sigpb_esc` | +10.1 | −32.7 | +9.1 | 68.7507 ± 0.0166 | 68.7521 ± 0.0158 | 10.2 | 68.7521 ± 0.0188 |
| `tau20_sigpb_esc_oxfer` | +10.1 | −22.6 | +20.9 | 68.7553 ± 0.0166 | 68.7563 ± 0.0158 | 8.5 | 68.7563 ± 0.0179 |
| `tau20_sigpb_esc_pxday` | +14.0 | −21.8 | −8.5 | 68.7529 ± 0.0165 | 68.7556 ± 0.0156 | 7.4 | 68.7556 ± 0.0173 |
| `tau20_sigpb_esc_oxfer_pxday` | +14.0 | −12.2 | +3.3 | 68.7574 ± 0.0165 | 68.7596 ± 0.0156 | 5.8 | 68.7596 ± 0.0167 |
<!---->
The table predates step 15. The quoted 68.7506 keV is `tau20_sigpb_esc_pxday` combined with its transfer-model and escape-fraction partners as correlated terms (`STATUS.md`).
<!---->
## CZT corrections
<!---->
### Day pixel map
<!---->
The pixel gains are one table: 558.456 keV over each pixel's 558 keV centroid in the 01/16 graphite-castle run, at every depth. Near 70–96 keV, after those gains, the pixels do not agree (`czt/PixelGainStability.cpp`, plots in `plots/czt/pixel_gain`). Each pixel's line position against its run's summed spectrum, in eV at 68.75 keV:
<!---->
- **The crystals read 80–300 eV apart.** Crystal means (x < 0, y < 0 / x > 0, y < 0 / x < 0, y > 0 / x > 0, y > 0): Pb Kα on 01/13 +63 / +22 / −81 / +8; the 95.902 keV line on 01/13 +140 / +53 / −138 / +22 and on 01/15 +131 / +40 / −154 / +16. At 558 keV the four agree within about ±10 eV.
- **The pixels within a crystal spread by about 110 eV rms** beyond their errors.
- **The 558 keV gains do not see this pattern:** on the same day the 95.902 keV and 558 keV maps correlate at 0.05 (01/13), −0.21 (01/14), 0.07 (01/15) and 0.00 (01/16).
- **It is mostly fixed.** The 95.902 keV maps of different days correlate at 0.65–0.93 (01/15 against 01/16: 0.93), and differ by 40–115 eV rms beyond their errors. The 558 keV pattern changes by 7–27 eV.
<!---->
A single run measures a pixel to 90–300 eV, about the size of the pattern; the day's beam runs together measure it to about 50 eV (01/13, eight runs) and 28 eV (01/15). So every event is aligned on its day's map: Pb Kα on 01/13, behind the Cu or Cd shield, where it stands twice the continuum; the 95.902 keV line on 01/15, where Pb Kα stands 2–5% above it. Both come from the same filtered events; the 95.902 keV events sit closer to the cathode (61–66% within 0.5 mm of it on 01/13–15, against 43–49% for Pb Kα), and on the same day its map and Pb Kα's correlate at 0.45. A pixel measured to 100 eV or better gets the factor μ_run/μ_pixel, which leaves each run's summed scale where it was; the source runs take their day's map.
<!---->
- **Effect.** Cu +3.8, Cd 10% +10.9, Cd 25% −18.2, rate-sub +6.4 eV (its fit error 42.7 → 37.5 eV, χ²/ndf 1.43 → 1.12); the quoted value 68.7511 → 68.7543 keV, its outlier term 11.2 → 7.1 eV.
- **Alternatives.** Each run on its own line moved the quoted value +6.9 eV, each signal run on its background run's −2.9 eV; both measure the pixels from one run, to 90–300 eV. A calibration per crystal (Am-241 and Pb Kα on 01/13, Am-241 and Ba-133 on 01/15) moved it +4.6 eV but fitted the in-situ pairs worse (χ²/ndf up by 0.1–0.4); the Am source runs light the crystals differently from the beam, most clearly on 01/15.
<!---->
### Signal-run Pb anchor
<!---->
Each 01/13 signal run is shifted onto its own Pb Kα, so it inherits only its pair's slope.
<!---->
- Ge follows its own run's Pb: in template fits Ge − Pb within a signal run is the same in the Cd and Cu runs to ±45 eV while Pb moves up to 213 eV.
- The Cd 10% signal run reads +36 ± 11 eV above its background run within the same crystal, pixel ring and depth slice (where the Pb lands contributes +0.2 eV): a shift of the whole run. Cu 01/13 and Cd 25% show none (+5 ± 12, +11 ± 25 eV). The cause is unknown.
- Pb Kβ, 10 keV above the window and not used by the anchor, follows: on the anchored calibration the Cd 10% signal run's Kβ reads +46 ± 30 eV against its background run's, where a shift of the Kα lines alone would leave it at −31 eV (2.6σ away) and a shift of the whole run at 0 (1.5σ). Cu 01/13 and Cd 25% read −28 ± 37 and +30 ± 77 eV, Kα itself 0 to +7 eV in all three (comb-period bins at 20 starting points; plots in `plots/czt/studies/pb_kbeta_check`).
- So Cd 10%'s low Ge, −35 eV from the mean and −1.0σ, is not an artefact of the anchor.
<!---->
### Gain transfer to the Cd pairs
<!---->
The Cd pairs inherit the Cu 10% calibration through their background runs' Pb Kα, as a gain (production) or as an offset (`TRANSFER=offset`). Pb Kβ lies 10 keV above the calibrated range and the Cd background runs' in-run Kα1 → Kβ1 slopes split between the two; no line inside the range resolves the 0.15% slope difference. Both can be defended, so the difference (+10.0 and +12.1 eV on the Cd rows) is a correlated systematic in the combination.
<!---->
### 95.902 keV escapes
<!---->
- **The line** is 113Cd(n,γ)'s E2 transition in 114Cd, converted in the K shell 47% of the time and emitted as a photon 37% (BrIcc: αK 1.288, αT 1.735). Its escapes come from the Cd K X-rays of conversions at captures near a surface, and from Cd and Te K X-rays when its photon is absorbed: Cd Kα at 72.73 / 72.92 keV under Pb Kα2, Cd Kβ at 69.25–69.81 keV beside Ge, Te Kα at 68.43 / 68.70 keV under Ge, Te Kβ at 64.2–64.9 keV.
- **Captures inside the crystal, every run** (2026-11-24). `czt/CaptureEscapeFraction.cpp` fits the 01/15 beam background run, which has neither a shield nor Pb, unbinned over 66–74.5 keV (1.29 million events): both Cd escape pairs as one yield on a linear continuum, Cd Kα held to Cd Kβ (3.671) and the Te escapes to Cd Kα (0.0534, then 0.300) by Geant4's ratios for captures in the crystal (`czt96_surface`, printed by `geant4/Escape96Plots.cpp`), at the line's centroid minus the SRD 128 X-ray energies, with the line's width (0.837 keV) and the tail interpolated between the day's Am-241 and Ba-133 fits. Over the line's Gaussian yield (production's capture-line fit, 225,065 ± 1,254 counts): Cd Kβ 0.01405 ± 0.0016 (statistics) ± 0.0008 (the largest shift over the photopeak width, 0.810 keV, and a quadratic continuum), 12.5%; the set 0.05158, 0.01405, 0.00275 and 0.00083 (`results/escape_fractions.result`, `ESCAPES_INSIDE`). Fitted apart, Cd Kα reads 0.0516 ± 0.0058 and Cd Kβ 0.0143 ± 0.0042. Pb Kα1 is absent (−2.3k ± 2.2k counts against 12k in the Cd Kα escapes), so no Pb Kα2 lies under them. Binned fits of the same window move Cd Kβ by 0.001 rms with where the bins start; the unbinned fit has none. Tied instead by the Table of Isotopes intensities times slab escape probabilities (Kα/Kβ 4.187), the fit gives 0.0122, the value used until 2026-11-24 from a script that no longer exists; nothing reproduces those slab probabilities, so Geant4's ratio is used (the user's call).
- **Photons from the Cd shield** (Geant4 `czt96_beam`, the line's photons onto the face): 0.0115, 0.0033, 0.0057 and 0.0017. No run uses them. The Cd-shield background runs cannot tell them from the in-crystal set: their 95.902 keV lines hold 6,947 and 9,569 counts, their Cd Kβ escape fitted beside free Pb Kα comes out to ± 0.03–0.09, and the two complete sets fit alike (Baker–Cousins χ² of the in-crystal set minus the shield's −0.4 and −1.1). The "no Cd Kβ escape, 0 ± 0.004" once quoted here came from the same lost script and does not reproduce.
- **Where the line is detected** (`czt/Capture96Positions.cpp`, single-interaction events, `plots/czt/studies/capture96_positions`): within 1 mm of the cathode for 82% (Cu), 87% and 85% (Cd 10% and 25%) and 84% (01/15) of the line's events, against 72–77% for Pb Kβ in the same runs and 67% for Am-241 entering the cathode, with no excess at the side faces (the edge-ring profile is alike in every configuration). So the line is mostly captures in the CZT's own cathode layer, in the Cd-foil runs as well: the foil covered the four sides and left the cathode face (toward the Ge) and the anode face open. The whole profiles, compared bin by bin (χ² shape tests, the macro's `SHAPE` lines): the Cd runs' depth profiles differ from the Cu run's (79/37 and 50/37) and Cd 25%'s ring profile too (22/5), but the Cu and no-shield runs, neither behind Cd, differ as much (83/37 in depth, 39/5 in rings), and the share in the last 0.25 mm (0.43, 0.46, 0.43, 0.42) and the edge ring show no Cd pattern (Cd 25%'s edge ring is the lowest). So every run takes the in-crystal fractions. In 01/16 the line (62%) and Pb Kβ (45%) reach deeper.
- **The 558 keV line** (114Cd's 2+ → 0+, same macro), which crosses the crystal, against two lines from outside it at a similar energy in the same runs: 511 keV, and 595.8 keV from captures in the Ge target, entering through the cathode face. Within 1 mm of the cathode lie 0.357 (Cu), 0.311 and 0.331 (Cd 10% and 25%) and 0.374 (01/15) of the 558 keV events, against 0.20–0.25 for the two outside lines in every run, so captures in the cathode's surface layer happen in every configuration, the Cd-foil runs included. The Cd runs' cathode excess is 0.71–0.82 of the Cu run's and their 558 keV line is 1.24 and 1.08 times as strong relative to 95.9 keV, so part of their 558 keV (about a fifth at Cd 10%) comes from outside, though no edge-ring excess shows where it enters. Scaled by its inside share, the 95.9 keV line is inside-type to 0.88–0.92 (Cd 10%), 0.89 (Cd 25%) and 0.85–0.92 (01/15) of the Cu run's.
- **In the fits** the escapes are fixed lines (`RooFitUtils::AddFixedLines`) at the run's own 95.902 keV centroid minus the X-ray energy, with that line's Gaussian yield times the fractions. In the in-situ fits each channel gets its own run's lines, drawn with the Pb Kα1 shape; the line's yield and centroid come from its fits in comb-period bins at five starting points over 91.5–100.5 keV, with the in-situ fit's converged tail (the raw fit first runs once without the escapes). In the rate-sub fit they are a fixed term with Ge's width and tail, scaled by the residual's own 95.902 keV line, which is 1.08 times stronger in the 01/15 signal run than live time predicts (1.2 with the old live times).
- **Effect** (production − `ESCAPES=off`, printed, not charged): Cu 01/13 −14.1 eV, Cd 10% −15.0, Cd 25% −16.4, 01/15 −3.3 eV.
- **Systematic.** Production rerun with every fraction moved up and down by 12.5% (`ESCAPE_FRACTIONS=up`, `down`); half the difference, correlated between all runs: Cu −1.9, Cd −1.7 and −2.0, 01/15 −0.4 eV, 1.7 eV on the quoted value. The full correction is no longer charged: the escapes are measured and plainly visible (`plots/czt/figures/czt_escapes_0115`), so "no escapes" is not a defensible alternative. The ratios between the series are Geant4's and carry no uncertainty of their own; the data cannot fix Te Kα/Cd Kα (fitted free, the Te Kα fraction reads −0.002 ± 0.004), and there is no second value to set against it.
- **Validation.** The fixed lines reproduce a testbed's doublet without bias (−0.1 ± 2.1 and −0.3 ± 1.4 eV; +12.4 and +7.6 eV without them), and with the escapes switched off production reproduces `tau20_sigpb` exactly. In production's own fit, on a synthetic replica of Cu 01/13 with the escapes in the data, the fit with the escape lines recovers the injected 68.750 keV (+1.8 ± 2.7 eV, 50 trials); the same events fitted without them read +11.1 ± 0.1 eV higher in every trial, the size of production's Cu correction.
<!---->
### Rate subtraction (Method 2)
<!---->
- **Without the subtraction** the fits read +91 (01/15) and +134 eV (01/16), and a template with a free scale degenerates (it fits 0.3–1.3 against live-time ratios of 3.6 and 1.8): the subtraction removes real structure near Ge.
- **The energy comb.** Every pixel's energies sit on their own 1–8 eV grid; summed, they show a comb with a 264 eV period. Bins near that period beat against it: 255 eV bins move rate-sub Ge by 13–16 eV rms with where they start, comb-period bins by 8–10 eV, hence comb-period bins at five starting points, averaged, with their spread in the systematic.
- **Overdispersion.** The signal and background runs fill the grids differently, so signal − α × background is 7–16 times Poisson at 0.5–10 eV bins, about 4 at 50 eV, 1.4–1.9 at 255 eV; unbinned or fine-binned two-run fits fail, and from ~20 eV bins up an on/off likelihood agrees with the subtraction (plots in `plots/czt/studies/method2_unbinned`). The fit error is scaled by √(χ²/ndf).
- **Live time.** Each crystal's rate is its counts over its own live time: the `liveTime` of all its events (the live interval since its previous one), summed in double before any cut. Against the seconds that hold data the live fraction falls with rate at about 21 µs per event in every run (0.998 at 120 events/s per crystal, 0.888 and 0.870 for the 01/15 background and signal runs at 5.1 and 6.0 kHz, 0.61 at 17 kHz); the N42 real time does not cover the same stretch of data and is not used.
- **Source runs to beam runs.** 01/15 is calibrated on source runs and no line appears in both, so that step cannot be checked there; no uncertainty is charged for it.
<!---->
### Model choices not charged
<!---->
- **The lineshape**: one shared low tail on the Pb lines and Ge. Every line here interacts at the cathode surface, so a common tail is the one defensible choice.
- **The in-situ continuum** stays flat. A sloped one takes up shape that belongs to the lines' tails: with it, two of the three pairs end with a forced covariance (Cu post-cal, Cd 10% pre-cal), and the spectra look flat across the window.
- **The Cd pairs across the Cd → Cu step** (130 eV): the agreement across that step, compared with 01/13 → 01/14, suggests the bias was not cycled there, or only briefly; the transfer systematic covers it.
<!---->
### Audit (2026-10-27)
<!---->
Four read-only reviews of the CZT chain (selection, calibration, fit core, Analysis-Utilities' fitting) found no defect in AU that biases Ge, and these in the CZT code:
<!---->
- **The rate-sub live time** (`czt/Filter.cpp`) was summed over accepted events only, so it carried each run's acceptance (0.96–0.98 between the 01/15 runs per crystal), and in Float_t, which drops 1e-4 s steps once the sum reaches thousands of seconds (01/15 signal 6–11% low). The background was scaled by 3.53–3.71 instead of 4.10–4.13, leaving 10–14% of it in the residual and most of the 95.902 keV line's apparent 1.2× excess. Method 2 68.7723 → 68.7771 keV: −1.3 eV without the escapes, the rest the smaller residual line's escapes (−8.6 → −2.6 eV).
- **The pile-up rule** compared with the previous accepted event, not the previous event: 0.11% of the events.
- **The reference lines** now come from one source per quantity (`STATUS.md`). Pb Kα were Bearden's values (72.8042, 74.9694 keV), 1.2 and 0.7 eV below SRD 128's; the Cd and Te Kβ escape groups were 6–12 eV high; the Cd Kα escape fraction became 0.0511 (since remeasured, above). With the pile-up rule they moved the in-situ rows +0.1, 0.0 and +0.7 eV.
- **The rate-sub bin-phase spread** came out 0: its rms was E² − mean² in Float_t at 68.8 keV. It is 8.2 eV. The 20 eV once charged for the source runs is gone (above).
- **The gain transfer's error** on the anchored Cd rows was set to 0; it acts on the slope over the Pb-to-Ge distance: 0.9 and 0.8 eV.
- **Guards.** A fit is used only from a clean minimum (Minuit status 0, covariance quality 3) with its centroid and width off their limits; a missing escape line, anchor, reference or live time stops the run instead of falling back to another model; a missing partner row stops the combination.
- **Method 1 against Method 2** was printed as independent; the shared escape fractions correlate them (29 ± 50 eV).
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
