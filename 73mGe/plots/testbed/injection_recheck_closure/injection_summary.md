# Injection (closure) test — recheck_closure

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | ratesub | 3 | 38.81 | 9.54 | 0.666 | 0.703 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cd 10% (01/13) | model | sim | 3 | 112.89 | 27.65 | 2.545 | 2.674 | 0.0 | 100.0 | 0.0 | -0.149 |
| Cd 10% (01/13) | model | sim_noesc | 3 | 112.65 | 27.56 | 2.649 | 2.796 | 0.0 | 100.0 | 0.0 | -0.150 |
| Cd 25% (01/13) | model | ratesub | 3 | -9.17 | 35.19 | -0.369 | 1.307 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cd 25% (01/13) | model | sim | 3 | 18.48 | 26.53 | 0.781 | 1.653 | 0.0 | 100.0 | 0.0 | 0.843 |
| Cd 25% (01/13) | model | sim_noesc | 3 | 18.65 | 26.38 | 0.835 | 1.690 | 0.0 | 100.0 | 0.0 | 0.844 |
| Cu 10% (01/13) | model | ratesub | 3 | -24.81 | 6.82 | -0.349 | 0.372 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cu 10% (01/13) | model | sim | 3 | 10.90 | 9.46 | 0.229 | 0.366 | 0.0 | 100.0 | 0.0 | 0.251 |
| Cu 10% (01/13) | model | sim_noesc | 3 | 13.90 | 8.89 | 0.276 | 0.371 | 0.0 | 100.0 | 0.0 | 0.215 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | ratesub | 68.750 | 3 | 38.81 | 9.54 | 0.666 | 0.703 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cd 10% (01/13) | model | sim | 68.750 | 3 | 112.89 | 27.65 | 2.545 | 2.674 | 0.0 | 100.0 | 0.0 | -0.149 |
| Cd 10% (01/13) | model | sim_noesc | 68.750 | 3 | 112.65 | 27.56 | 2.649 | 2.796 | 0.0 | 100.0 | 0.0 | -0.150 |
| Cd 25% (01/13) | model | ratesub | 68.750 | 3 | -9.17 | 35.19 | -0.369 | 1.307 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cd 25% (01/13) | model | sim | 68.750 | 3 | 18.48 | 26.53 | 0.781 | 1.653 | 0.0 | 100.0 | 0.0 | 0.843 |
| Cd 25% (01/13) | model | sim_noesc | 68.750 | 3 | 18.65 | 26.38 | 0.835 | 1.690 | 0.0 | 100.0 | 0.0 | 0.844 |
| Cu 10% (01/13) | model | ratesub | 68.750 | 3 | -24.81 | 6.82 | -0.349 | 0.372 | 0.0 | 100.0 | 0.0 | 0.000 |
| Cu 10% (01/13) | model | sim | 68.750 | 3 | 10.90 | 9.46 | 0.229 | 0.366 | 0.0 | 100.0 | 0.0 | 0.251 |
| Cu 10% (01/13) | model | sim_noesc | 68.750 | 3 | 13.90 | 8.89 | 0.276 | 0.371 | 0.0 | 100.0 | 0.0 | 0.215 |

## Escape effect (fit without the escape lines minus with, paired by trial)

| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |
|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | 68.750 | 3 | -0.23 | 0.11 | 0.15 |
| Cd 25% (01/13) | model | 68.750 | 3 | 0.17 | 0.15 | 0.21 |
| Cu 10% (01/13) | model | 68.750 | 3 | 3.00 | 0.66 | 0.94 |

![escape bias by pair](escape_bias_by_pair.png)

## Interpretation


### model / ratesub

![bias vs energy](bias_vs_energy_model_ratesub.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | -24.81 +/- 6.82 eV | 0.2500 +/- 6958505.3814 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 38.81 +/- 9.54 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | -9.17 +/- 35.19 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim

![bias vs energy](bias_vs_energy_model_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 10.90 +/- 9.46 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 112.89 +/- 27.65 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | 18.48 +/- 26.53 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim_noesc

![bias vs energy](bias_vs_energy_model_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 13.90 +/- 8.89 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 112.65 +/- 27.56 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | 18.65 +/- 26.38 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### empirical / ratesub

![bias vs energy](bias_vs_energy_empirical_ratesub.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|

### empirical / sim

![bias vs energy](bias_vs_energy_empirical_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|

### empirical / sim_noesc

![bias vs energy](bias_vs_energy_empirical_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|

## Pre-fixed interpretation

- Zero-consistent bias: unbiased; quote the bias uncertainty as the systematic.
- Significant energy-independent bias: correct and propagate.
- Energy-dependent bias (line slope > 2 sigma): stop and report; do not correct.
