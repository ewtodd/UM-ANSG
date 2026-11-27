# Injection (closure) test — final_insitu_synthetic_halves

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 80 | 4.57 | 4.56 | 0.119 | 1.056 | 0.0 | 100.0 | 0.0 | -0.012 |
| Cd 10% (01/13) | model | sim_noesc | 79 | 6.12 | 4.50 | 0.158 | 1.043 | 1.2 | 100.0 | 0.0 | -0.019 |
| Cd 25% (01/13) | model | sim | 80 | -7.28 | 3.33 | -0.302 | 1.216 | 0.0 | 100.0 | 0.0 | -0.340 |
| Cd 25% (01/13) | model | sim_noesc | 80 | -5.17 | 3.32 | -0.219 | 1.197 | 0.0 | 100.0 | 0.0 | -0.340 |
| Cu 10% (01/13) | model | sim | 80 | 10.12 | 5.45 | 0.229 | 1.111 | 0.0 | 100.0 | 0.0 | 0.045 |
| Cu 10% (01/13) | model | sim_noesc | 80 | 12.99 | 5.42 | 0.295 | 1.127 | 0.0 | 100.0 | 0.0 | 0.045 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 68.750 | 80 | 4.57 | 4.56 | 0.119 | 1.056 | 0.0 | 100.0 | 0.0 | -0.012 |
| Cd 10% (01/13) | model | sim_noesc | 68.750 | 79 | 6.12 | 4.50 | 0.158 | 1.043 | 1.2 | 100.0 | 0.0 | -0.019 |
| Cd 25% (01/13) | model | sim | 68.750 | 80 | -7.28 | 3.33 | -0.302 | 1.216 | 0.0 | 100.0 | 0.0 | -0.340 |
| Cd 25% (01/13) | model | sim_noesc | 68.750 | 80 | -5.17 | 3.32 | -0.219 | 1.197 | 0.0 | 100.0 | 0.0 | -0.340 |
| Cu 10% (01/13) | model | sim | 68.750 | 80 | 10.12 | 5.45 | 0.229 | 1.111 | 0.0 | 100.0 | 0.0 | 0.045 |
| Cu 10% (01/13) | model | sim_noesc | 68.750 | 80 | 12.99 | 5.42 | 0.295 | 1.127 | 0.0 | 100.0 | 0.0 | 0.045 |

## Escape effect (fit without the escape lines minus with, paired by trial)

| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |
|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | 68.750 | 10 | 2.36 | 0.11 | 0.32 |
| Cd 25% (01/13) | model | 68.750 | 10 | 2.23 | 0.09 | 0.27 |
| Cu 10% (01/13) | model | 68.750 | 10 | 2.84 | 0.23 | 0.68 |

![escape bias by pair](escape_bias_by_pair.png)

## Interpretation


### model / sim

![bias vs energy](bias_vs_energy_model_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 10.12 +/- 5.45 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 4.57 +/- 4.56 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | -7.28 +/- 3.33 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim_noesc

![bias vs energy](bias_vs_energy_model_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 12.99 +/- 5.42 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 6.12 +/- 4.50 eV | 0.0000 +/- 5271766.6551 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | -5.17 +/- 3.32 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

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
