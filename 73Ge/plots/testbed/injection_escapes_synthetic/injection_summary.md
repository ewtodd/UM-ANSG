# Injection (closure) test — escapes_synthetic

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 20 | -11.37 | 8.32 | -0.348 | 1.252 | 0.0 | 95.0 | 0.0 | 0.252 |
| Cd 10% (01/13) | model | sim_noesc | 20 | -10.95 | 8.22 | -0.356 | 1.239 | 0.0 | 95.0 | 0.0 | 0.250 |
| Cd 25% (01/13) | model | sim | 20 | -4.34 | 9.23 | -0.127 | 1.205 | 0.0 | 100.0 | 0.0 | -0.261 |
| Cd 25% (01/13) | model | sim_noesc | 20 | -4.27 | 9.13 | -0.126 | 1.202 | 0.0 | 100.0 | 0.0 | -0.261 |
| Cu 10% (01/13) | model | sim | 50 | 1.82 | 2.66 | 0.075 | 0.886 | 0.0 | 84.0 | 0.0 | -0.023 |
| Cu 10% (01/13) | model | sim_noesc | 50 | 12.94 | 2.64 | 0.623 | 1.091 | 0.0 | 98.0 | 0.0 | -0.020 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 68.750 | 20 | -11.37 | 8.32 | -0.348 | 1.252 | 0.0 | 95.0 | 0.0 | 0.252 |
| Cd 10% (01/13) | model | sim_noesc | 68.750 | 20 | -10.95 | 8.22 | -0.356 | 1.239 | 0.0 | 95.0 | 0.0 | 0.250 |
| Cd 25% (01/13) | model | sim | 68.750 | 20 | -4.34 | 9.23 | -0.127 | 1.205 | 0.0 | 100.0 | 0.0 | -0.261 |
| Cd 25% (01/13) | model | sim_noesc | 68.750 | 20 | -4.27 | 9.13 | -0.126 | 1.202 | 0.0 | 100.0 | 0.0 | -0.261 |
| Cu 10% (01/13) | model | sim | 68.750 | 50 | 1.82 | 2.66 | 0.075 | 0.886 | 0.0 | 84.0 | 0.0 | -0.023 |
| Cu 10% (01/13) | model | sim_noesc | 68.750 | 50 | 12.94 | 2.64 | 0.623 | 1.091 | 0.0 | 98.0 | 0.0 | -0.020 |

## Escape effect (fit without the escape lines minus with, paired by trial)

| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |
|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | 68.750 | 20 | 0.43 | 0.11 | 0.48 |
| Cd 25% (01/13) | model | 68.750 | 20 | 0.07 | 0.11 | 0.48 |
| Cu 10% (01/13) | model | 68.750 | 50 | 11.12 | 0.07 | 0.49 |

![escape bias by pair](escape_bias_by_pair.png)

## Interpretation


### model / sim

![bias vs energy](bias_vs_energy_model_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 1.82 +/- 2.66 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | -11.37 +/- 8.32 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | -4.34 +/- 9.23 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim_noesc

![bias vs energy](bias_vs_energy_model_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 12.94 +/- 2.64 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | -10.95 +/- 8.22 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | -4.27 +/- 9.13 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

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
