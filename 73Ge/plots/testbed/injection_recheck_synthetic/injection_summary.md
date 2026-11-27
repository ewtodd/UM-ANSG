# Injection (closure) test — recheck_synthetic

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 3 | 4.46 | 16.09 | 0.075 | 0.749 | 0.0 | 100.0 | 0.0 | -0.885 |
| Cd 10% (01/13) | model | sim_noesc | 3 | 4.79 | 15.98 | 0.113 | 0.776 | 0.0 | 100.0 | 0.0 | -0.887 |
| Cd 25% (01/13) | model | sim | 3 | 20.12 | 43.86 | 0.588 | 1.795 | 0.0 | 100.0 | 0.0 | -0.925 |
| Cd 25% (01/13) | model | sim_noesc | 3 | 20.02 | 43.42 | 0.589 | 1.815 | 0.0 | 100.0 | 0.0 | -0.926 |
| Cu 10% (01/13) | model | sim | 3 | -9.62 | 14.88 | -0.215 | 0.966 | 0.0 | 100.0 | 0.0 | 0.237 |
| Cu 10% (01/13) | model | sim_noesc | 2 | 11.30 | 20.01 | 0.573 | 1.119 | 33.3 | 100.0 | 0.0 | 0.000 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | sim | 68.750 | 3 | 4.46 | 16.09 | 0.075 | 0.749 | 0.0 | 100.0 | 0.0 | -0.885 |
| Cd 10% (01/13) | model | sim_noesc | 68.750 | 3 | 4.79 | 15.98 | 0.113 | 0.776 | 0.0 | 100.0 | 0.0 | -0.887 |
| Cd 25% (01/13) | model | sim | 68.750 | 3 | 20.12 | 43.86 | 0.588 | 1.795 | 0.0 | 100.0 | 0.0 | -0.925 |
| Cd 25% (01/13) | model | sim_noesc | 68.750 | 3 | 20.02 | 43.42 | 0.589 | 1.815 | 0.0 | 100.0 | 0.0 | -0.926 |
| Cu 10% (01/13) | model | sim | 68.750 | 3 | -9.62 | 14.88 | -0.215 | 0.966 | 0.0 | 100.0 | 0.0 | 0.237 |
| Cu 10% (01/13) | model | sim_noesc | 68.750 | 2 | 11.30 | 20.01 | 0.573 | 1.119 | 33.3 | 100.0 | 0.0 | 0.000 |

## Escape effect (fit without the escape lines minus with, paired by trial)

| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |
|---|---|---|---|---|---|---|
| Cd 10% (01/13) | model | 68.750 | 3 | 0.33 | 0.13 | 0.18 |
| Cd 25% (01/13) | model | 68.750 | 3 | -0.10 | 0.44 | 0.62 |
| Cu 10% (01/13) | model | 68.750 | 2 | 11.98 | 0.61 | 0.61 |

![escape bias by pair](escape_bias_by_pair.png)

## Interpretation


### model / sim

![bias vs energy](bias_vs_energy_model_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | -9.62 +/- 14.88 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 4.46 +/- 16.09 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | 20.12 +/- 43.86 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim_noesc

![bias vs energy](bias_vs_energy_model_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 11.30 +/- 20.01 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |
| Cd 10% (01/13) | 4.79 +/- 15.98 eV | 0.0000 +/- 16797706.1840 | 0.000 | 0.000 | insufficient energy points |
| Cd 25% (01/13) | 20.02 +/- 43.42 eV | 0.0000 +/- 49454472.9362 | 0.000 | 0.000 | insufficient energy points |

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
