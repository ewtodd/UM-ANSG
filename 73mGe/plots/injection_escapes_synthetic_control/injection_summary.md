# Injection (closure) test — escapes_synthetic_control

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cu 10% (01/13) | model | sim | 20 | 2.03 | 4.46 | 0.086 | 0.938 | 0.0 | 60.0 | 0.0 | -0.367 |
| Cu 10% (01/13) | model | sim_noesc | 20 | 2.03 | 4.46 | 0.086 | 0.938 | 0.0 | 60.0 | 0.0 | -0.367 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cu 10% (01/13) | model | sim | 68.750 | 20 | 2.03 | 4.46 | 0.086 | 0.938 | 0.0 | 60.0 | 0.0 | -0.367 |
| Cu 10% (01/13) | model | sim_noesc | 68.750 | 20 | 2.03 | 4.46 | 0.086 | 0.938 | 0.0 | 60.0 | 0.0 | -0.367 |

## Escape effect (fit without the escape lines minus with, paired by trial)

| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |
|---|---|---|---|---|---|---|
| Cu 10% (01/13) | model | 68.750 | 20 | 0.00 | 0.00 | 0.00 |

![escape bias by pair](escape_bias_by_pair.png)

## Interpretation


### model / sim

![bias vs energy](bias_vs_energy_model_sim.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 2.03 +/- 4.46 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

### model / sim_noesc

![bias vs energy](bias_vs_energy_model_sim_noesc.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| Cu 10% (01/13) | 2.03 +/- 4.46 eV | 0.0000 +/- 0.0000 | 0.000 | 0.000 | insufficient energy points |

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
