# Injection (closure) test — final_direct_synthetic

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| No shield 5% (01/15) | model | direct | 181 | 122.72 | 87.99 | 0.938 | 11.317 | 4.7 | 100.0 | 0.5 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 60 | 177.01 | 182.03 | 0.902 | 9.009 | 68.4 | 100.0 | 0.0 | 0.000 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| No shield 5% (01/15) | model | direct | 66.500 | 2 | 11288.22 | 202.26 | 95.811 | 107.282 | 80.0 | 100.0 | 10.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.650 | 10 | -1.06 | 4.81 | -0.061 | 0.886 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.700 | 10 | 1.67 | 5.28 | 0.095 | 0.964 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.725 | 10 | -0.93 | 5.48 | -0.054 | 1.005 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.750 | 110 | -1.93 | 1.37 | -0.117 | 0.878 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.775 | 10 | -1.19 | 6.84 | -0.076 | 1.240 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.800 | 10 | -11.30 | 3.82 | -0.690 | 0.980 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.850 | 10 | 0.45 | 3.95 | 0.032 | 0.719 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 70.500 | 9 | -3.03 | 9.11 | -0.149 | 1.437 | 10.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 66.500 | 1 | 10915.80 | 0.00 | 69.380 | 69.380 | 90.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.650 | 0 | 0.00 | 0.00 | 0.000 | 0.000 | 100.0 | 0.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.700 | 3 | -12.61 | 10.60 | -0.645 | 1.011 | 70.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.725 | 3 | -19.24 | 13.59 | -0.900 | 1.222 | 70.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.750 | 35 | -0.21 | 2.95 | -0.010 | 0.892 | 68.2 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.775 | 3 | -7.03 | 21.55 | -0.376 | 1.621 | 70.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.800 | 4 | -19.44 | 4.07 | -1.080 | 1.152 | 60.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.850 | 7 | -8.00 | 6.20 | -0.412 | 0.868 | 30.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 70.500 | 4 | -9.35 | 6.16 | -0.482 | 0.735 | 60.0 | 100.0 | 0.0 | 0.000 |

## Interpretation


### model / direct

![bias vs energy](bias_vs_energy_model_direct.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| No shield 5% (01/15) | -1.94 +/- 1.11 eV | -17.8078 +/- 5.1508 | 390.368 | 444.427 | ENERGY-DEPENDENT (non-linear): stop and report; do not correct |

### model / direct_free_tail

![bias vs energy](bias_vs_energy_model_direct_free_tail.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| No shield 5% (01/15) | 8774.10 +/- 0.90 eV | -4191.5496 +/- 0.8919 | 3342130.956 | 218318.776 | ENERGY-DEPENDENT (non-linear): stop and report; do not correct |

### empirical / direct

![bias vs energy](bias_vs_energy_empirical_direct.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|

### empirical / direct_free_tail

![bias vs energy](bias_vs_energy_empirical_direct_free_tail.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|

## Pre-fixed interpretation

- Zero-consistent bias: unbiased; quote the bias uncertainty as the systematic.
- Significant energy-independent bias: correct and propagate.
- Energy-dependent bias (line slope > 2 sigma): stop and report; do not correct.
