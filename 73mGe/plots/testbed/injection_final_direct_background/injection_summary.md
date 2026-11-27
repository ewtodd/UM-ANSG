# Injection (closure) test — final_direct_background

Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS should be ~1 if the fit errors are correct.

## Aggregate

| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| No shield 5% (01/15) | model | direct | 160 | 14.00 | 3.02 | 0.389 | 1.163 | 11.1 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 145 | 101.03 | 67.16 | 0.979 | 2.486 | 19.4 | 100.0 | 0.0 | 0.000 |

## Per injection energy

| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| No shield 5% (01/15) | model | direct | 66.500 | 1 | -22.70 | 0.00 | -0.782 | 0.782 | 95.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.650 | 20 | 15.89 | 8.62 | 0.441 | 1.171 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.700 | 20 | 15.74 | 7.70 | 0.446 | 1.062 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.725 | 20 | 22.10 | 5.00 | 0.627 | 0.900 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.750 | 20 | 23.46 | 8.09 | 0.671 | 1.193 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.775 | 20 | 22.12 | 8.25 | 0.626 | 1.179 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.800 | 20 | 16.13 | 9.77 | 0.460 | 1.336 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 68.850 | 20 | 22.97 | 8.40 | 0.643 | 1.199 | 0.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct | 70.500 | 19 | -26.64 | 7.49 | -0.803 | 1.230 | 5.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 66.500 | 5 | 1892.05 | 1964.18 | 3.507 | 10.794 | 75.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.650 | 19 | 29.75 | 12.93 | 0.861 | 1.367 | 5.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.700 | 18 | 48.19 | 7.72 | 1.150 | 1.395 | 10.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.725 | 16 | 30.61 | 13.84 | 0.787 | 1.494 | 20.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.750 | 16 | 53.41 | 16.67 | 1.280 | 2.053 | 20.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.775 | 18 | 56.72 | 11.36 | 1.266 | 1.624 | 10.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.800 | 17 | 39.75 | 15.55 | 0.862 | 1.647 | 15.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 68.850 | 17 | 34.79 | 16.20 | 0.805 | 1.421 | 15.0 | 100.0 | 0.0 | 0.000 |
| No shield 5% (01/15) | model | direct_free_tail | 70.500 | 19 | 6.54 | 6.83 | 0.168 | 0.784 | 5.0 | 100.0 | 0.0 | 0.000 |

## Interpretation


### model / direct

![bias vs energy](bias_vs_energy_model_direct.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| No shield 5% (01/15) | -18.05 +/- 0.94 eV | 12.6017 +/- 1.1138 | 26.034 | 11.467 | ENERGY-DEPENDENT (non-linear): stop and report; do not correct |

### model / direct_free_tail

![bias vs energy](bias_vs_energy_model_direct_free_tail.png)

| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | chi2/ndf (line) | Classification |
|---|---|---|---|---|---|
| No shield 5% (01/15) | 32.25 +/- 3.82 eV | -21.0018 +/- 4.6548 | 3.235 | 0.789 | ENERGY-DEPENDENT (non-linear): stop and report; do not correct |

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
