"""Noise robustness study for the trained XGBoost regressor.

Adds white (Gaussian) and shot (sqrt(signal) Gaussian-approximated) noise to
the held-out test waveforms before normalization and reports AUC. Runs two
experimental conditions:

1. **Test-only noisy** -- the production XGBoost model (trained on the
   unmodified waveforms) is evaluated on test data with extra noise added.
   Tests robustness of an already-deployed model to unforeseen experimental
   noise.
2. **Matched noisy** -- a fresh XGBoost is trained on noisy training data
   and evaluated on noisy test data at the same level. Tests the ceiling on
   what's achievable if the noise conditions were known at training time.

White-noise sigma is expressed in ADC counts (1 count = 1 LSB) so sub-LSB
and super-LSB regimes are directly readable. Shot-noise sigma at each
sample = ``multiplier * sqrt(|sample_ADC|)``, the Gaussian approximation to
per-sample Poisson statistics; the multiplier isolates the noise amplitude
from the natural sqrt(signal) shape.
"""
import os
import pickle
import numpy as np
import pandas as pd
import ROOT
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, iirnotch, sosfiltfilt, tf2sos
from sklearn.base import clone
from analysis_utilities.io import load_tree_data
from psd_utils import (process_waveforms, ANALYSIS_CACHE_DIR, ROOT_FILES_DIR,
                       per_lo_auc_from_scores, error_weighted_auc,
                       cumulative_shap, style_residual_pad_axes,
                       compute_mean_abs_shap)
from psd_utils import _get_training_indices
from regressors import get_default_regressors
import analysis_utilities

analysis_utilities.load_cpp_library()
ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPDF)

CACHE_DIR = "noise_study_cache"
MODEL_PATH = os.path.join(ANALYSIS_CACHE_DIR, "xgb_regressor.pkl")

# Cached XGBoost SHAP feature-importance arrays from parameter_study.py
# (one numpy array per random seed, mean |SHAP| per input feature).
SWEEP_CACHE_DIR = "sweep_cache"
SHAP_SEEDS = [42, 123, 256]
# Number of test events to compute SHAP on. Matches the n=1000 limit in
# parameter_study.py for consistency.
SHAP_N_EXPLAIN = 1000

# Per-sample variance vs mean polynomial fit, sigma^2(mu) = p0 + p1*mu +
# p2*mu^2, obtained interactively from variance_vs_mean.root via FitPanel.
# Used to subtract the (electronic floor + shot + multiplicative) trend
# from per-sample variance so the residual isolates jitter / other effects.
ALPHA_VAR_QUAD_P0 = 168
ALPHA_VAR_QUAD_P1 = 2.46797
ALPHA_VAR_QUAD_P2 = 0.0172492
GAMMA_VAR_QUAD_P0 = 1022
GAMMA_VAR_QUAD_P1 = 1.08675
GAMMA_VAR_QUAD_P2 = 0.0184931

# Approximate electronic-noise RMS measured on the raw waveforms (pre-trigger
# baseline standard deviation). Marked on the white-noise AUC plot as a
# reference for "extra noise equal to / above the existing noise floor".
BASELINE_RMS_ADC = 10.0
WHITE_SIGMAS_ADC = [0.0, 10.0, 15.0, 20.0, 30.0]
SHOT_MULTIPLIERS = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
NOISE_SEED = 42

# Filter sweeps for the train+test filtered comparison. Both sweeps share
# an effective -3dB cutoff x-axis for plotting; the parameter that's
# actually swept (and cached) is the filter's native knob.
#
# Sample rate: 500 MHz (2 ns/sample) -> Nyquist = 250 MHz.
SAMPLE_RATE_HZ = 500e6
LPF_ORDER = 4
LPF_CUTOFFS_MHZ = [25.0, 50.0, 75.0, 100.0, 150.0, 249.0]
# Gaussian smoothing sigma in samples. Effective -3dB cutoff f =
# sqrt(ln 2) / (2 pi sigma_t), sigma_t = sigma_samples * 2 ns. So
# sigma_samples = 1 -> ~66 MHz, 2 -> ~33 MHz, 5 -> ~13 MHz.
GAUSSIAN_SIGMAS_SAMPLES = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

# Discrete coherent-noise peaks identified in the per-event residual PSD
# (likely sample-clock / PLL feedthrough: fs/5 = 100 MHz, fs/4 = 125 MHz,
# fs/2.5 = 200 MHz at 500 MS/s). Notch experiment retrains XGBoost with
# these frequencies suppressed to test whether the model is exploiting
# class-asymmetric hardware artifacts vs. true scintillation physics.
NOTCH_FREQS_MHZ = [100.0, 125.0, 200.0]
NOTCH_Q = 30.0

# Balanced subsampling to avoid class-imbalance effects in AUC computation.
N_TRAIN_PER_CLASS = 10000
N_TEST_PER_CLASS = 100000
SUBSAMPLE_SEED = 123

# Long-gate integral parameters matching the C++ pipeline
# (Analysis-Utilities/src/WaveformProcessingUtils.cpp + YAP-PSD Constants.hpp).
# Stored waveforms are already baseline-subtracted and polarity-inverted, so
# the long_integral feature is reproduced by summing samples from
# INTEGRATION_START to INTEGRATION_START+LONG_GATE clamped to the array end.
PRE_SAMPLES = 20
PRE_GATE = 5
LONG_GATE = 250
INTEGRATION_START = PRE_SAMPLES - PRE_GATE

# Restrict both classes to a common light-output window so the alpha/gamma
# samples have comparable peak amplitudes. This isolates the noise effects
# from the peak-amplitude-distribution confound introduced by max-normalization
# (without this restriction, low-LO gammas get noise amplified relative to
# alphas simply because they get divided by a smaller peak during
# process_waveforms).
LO_LOWER = {"alpha": 375, "gamma": 0}
LO_UPPER = {"alpha": 1575, "gamma": 1750}
LO_TAG = f"{LO_LOWER['alpha']}_{LO_UPPER['alpha']}"

RAW_DATA_CACHE = os.path.join(CACHE_DIR, f"raw_split_{LO_TAG}.npz")


def _add_white_noise(wf, sigma_adc, rng):
    if sigma_adc <= 0:
        return wf
    return wf + rng.normal(0.0, sigma_adc, size=wf.shape)


def _add_shot_noise(wf, multiplier, rng):
    if multiplier <= 0:
        return wf
    sigma = multiplier * np.sqrt(np.abs(wf))
    return wf + rng.normal(size=wf.shape) * sigma


def _balanced_subsample(arr, n, seed):
    """Pick ``n`` events uniformly at random from ``arr`` without replacement.
    Returns the full array unchanged if it already has <= n events."""
    if n >= len(arr):
        return arr
    rng = np.random.RandomState(seed)
    idx = rng.choice(len(arr), size=n, replace=False)
    return arr[idx]


def _balanced_subsample_pair(arr_a, arr_b, n, seed):
    """Subsample two aligned arrays with the same random indices."""
    if n >= len(arr_a):
        return arr_a, arr_b
    rng = np.random.RandomState(seed)
    idx = rng.choice(len(arr_a), size=n, replace=False)
    return arr_a[idx], arr_b[idx]


def _make_xgb(seed=None):
    """Return an unfitted XGBoost regressor matching regressors.py. If seed
    is given, set random_state for reproducible per-seed training."""
    for r in get_default_regressors():
        if r["name"] == "XGBoost":
            model = clone(r["model"])
            if seed is not None:
                model.set_params(random_state=seed)
            return model
    raise KeyError("XGBoost not found in get_default_regressors()")


def _load_or_build_split():
    """Load raw alpha/gamma waveforms split into the same train/test sets as
    the main analysis. Also returns per-test-event light-output arrays so
    downstream AUCs can be computed per LO bin. Caches to disk after the
    first call; older caches without LO arrays are rebuilt automatically."""
    if os.path.exists(RAW_DATA_CACHE):
        d = np.load(RAW_DATA_CACHE)
        if "alpha_test_lo" in d.files and "gamma_test_lo" in d.files:
            print(f"Loading raw split from {RAW_DATA_CACHE}")
            return (d["alpha_train"], d["gamma_train"], d["alpha_test"],
                    d["gamma_test"], d["alpha_test_lo"], d["gamma_test_lo"])
        print(f"Cache {RAW_DATA_CACHE} predates LO arrays; rebuilding.")

    print("Loading alpha data (Am-241)...")
    alpha_feat, alpha_wf = load_tree_data(ROOT_FILES_DIR + "Am241.root",
                                          array_branch="Samples")
    print("Loading gamma data (Na-22)...")
    gamma_feat, gamma_wf = load_tree_data(ROOT_FILES_DIR + "Na22.root",
                                          array_branch="Samples")

    a_train_idx, g_train_idx = _get_training_indices(alpha_feat, gamma_feat)

    a_all = np.arange(len(alpha_wf))
    g_all = np.arange(len(gamma_wf))
    a_test_idx = np.setdiff1d(a_all, a_train_idx)
    g_test_idx = np.setdiff1d(g_all, g_train_idx)

    # Light-output filter on both training and test sets so the alpha and
    # gamma populations have matched peak-amplitude distributions.
    a_test_lo_mask = (
        (alpha_feat["light_output"].iloc[a_test_idx] >= LO_LOWER["alpha"])
        & (alpha_feat["light_output"].iloc[a_test_idx]
           <= LO_UPPER["alpha"])).values
    g_test_lo_mask = (
        (gamma_feat["light_output"].iloc[g_test_idx] >= LO_LOWER["gamma"])
        & (gamma_feat["light_output"].iloc[g_test_idx]
           <= LO_UPPER["gamma"])).values
    a_train_lo_mask = (
        (alpha_feat["light_output"].iloc[a_train_idx] >= LO_LOWER["alpha"])
        & (alpha_feat["light_output"].iloc[a_train_idx]
           <= LO_UPPER["alpha"])).values
    g_train_lo_mask = (
        (gamma_feat["light_output"].iloc[g_train_idx] >= LO_LOWER["gamma"])
        & (gamma_feat["light_output"].iloc[g_train_idx]
           <= LO_UPPER["gamma"])).values

    a_test = alpha_wf[a_test_idx][a_test_lo_mask]
    g_test = gamma_wf[g_test_idx][g_test_lo_mask]
    a_test_lo = alpha_feat["light_output"].iloc[a_test_idx].values[
        a_test_lo_mask]
    g_test_lo = gamma_feat["light_output"].iloc[g_test_idx].values[
        g_test_lo_mask]

    a_train = alpha_wf[a_train_idx][a_train_lo_mask]
    g_train = gamma_wf[g_train_idx][g_train_lo_mask]

    np.savez(RAW_DATA_CACHE,
             alpha_train=a_train,
             gamma_train=g_train,
             alpha_test=a_test,
             gamma_test=g_test,
             alpha_test_lo=a_test_lo,
             gamma_test_lo=g_test_lo)
    print(f"Raw split cached to {RAW_DATA_CACHE}")
    return a_train, g_train, a_test, g_test, a_test_lo, g_test_lo


def _score_clean_trained(model, alpha_test_noisy, gamma_test_noisy,
                         alpha_test_lo, gamma_test_lo):
    """Score model on test set; return LO-bin-averaged AUC.

    AUC is computed per LO bin with balanced subsampling, then combined
    across bins via inverse-variance weighting so noise-suppression on the
    minority class doesn't dominate the reported value.
    """
    X = np.vstack([
        process_waveforms(alpha_test_noisy),
        process_waveforms(gamma_test_noisy)
    ])
    y_pred = model.predict(X)
    alpha_scores = y_pred[:len(alpha_test_noisy)]
    gamma_scores = y_pred[len(alpha_test_noisy):]
    aucs, errs = per_lo_auc_from_scores(alpha_scores, gamma_scores,
                                        alpha_test_lo, gamma_test_lo)
    return error_weighted_auc(aucs, errs)


def _train_and_score(alpha_train,
                     gamma_train,
                     alpha_test,
                     gamma_test,
                     alpha_test_lo,
                     gamma_test_lo,
                     compute_shap=False):
    """Train one XGBoost per seed in SHAP_SEEDS, score each on the test set
    via LO-bin-averaged AUC, and aggregate across seeds.

    AUC mean = mean over seeds of per-seed LO-averaged AUCs. AUC error
    combines per-seed combined-bin uncertainty (averaged across seeds) and
    between-seed variance in quadrature so it reflects both sources.

    If compute_shap is True, also compute mean |SHAP| per input feature
    per seed and average across seeds, using interventional TreeSHAP on
    the noisy training set.

    Returns:
        compute_shap=False -> (auc_mean, auc_err)
        compute_shap=True  -> (auc_mean, auc_err, mean_shap)
    """
    X_train = np.vstack(
        [process_waveforms(alpha_train),
         process_waveforms(gamma_train)])
    y_train = np.array([0] * len(alpha_train) + [1] * len(gamma_train))

    seed_auc_means = []
    seed_auc_stds = []
    seed_shaps = []
    for seed in SHAP_SEEDS:
        model = _make_xgb(seed=seed)
        model.fit(X_train, y_train)
        auc_mean, auc_std = _score_clean_trained(model, alpha_test,
                                                 gamma_test, alpha_test_lo,
                                                 gamma_test_lo)
        seed_auc_means.append(auc_mean)
        seed_auc_stds.append(auc_std)
        if compute_shap:
            seed_shaps.append(
                compute_mean_abs_shap(model,
                                      X_train,
                                      seed,
                                      explainer="tree",
                                      n_explain=SHAP_N_EXPLAIN))

    auc_mean = float(np.mean(seed_auc_means))
    auc_err = float(
        np.sqrt(
            np.mean(np.square(seed_auc_stds)) +
            np.var(seed_auc_means, ddof=0)))
    if compute_shap:
        mean_shap = np.mean(np.stack(seed_shaps, axis=0), axis=0)
        return auc_mean, auc_err, mean_shap
    return auc_mean, auc_err


def _load_sweep_cache(path, want_shap_key):
    """Load a per-level sweep cache file. Returns
    (auc_mean, auc_std, mean_shap_or_None) or None if the file is missing
    or doesn't contain the expected keys. want_shap_key is "shap" if a
    SHAP array is required, else None."""
    if not os.path.exists(path):
        return None
    data = np.load(path)
    if "auc_mean" not in data or "auc_std" not in data:
        return None
    if want_shap_key is not None and want_shap_key not in data:
        return None
    mean_shap = (np.asarray(data[want_shap_key])
                 if want_shap_key is not None else None)
    return float(data["auc_mean"]), float(data["auc_std"]), mean_shap


def _save_sweep_cache(path, auc_mean, auc_std, mean_shap):
    """Save per-level sweep result. mean_shap may be None."""
    kwargs = dict(auc_mean=np.float64(auc_mean), auc_std=np.float64(auc_std))
    if mean_shap is not None:
        kwargs["shap"] = np.asarray(mean_shap)
    np.savez(path, **kwargs)


def _sweep_test_only(clean_model,
                     alpha_test,
                     gamma_test,
                     alpha_test_lo,
                     gamma_test_lo,
                     levels,
                     noise_fn,
                     label,
                     compute_shap=False):
    """Score the unchanged clean_model on noisy test data at each level.
    If compute_shap=True, also compute mean |SHAP| from clean_model on the
    noisy test set per level so feature usage of the unchanged model can be
    overlaid across noise levels.

    Each (label, level) result is cached to noise_study_cache/ so adding
    or removing levels only triggers work on the new ones. Delete the
    matching .npz files to force recomputation (e.g. after changing
    SHAP_N_EXPLAIN or SHAP_SEEDS).

    Returns (aucs, aucs_err) when compute_shap is False, else
    (aucs, aucs_err, shaps_per_level)."""
    aucs = []
    aucs_err = []
    shaps_per_level = [] if compute_shap else None
    for level in levels:
        cache_path = os.path.join(
            CACHE_DIR, f"test_only_{LO_TAG}_{label}_lvl{level:g}.npz")
        want_shap_key = "shap" if compute_shap else None
        cached = _load_sweep_cache(cache_path, want_shap_key)
        if cached is not None:
            auc_mean, auc_std, mean_shap = cached
            aucs.append(auc_mean)
            aucs_err.append(auc_std)
            if compute_shap:
                shaps_per_level.append(mean_shap)
            print(f"  test-only  {label}={level}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  (cached)")
            continue

        rng = np.random.RandomState(NOISE_SEED)
        a = noise_fn(alpha_test, level, rng)
        g = noise_fn(gamma_test, level, rng)
        auc_mean, auc_std = _score_clean_trained(clean_model, a, g,
                                                 alpha_test_lo,
                                                 gamma_test_lo)
        aucs.append(auc_mean)
        aucs_err.append(auc_std)
        mean_shap = None
        if compute_shap:
            X_test = np.vstack([process_waveforms(a), process_waveforms(g)])
            mean_shap = compute_mean_abs_shap(clean_model,
                                              X_test,
                                              SHAP_SEEDS[0],
                                              explainer="tree",
                                              n_explain=SHAP_N_EXPLAIN)
            shaps_per_level.append(mean_shap)
            print(f"  test-only  {label}={level}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  (SHAP computed)")
        else:
            print(f"  test-only  {label}={level}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}")
        _save_sweep_cache(cache_path, auc_mean, auc_std, mean_shap)
    if compute_shap:
        return aucs, aucs_err, shaps_per_level
    return aucs, aucs_err


def _sweep_matched(alpha_train, gamma_train, alpha_test, gamma_test,
                   alpha_test_lo, gamma_test_lo, levels, noise_fn, label):
    """Train matched-noise XGBoost (3 seeds) at every level and collect
    AUC mean + combined error + mean |SHAP| per level. Returns
    (aucs, aucs_err, shaps_per_level), where shaps_per_level is a list of
    1D numpy arrays of length n_features.

    Each (label, level) result is cached to noise_study_cache/ so adding
    or removing levels only triggers work on the new ones. Delete the
    matching .npz files to force recomputation (e.g. after changing
    SHAP_N_EXPLAIN or SHAP_SEEDS)."""
    aucs = []
    aucs_err = []
    shaps_per_level = []
    for level in levels:
        cache_path = os.path.join(
            CACHE_DIR, f"matched_{LO_TAG}_{label}_lvl{level:g}.npz")
        cached = _load_sweep_cache(cache_path, "shap")
        if cached is not None:
            auc_mean, auc_std, mean_shap = cached
            aucs.append(auc_mean)
            aucs_err.append(auc_std)
            shaps_per_level.append(mean_shap)
            print(f"  matched    {label}={level}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  (cached)")
            continue

        rng = np.random.RandomState(NOISE_SEED)
        a_tr = noise_fn(alpha_train, level, rng)
        g_tr = noise_fn(gamma_train, level, rng)
        a_te = noise_fn(alpha_test, level, rng)
        g_te = noise_fn(gamma_test, level, rng)
        auc_mean, auc_std, mean_shap = _train_and_score(a_tr,
                                                        g_tr,
                                                        a_te,
                                                        g_te,
                                                        alpha_test_lo,
                                                        gamma_test_lo,
                                                        compute_shap=True)
        aucs.append(auc_mean)
        aucs_err.append(auc_std)
        shaps_per_level.append(mean_shap)
        print(f"  matched    {label}={level}  "
              f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  "
              f"(SHAP averaged over {len(SHAP_SEEDS)} seeds)")
        _save_sweep_cache(cache_path, auc_mean, auc_std, mean_shap)
    return aucs, aucs_err, shaps_per_level


def _plot_two_spectra(freqs_a_mhz,
                      mag_a,
                      freqs_g_mhz,
                      mag_g,
                      output_name,
                      name_tag,
                      y_title,
                      log_x=False,
                      f_min_mhz=None):
    """Overlay two max-normalized spectra on log-y vs frequency in MHz.
    Inputs are already in the frequency domain. name_tag goes into the
    TGraph names so the output .root file is searchable per quantity.
    f_min_mhz: if set, x range starts at f_min_mhz and y range is fit to
    only the data in that window (zoomed view, still max-normalized over
    the full spectrum)."""
    freqs_a = np.asarray(freqs_a_mhz, dtype=np.float64)
    freqs_g = np.asarray(freqs_g_mhz, dtype=np.float64)
    a_norm = (np.asarray(mag_a, dtype=np.float64) /
              float(np.asarray(mag_a).max()))
    g_norm = (np.asarray(mag_g, dtype=np.float64) /
              float(np.asarray(mag_g).max()))

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas.SetLogy(True)
    if log_x:
        canvas.SetLogx(True)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)

    g_alpha = ROOT.TGraph(len(freqs_a), freqs_a.astype(np.float64), a_norm)
    g_alpha.SetName(f"alpha_{name_tag}")
    g_alpha.SetTitle("")
    g_alpha.SetLineColor(ROOT.kRed + 2)
    g_alpha.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_alpha.GetXaxis().SetTitle("Frequency [MHz]")
    g_alpha.GetYaxis().SetTitle(y_title)
    g_alpha.GetYaxis().SetTitleOffset(1.2)
    if f_min_mhz is not None:
        x_lo = float(f_min_mhz)
        mask_a = freqs_a >= f_min_mhz
        mask_g = freqs_g >= f_min_mhz
        vis_a = a_norm[mask_a]
        vis_g = g_norm[mask_g]
        vis_pos = np.concatenate([vis_a[vis_a > 0], vis_g[vis_g > 0]])
        y_floor = float(vis_pos.min()) if len(vis_pos) else 1e-12
        y_top = float(max(vis_a.max(), vis_g.max())) * 1.5
    else:
        x_lo = float(freqs_a[1]) * 0.5 if log_x else 0.0
        y_floor = float(min(a_norm[a_norm > 0].min(),
                            g_norm[g_norm > 0].min()))
        y_top = 1.5
    g_alpha.GetXaxis().SetLimits(x_lo, float(freqs_a[-1]))
    g_alpha.GetYaxis().SetRangeUser(y_floor * 0.5, y_top)
    g_alpha.Draw("AL")

    g_gamma = ROOT.TGraph(len(freqs_g), freqs_g.astype(np.float64), g_norm)
    g_gamma.SetName(f"gamma_{name_tag}")
    g_gamma.SetLineColor(ROOT.kBlue + 2)
    g_gamma.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_gamma.Draw("L SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.55, 0.8, 0.55, 0.88)
    leg.AddEntry(g_alpha, "Am-241 (#alpha)", "l")
    leg.AddEntry(g_gamma, "Na-22 (#gamma)", "l")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLOG)
    canvas.Close()
    print(f"Saved {output_name}")


def _signal_fft(signal_1d, power):
    """rfft of a 1D time-domain signal. Returns (freqs_mhz, values) where
    values is |FFT| if power=False, |FFT|^2 if power=True."""
    n = len(signal_1d)
    mag = np.abs(np.fft.rfft(signal_1d))
    if power:
        mag = mag**2
    freqs_mhz = np.fft.rfftfreq(n, d=1.0 / SAMPLE_RATE_HZ) * 1e-6
    return freqs_mhz, mag


def _plot_avg_waveform_fft(alpha_test,
                           gamma_test,
                           output_name,
                           power=False,
                           f_min_mhz=None):
    """FFT of the class-average alpha and gamma waveforms."""
    freqs_a, mag_a = _signal_fft(np.mean(alpha_test, axis=0), power)
    freqs_g, mag_g = _signal_fft(np.mean(gamma_test, axis=0), power)
    y_title = ("Normalized |FFT|^{2} [a.u.]"
               if power else "Normalized |FFT| [a.u.]")
    _plot_two_spectra(freqs_a,
                      mag_a,
                      freqs_g,
                      mag_g,
                      output_name,
                      f"avg_waveform_{'power' if power else 'mag'}",
                      y_title,
                      log_x=power and f_min_mhz is None,
                      f_min_mhz=f_min_mhz)


def _residual_psd(events):
    """Ensemble-averaged single-event noise PSD: mean over events of
    |FFT(event - <event>)|^2. Returns (freqs_mhz, psd)."""
    residuals = events - np.mean(events, axis=0, keepdims=True)
    spec = np.abs(np.fft.rfft(residuals, axis=-1))**2
    psd = np.mean(spec, axis=0)
    freqs_mhz = np.fft.rfftfreq(events.shape[-1],
                                d=1.0 / SAMPLE_RATE_HZ) * 1e-6
    return freqs_mhz, psd


def _plot_residual_psd(alpha_test,
                       gamma_test,
                       output_name,
                       log_x=True,
                       f_min_mhz=None):
    """Ensemble-averaged per-event noise PSD. Reveals periodic noise
    sources (mains, switching supplies), 1/f drift, and ringing — none
    of which the across-event variance FFT can separate from the
    pulse-envelope structure."""
    freqs_a, psd_a = _residual_psd(alpha_test)
    freqs_g, psd_g = _residual_psd(gamma_test)
    _plot_two_spectra(freqs_a,
                      psd_a,
                      freqs_g,
                      psd_g,
                      output_name,
                      "residual_psd",
                      "Normalized noise PSD [a.u.]",
                      log_x=log_x and f_min_mhz is None,
                      f_min_mhz=f_min_mhz)


def _plot_sample_value_distributions(alpha_test,
                                     gamma_test,
                                     output_name,
                                     times_ns=(50, 75, 100, 150)):
    """Per-class histograms of the per-event max-normalized sample value
    at fixed time indices post-trigger. The classical-PSD / matched-
    filter family is sensitive to per-class differences in the *mean*
    waveform shape; XGBoost can additionally exploit per-class
    differences in the *width* of the sample distribution at fixed time
    (the signature of class-dependent pulse-shape jitter). If the
    histograms have identical widths at every panel, the non-linear ML
    lift isn't coming from jitter exploitation; if widths visibly
    differ, that's the smoking gun."""
    a_norm = process_waveforms(alpha_test)
    g_norm = process_waveforms(gamma_test)

    canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "", 1200,
                          900)
    canvas.Divide(2, 2)
    keep = []
    for k, t_ns in enumerate(times_ns):
        canvas.cd(k + 1)
        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetRightMargin(0.04)
        ROOT.gPad.SetTopMargin(0.1)
        ROOT.gPad.SetBottomMargin(0.14)
        idx = PRE_SAMPLES + t_ns // 2
        if idx >= a_norm.shape[1] or idx >= g_norm.shape[1]:
            continue
        a_vals = a_norm[:, idx].astype(np.float64)
        g_vals = g_norm[:, idx].astype(np.float64)
        all_vals = np.concatenate([a_vals, g_vals])
        vmin = float(np.percentile(all_vals, 0.5))
        vmax = float(np.percentile(all_vals, 99.5))
        if vmax <= vmin:
            continue
        nbins = 120
        ha = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", nbins,
                       vmin, vmax)
        hg = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", nbins,
                       vmin, vmax)
        ha.SetDirectory(0)
        hg.SetDirectory(0)
        ha.FillN(len(a_vals), a_vals, np.ones(len(a_vals), dtype=np.float64))
        hg.FillN(len(g_vals), g_vals, np.ones(len(g_vals), dtype=np.float64))
        if ha.Integral() > 0:
            ha.Scale(1.0 / ha.Integral())
        if hg.Integral() > 0:
            hg.Scale(1.0 / hg.Integral())
        ha.SetLineColor(ROOT.kRed + 2)
        ha.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        ha.SetFillColorAlpha(ROOT.kRed + 2, 0.3)
        hg.SetLineColor(ROOT.kBlue + 2)
        hg.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        hg.SetFillColorAlpha(ROOT.kBlue + 2, 0.3)
        ha.SetTitle(f"t = {t_ns} ns post-trigger")
        ha.GetXaxis().SetTitle("Normalized sample value [a.u.]")
        ha.GetYaxis().SetTitle("Density")
        ha.GetYaxis().SetTitleOffset(1.3)
        y_max = max(ha.GetMaximum(), hg.GetMaximum()) * 1.25
        ha.SetMaximum(y_max)
        ha.Draw("HIST")
        hg.Draw("HIST SAME")
        a_mu = float(np.mean(a_vals))
        a_sd = float(np.std(a_vals))
        g_mu = float(np.mean(g_vals))
        g_sd = float(np.std(g_vals))
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.038)
        latex.SetTextColor(ROOT.kRed + 2)
        latex.DrawLatex(
            0.50, 0.83,
            f"#mu_{{#alpha}}={a_mu:.3f}, #sigma_{{#alpha}}={a_sd:.3f}")
        latex2 = ROOT.TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.038)
        latex2.SetTextColor(ROOT.kBlue + 2)
        latex2.DrawLatex(
            0.50, 0.78,
            f"#mu_{{#gamma}}={g_mu:.3f}, #sigma_{{#gamma}}={g_sd:.3f}")
        ratio = g_sd / a_sd if a_sd > 0 else float("nan")
        latex3 = ROOT.TLatex()
        latex3.SetNDC()
        latex3.SetTextSize(0.038)
        latex3.SetTextColor(ROOT.kBlack)
        latex3.DrawLatex(0.50, 0.73,
                         f"#sigma_{{#gamma}}/#sigma_{{#alpha}}={ratio:.3f}")
        keep.extend([ha, hg, latex, latex2, latex3])

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _apply_butter_lpf(wf, cutoff_mhz):
    """Zero-phase 4th-order Butterworth low-pass at cutoff_mhz. Operates
    on the last axis so 1D or 2D (n_events, n_samples) arrays both work."""
    wn = cutoff_mhz * 1e6 / (SAMPLE_RATE_HZ / 2.0)
    sos = butter(LPF_ORDER, wn, btype="lowpass", output="sos")
    return sosfiltfilt(sos, wf, axis=-1)


def _apply_gaussian(wf, sigma_samples):
    """Gaussian smoothing with kernel sigma in samples. Reflect-mode at
    edges to avoid baseline distortion."""
    return gaussian_filter1d(wf, sigma=sigma_samples, axis=-1, mode="reflect")


def _apply_notches(wf, freqs_mhz, Q=NOTCH_Q):
    """Cascade narrow IIR notches at each freq in freqs_mhz, applied
    zero-phase via sosfiltfilt. Q controls width: bandwidth = f0 / Q, so
    Q = 30 at 100 MHz -> ~3.3 MHz notch. Operates on the last axis."""
    out = wf
    for f0 in freqs_mhz:
        b, a = iirnotch(f0 * 1e6, Q, fs=SAMPLE_RATE_HZ)
        sos = tf2sos(b, a)
        out = sosfiltfilt(sos, out, axis=-1)
    return out


def _run_notch_experiment(alpha_train, gamma_train, alpha_test, gamma_test,
                          alpha_test_lo, gamma_test_lo):
    """Train+test XGBoost (3 seeds) with each of: no notch, each single
    notch in NOTCH_FREQS_MHZ, and all notches combined. Per-config result
    cached to noise_study_cache/. Returns list of (label, auc, auc_err)."""
    configs = [("none", [])]
    for f in NOTCH_FREQS_MHZ:
        configs.append((f"f{f:g}", [f]))
    configs.append(("all", list(NOTCH_FREQS_MHZ)))

    results = []
    for label, freqs in configs:
        cache_path = os.path.join(CACHE_DIR, f"notch_{LO_TAG}_{label}.npz")
        cached = _load_sweep_cache(cache_path, None)
        if cached is not None:
            auc_mean, auc_std, _ = cached
            print(f"  notch    {label}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  (cached)")
        else:
            a_tr = (_apply_notches(alpha_train, freqs)
                    if freqs else alpha_train)
            g_tr = (_apply_notches(gamma_train, freqs)
                    if freqs else gamma_train)
            a_te = _apply_notches(alpha_test, freqs) if freqs else alpha_test
            g_te = _apply_notches(gamma_test, freqs) if freqs else gamma_test
            auc_mean, auc_std = _train_and_score(a_tr,
                                                 g_tr,
                                                 a_te,
                                                 g_te,
                                                 alpha_test_lo,
                                                 gamma_test_lo,
                                                 compute_shap=False)
            print(f"  notch    {label}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}")
            _save_sweep_cache(cache_path, auc_mean, auc_std, None)
        results.append((label, auc_mean, auc_std))
    return results


def _plot_notch_results(results, output_name):
    """Bar plot of AUC per notch configuration. results is the list
    returned by _run_notch_experiment."""
    n = len(results)
    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.18)

    frame = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", n, 0.5,
                      n + 0.5)
    frame.SetStats(False)
    aucs = [r[1] for r in results]
    errs = [r[2] for r in results]
    y_lo = min(a - e for a, e in zip(aucs, errs))
    y_hi = max(a + e for a, e in zip(aucs, errs))
    pad = max(0.005, 0.2 * (y_hi - y_lo))
    frame.GetYaxis().SetRangeUser(y_lo - pad, y_hi + pad)
    frame.GetYaxis().SetTitle("ROC AUC")
    frame.GetXaxis().SetLabelSize(0.04)
    for i, (label, _, _) in enumerate(results):
        display = ("baseline" if label == "none" else
                   "all notches" if label == "all" else f"{label[1:]} MHz")
        frame.GetXaxis().SetBinLabel(i + 1, display)
    frame.Draw()

    x = np.arange(1, n + 1, dtype=np.float64)
    ex = np.zeros(n, dtype=np.float64)
    g = ROOT.TGraphErrors(n, x, np.array(aucs, dtype=np.float64), ex,
                          np.array(errs, dtype=np.float64))
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.3)
    g.SetMarkerColor(ROOT.kBlue + 2)
    g.SetLineColor(ROOT.kBlue + 2)
    g.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g.Draw("P SAME")

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _sweep_filter(alpha_train, gamma_train, alpha_test, gamma_test,
                  alpha_test_lo, gamma_test_lo, levels, filter_fn,
                  filter_tag):
    """Apply filter_fn(wf, level) to alpha+gamma train and test waveforms
    independently, train fresh XGBoost (3 seeds), bootstrap-score on the
    filtered test set. Returns (aucs, aucs_err) per level.

    Each (filter_tag, level) result is cached to noise_study_cache/ so
    adding or removing levels only triggers work on the new ones."""
    aucs = []
    aucs_err = []
    for level in levels:
        cache_path = os.path.join(
            CACHE_DIR, f"filter_{LO_TAG}_{filter_tag}_lvl{level:g}.npz")
        cached = _load_sweep_cache(cache_path, None)
        if cached is not None:
            auc_mean, auc_std, _ = cached
            aucs.append(auc_mean)
            aucs_err.append(auc_std)
            print(f"  filter   {filter_tag}={level:g}  "
                  f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}  (cached)")
            continue

        a_tr = filter_fn(alpha_train, level)
        g_tr = filter_fn(gamma_train, level)
        a_te = filter_fn(alpha_test, level)
        g_te = filter_fn(gamma_test, level)
        auc_mean, auc_std = _train_and_score(a_tr,
                                             g_tr,
                                             a_te,
                                             g_te,
                                             alpha_test_lo,
                                             gamma_test_lo,
                                             compute_shap=False)
        aucs.append(auc_mean)
        aucs_err.append(auc_std)
        print(f"  filter   {filter_tag}={level:g}  "
              f"AUC = {auc_mean:.4f} +/- {auc_std:.4f}")
        _save_sweep_cache(cache_path, auc_mean, auc_std, None)
    return aucs, aucs_err


def _plot_filter_auc(levels, aucs, errs, baseline_auc, baseline_auc_err,
                     x_title, curve_label, curve_color, marker_style,
                     output_name):
    """AUC vs filter knob with the clean baseline AUC drawn as a +/-1 sigma
    band. levels are plotted in their native units (no conversion)."""
    canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "", 1200,
                          600)
    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.72, 1.0)
    pad_plot.SetRightMargin(0.02)
    pad_plot.Draw()
    pad_leg = ROOT.TPad("pad_leg", "", 0.72, 0.0, 1.0, 1.0)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_leg.Draw()
    pad_plot.cd()

    x_arr = np.array(levels, dtype=np.float64)
    ex_arr = np.zeros(len(x_arr), dtype=np.float64)
    g = ROOT.TGraphErrors(len(x_arr), x_arr, np.array(aucs, dtype=np.float64),
                          ex_arr, np.array(errs, dtype=np.float64))
    g.SetLineColor(curve_color)
    g.SetMarkerColor(curve_color)
    g.SetMarkerStyle(marker_style)
    g.SetMarkerSize(1.2)
    g.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

    all_lo = ([a - e for a, e in zip(aucs, errs)] +
              [baseline_auc - baseline_auc_err])
    all_hi = ([a + e for a, e in zip(aucs, errs)] +
              [baseline_auc + baseline_auc_err])
    y_min = min(all_lo)
    y_max = max(all_hi)
    pad = max(0.005, 0.2 * (y_max - y_min))

    x_range = float(x_arr.max() - x_arr.min())
    x_lo = float(x_arr.min()) - 0.05 * x_range
    x_hi = float(x_arr.max()) + 0.05 * x_range

    g.SetTitle("")
    g.GetXaxis().SetTitle(x_title)
    g.GetYaxis().SetTitle("ROC AUC")
    g.GetYaxis().SetRangeUser(y_min - pad, y_max + pad)
    g.GetXaxis().SetLimits(x_lo, x_hi)
    g.Draw("APL")

    band_x = np.array([x_lo, x_hi], dtype=np.float64)
    band_ex = np.zeros(2, dtype=np.float64)
    base_band = ROOT.TGraphErrors(
        2, band_x, np.array([baseline_auc, baseline_auc], dtype=np.float64),
        band_ex,
        np.array([baseline_auc_err, baseline_auc_err], dtype=np.float64))
    base_band.SetFillColorAlpha(ROOT.kGray + 2, 0.3)
    base_band.SetLineColor(ROOT.kGray + 2)
    base_band.SetLineStyle(2)
    base_band.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    base_band.Draw("L3 SAME")

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.0, 0.95, 0.25, 0.85)
    leg.SetMargin(0.15)
    leg.AddEntry(g, curve_label, "lpe")
    leg.AddEntry(base_band, "#splitline{Unmodified}{baseline}", "lf")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _plot_waveform_panel(waveforms,
                         level_label,
                         output_name,
                         output_subdir,
                         n_samples=5):
    """Overlay n_samples normalized waveforms on one canvas."""
    n_points = waveforms.shape[1]
    x_values = np.arange(n_points, dtype=np.float64) * 2

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    colors = [
        ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kOrange + 2,
        ROOT.kMagenta + 2
    ]

    graphs = []
    for i in range(min(n_samples, len(waveforms))):
        graph = ROOT.TGraph(n_points, x_values,
                            waveforms[i].astype(np.float64))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        if i == 0:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Time [ns]")
            graph.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
            graph.GetXaxis().SetRangeUser(0, x_values[-1])
            graph.GetYaxis().SetRangeUser(-0.2, 1.2)
            graph.Draw("AL")
        else:
            graph.Draw("L SAME")
        graphs.append(graph)

    _ = ROOT.PlottingUtils.AddText(level_label, 0.78, 0.85)

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, output_subdir,
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()


def _plot_noise_sample_waveforms(alpha_test,
                                 gamma_test,
                                 noise_levels,
                                 noise_fn,
                                 level_name,
                                 level_unit,
                                 output_subdir,
                                 n_samples=5,
                                 pick_seed=7):
    """Plot 5 normalized waveforms each for alpha and gamma at every noise
    level (including 0 = unmodified baseline). The same events and the same
    noise-realization seed are used at every level, so visual comparison
    isolates the effect of changing the noise magnitude.
    """
    pick_rng = np.random.RandomState(pick_seed)
    a_idx = pick_rng.choice(len(alpha_test), size=n_samples, replace=False)
    g_idx = pick_rng.choice(len(gamma_test), size=n_samples, replace=False)
    a_sample = alpha_test[a_idx]
    g_sample = gamma_test[g_idx]

    for level in [0.0] + list(noise_levels):
        rng = np.random.RandomState(NOISE_SEED)
        a_noisy = noise_fn(a_sample, level, rng)
        g_noisy = noise_fn(g_sample, level, rng)
        a_norm = process_waveforms(a_noisy)
        g_norm = process_waveforms(g_noisy)

        level_tag = f"{level:g}".replace(".", "p")
        unit_str = f" {level_unit}" if level_unit else ""
        alpha_label = f"Am-241 (#alpha), {level_name}={level:g}{unit_str}"
        gamma_label = f"Na-22 (#gamma), {level_name}={level:g}{unit_str}"

        _plot_waveform_panel(a_norm, alpha_label,
                             f"sample_alpha_{level_name}_{level_tag}",
                             output_subdir, n_samples)
        _plot_waveform_panel(g_norm, gamma_label,
                             f"sample_gamma_{level_name}_{level_tag}",
                             output_subdir, n_samples)

    print(f"Saved {2 * (1 + len(noise_levels))} waveform plots to "
          f"plots/{output_subdir}/")


def _compute_long_integral(waveforms):
    """Reproduce the long_integral feature from the C++ pipeline.

    Matches WaveformProcessingUtils::ExtractFeatures: sum of samples from
    INTEGRATION_START to INTEGRATION_START+LONG_GATE, clamped to the array
    end. Waveforms are already baseline-subtracted and polarity-inverted
    upstream, so this is a direct port.
    """
    n = waveforms.shape[1]
    long_end = min(INTEGRATION_START + LONG_GATE, n)
    return np.sum(waveforms[:, INTEGRATION_START:long_end], axis=1)


def _make_hist(values, n_bins, lo, hi):
    """Create and fill a TH1F using per-value Fill, matching the pattern
    used elsewhere in the codebase (e.g. _plot_classified_spectra in
    proof_of_concept.py)."""
    h = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", n_bins, lo, hi)
    for v in values:
        h.Fill(v)
    return h


def _plot_lo_overlay(am_unmod,
                     am_noisy,
                     na_unmod,
                     na_noisy,
                     level_label,
                     output_name,
                     output_subdir,
                     lo_range,
                     n_bins=200):
    """Four overlaid LO histograms on one canvas:
        Am-241 unmodified (dashed red) vs noisy (solid red),
        Na-22  unmodified (dashed blue) vs noisy (solid blue).
    """
    canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "", 1200,
                          600)
    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.72, 1.0)
    pad_plot.SetLogy()
    pad_plot.SetLeftMargin(0.12)
    pad_plot.SetRightMargin(0.05)
    pad_plot.Draw()
    pad_leg = ROOT.TPad("pad_leg", "", 0.72, 0.0, 1.0, 1.0)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_leg.Draw()
    pad_plot.cd()

    lo_min, lo_max = lo_range
    h_am_u = _make_hist(am_unmod, n_bins, lo_min, lo_max)
    h_am_n = _make_hist(am_noisy, n_bins, lo_min, lo_max)
    h_na_u = _make_hist(na_unmod, n_bins, lo_min, lo_max)
    h_na_n = _make_hist(na_noisy, n_bins, lo_min, lo_max)

    ROOT.PlottingUtils.ConfigureHistogram(h_am_u, ROOT.kRed + 2)
    ROOT.PlottingUtils.ConfigureHistogram(h_am_n, ROOT.kRed + 2)
    ROOT.PlottingUtils.ConfigureHistogram(h_na_u, ROOT.kBlue + 2)
    ROOT.PlottingUtils.ConfigureHistogram(h_na_n, ROOT.kBlue + 2)
    h_am_u.SetLineStyle(2)
    h_na_u.SetLineStyle(2)

    h_am_u.GetXaxis().SetTitle("Long integral [ADC]")
    h_am_u.GetYaxis().SetTitle("Counts")
    h_am_u.GetYaxis().SetTitleOffset(1)
    h_am_u.SetTitle("")

    max_val = max(h_am_u.GetMaximum(), h_am_n.GetMaximum(),
                  h_na_u.GetMaximum(), h_na_n.GetMaximum())
    h_am_u.SetMaximum(max_val * 1.4)

    h_am_u.Draw("HIST")
    h_am_n.Draw("HIST SAME")
    h_na_u.Draw("HIST SAME")
    h_na_n.Draw("HIST SAME")

    _ = ROOT.PlottingUtils.AddText(level_label, 0.92, 0.85)

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.0, 0.95, 0.55, 0.85)
    leg.SetMargin(0.15)
    leg.AddEntry(h_am_u, "Am-241 unmodified", "l")
    leg.AddEntry(h_am_n, "Am-241 noisy", "l")
    leg.AddEntry(h_na_u, "Na-22 unmodified", "l")
    leg.AddEntry(h_na_n, "Na-22 noisy", "l")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, output_subdir,
                                  ROOT.PlotSaveOptions.kLOG)
    canvas.Close()
    h_am_u.Delete()
    h_am_n.Delete()
    h_na_u.Delete()
    h_na_n.Delete()


def _plot_noise_light_output(alpha_test, gamma_test, noise_levels, noise_fn,
                             level_name, level_unit, output_subdir):
    """For each noise level (plus the unmodified baseline), plot LO spectra
    overlaying the unmodified data against the same data with noise added.

    Noise is applied directly to the full waveform array; the long_integral
    is then recomputed using the same algorithm as the C++ pipeline.
    """
    am_unmod = _compute_long_integral(alpha_test)
    na_unmod = _compute_long_integral(gamma_test)

    all_lo = np.concatenate([am_unmod, na_unmod])
    lo_lo, lo_hi = np.percentile(all_lo, [0.1, 99.9])
    span = lo_hi - lo_lo
    lo_range = (lo_lo - 0.3 * span, lo_hi + 0.5 * span)

    for level in [0.0] + list(noise_levels):
        rng = np.random.RandomState(NOISE_SEED)
        am_noisy = _compute_long_integral(noise_fn(alpha_test, level, rng))
        na_noisy = _compute_long_integral(noise_fn(gamma_test, level, rng))

        level_tag = f"{level:g}".replace(".", "p")
        unit_str = f" {level_unit}" if level_unit else ""
        level_label = f"{level_name}={level:g}{unit_str}"
        output_name = f"lo_{level_name}_{level_tag}"
        _plot_lo_overlay(am_unmod, am_noisy, na_unmod, na_noisy, level_label,
                         output_name, output_subdir, lo_range)

    print(f"Saved {1 + len(noise_levels)} LO spectra to "
          f"plots/{output_subdir}/")


def _plot_variance_vs_mean(alpha_test,
                           gamma_test,
                           output_name,
                           min_amplitude_frac=0.03):
    """Per-sample sigma^2(t) vs mu(t) on log-log axes -- bare-data plot
    with no overlays. TGraphs are written to a ROOT file in CACHE_DIR
    with linear (mu, sigma^2) coordinates so the file can be opened in
    TBrowser and fitted interactively.
    """
    mean_a = np.mean(alpha_test, axis=0)
    var_a = np.var(alpha_test, axis=0)
    mean_g = np.mean(gamma_test, axis=0)
    var_g = np.var(gamma_test, axis=0)

    thresh_a = min_amplitude_frac * float(mean_a.max())
    thresh_g = min_amplitude_frac * float(mean_g.max())
    mask_a = (mean_a > thresh_a) & (var_a > 0)
    mask_g = (mean_g > thresh_g) & (var_g > 0)
    mu_a = mean_a[mask_a].astype(np.float64)
    v_a = var_a[mask_a].astype(np.float64)
    mu_g = mean_g[mask_g].astype(np.float64)
    v_g = var_g[mask_g].astype(np.float64)

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas.SetLogx(True)
    canvas.SetLogy(True)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)

    g_alpha = ROOT.TGraph(len(mu_a), mu_a, v_a)
    g_alpha.SetName("alpha_variance_vs_mean")
    g_alpha.SetTitle("")
    g_alpha.SetMarkerColor(ROOT.kRed + 2)
    g_alpha.SetMarkerStyle(20)
    g_alpha.SetMarkerSize(0.9)
    g_alpha.GetXaxis().SetTitle("Per-sample mean amplitude #mu(t) [ADC]")
    g_alpha.GetYaxis().SetTitle("Per-sample variance #sigma^{2}(t) [ADC^{2}]")
    g_alpha.GetYaxis().SetTitleOffset(1.2)

    all_mu = np.concatenate([mu_a, mu_g])
    all_v = np.concatenate([v_a, v_g])
    mu_lo = float(all_mu.min()) * 0.5
    mu_hi = float(all_mu.max()) * 2.0
    v_lo = float(all_v.min()) * 0.5
    v_hi = float(all_v.max()) * 2.0
    g_alpha.GetXaxis().SetLimits(mu_lo, mu_hi)
    g_alpha.GetYaxis().SetRangeUser(v_lo, v_hi)
    g_alpha.Draw("AP")

    g_gamma = ROOT.TGraph(len(mu_g), mu_g, v_g)
    g_gamma.SetName("gamma_variance_vs_mean")
    g_gamma.SetMarkerColor(ROOT.kBlue + 2)
    g_gamma.SetMarkerStyle(21)
    g_gamma.SetMarkerSize(0.9)
    g_gamma.Draw("P SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.17, 0.35, 0.65, 0.88)
    leg.AddEntry(g_alpha, "Am-241 (#alpha)", "p")
    leg.AddEntry(g_gamma, "Na-22 (#gamma)", "p")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLOG)
    canvas.Close()
    print(f"Saved {output_name}")

    root_path = os.path.join(CACHE_DIR, f"{output_name}.root")
    out_file = ROOT.TFile(root_path, "RECREATE")
    g_alpha.Write()
    g_gamma.Write()
    out_file.Close()
    print(f"Saved {root_path}")


def _load_xgb_shap_mean():
    """Average XGBoost SHAP importance across seeds from parameter_study.py.
    Returns None if no cached SHAP files are found."""
    arrs = []
    for seed in SHAP_SEEDS:
        path = os.path.join(SWEEP_CACHE_DIR,
                            f"xgb_tree_shap_values_seed{seed}.npy")
        if os.path.exists(path):
            arrs.append(np.load(path))
    if not arrs:
        return None
    return np.mean(np.stack(arrs, axis=0), axis=0)


def _plot_variance_vs_dmu_dt(alpha_test,
                             gamma_test,
                             output_name,
                             min_amplitude_frac=0.03,
                             n_top=7):
    """Per-sample sigma^2(t) vs (dmu/dt)^2(t) on linear axes (matching the
    quadratic-residual plot). Per-event timing jitter delta_t contributes
    Var(delta_t) * (dmu/dt)^2 to sigma^2, so jitter-dominated samples
    should fall on a line through origin with slope Var(delta_t). No
    subtraction or fit overlay -- pure data.
    """
    mean_a = np.mean(alpha_test, axis=0)
    var_a = np.var(alpha_test, axis=0)
    mean_g = np.mean(gamma_test, axis=0)
    var_g = np.var(gamma_test, axis=0)

    # dt = 2 ns/sample. np.gradient: centered differences interior,
    # one-sided at the endpoints.
    dmu_a = np.gradient(mean_a) / 2.0
    dmu_g = np.gradient(mean_g) / 2.0

    thresh_a = min_amplitude_frac * float(mean_a.max())
    thresh_g = min_amplitude_frac * float(mean_g.max())
    mask_a = (mean_a > thresh_a) & (var_a > 0)
    mask_g = (mean_g > thresh_g) & (var_g > 0)

    dx_a = (dmu_a[mask_a]**2).astype(np.float64)
    v_a = var_a[mask_a].astype(np.float64)
    dx_g = (dmu_g[mask_g]**2).astype(np.float64)
    v_g = var_g[mask_g].astype(np.float64)

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)

    g_alpha = ROOT.TGraph(len(dx_a), dx_a, v_a)
    g_alpha.SetName("alpha_variance_vs_dmu_dt_sq")
    g_alpha.SetTitle("")
    g_alpha.SetMarkerColor(ROOT.kRed + 2)
    g_alpha.SetMarkerStyle(20)
    g_alpha.SetMarkerSize(0.9)
    g_alpha.GetXaxis().SetTitle("(d#mu/dt)^{2}(t) [(ADC/ns)^{2}]")
    g_alpha.GetYaxis().SetTitle("Per-sample variance #sigma^{2}(t) [ADC^{2}]")
    g_alpha.GetYaxis().SetTitleOffset(1.2)
    g_alpha.GetXaxis().SetNdivisions(505)

    all_dx = np.concatenate([dx_a, dx_g])
    all_v = np.concatenate([v_a, v_g])
    x_hi = float(all_dx.max()) * 1.05
    v_hi = float(all_v.max()) * 1.05
    g_alpha.GetXaxis().SetLimits(0.0, x_hi)
    g_alpha.GetYaxis().SetRangeUser(0.0, v_hi)
    g_alpha.Draw("AP")

    g_gamma = ROOT.TGraph(len(dx_g), dx_g, v_g)
    g_gamma.SetName("gamma_variance_vs_dmu_dt_sq")
    g_gamma.SetMarkerColor(ROOT.kBlue + 2)
    g_gamma.SetMarkerStyle(21)
    g_gamma.SetMarkerSize(0.9)
    g_gamma.Draw("P SAME")

    # Annotate the n_top alpha samples with the largest (dmu/dt)^2.
    # order_in_dx is descending so labels are listed from largest down.
    sample_idx_a = np.where(mask_a)[0]
    order_in_dx = np.argsort(dx_a)[::-1][:n_top]
    top_orig_idx = sample_idx_a[order_in_dx]

    label_objs = []
    for j, orig_idx in zip(order_in_dx, top_orig_idx):
        lbl = ROOT.TLatex(float(dx_a[j]), float(v_a[j]), f" {int(orig_idx)}")
        lbl.SetTextSize(0.022)
        lbl.SetTextColor(ROOT.kRed + 2)
        lbl.Draw()
        label_objs.append(lbl)

    leg = ROOT.PlottingUtils.AddLegend(0.17, 0.35, 0.5, 0.88)
    leg.AddEntry(g_alpha, "Am-241 (#alpha)", "p")
    leg.AddEntry(g_gamma, "Na-22 (#gamma)", "p")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")

    # Print SHAP values (normalized so max over all features = 1) at the
    # top-n_top (dmu/dt)^2 alpha samples. SHAP_norm = 1 marks the single
    # most important feature; values near 1 flag samples nearly as important.
    shap_mean = _load_xgb_shap_mean()
    if shap_mean is None:
        print(
            f"  SHAP cache not found in {SWEEP_CACHE_DIR}/, skipping printout")
    else:
        shap_max = float(shap_mean.max())
        shap_norm = shap_mean / shap_max if shap_max > 0 else shap_mean
        print(f"  Top {n_top} (d#mu/dt)^2 alpha samples and their normalized "
              f"SHAP (max-normalized, SHAP length = {len(shap_mean)}):")
        print(f"    sample_idx   t_ns   (dmu/dt)^2 [(ADC/ns)^2]   SHAP_norm")
        for orig_idx in top_orig_idx:
            shap_val = (float(shap_norm[orig_idx])
                        if orig_idx < len(shap_norm) else float("nan"))
            dmu_sq = float(dmu_a[orig_idx])**2
            print(f"    {int(orig_idx):10d}   {int(orig_idx) * 2:4d}   "
                  f"{dmu_sq:22.4f}   {shap_val:.4f}")

    root_path = os.path.join(CACHE_DIR, f"{output_name}.root")
    out_file = ROOT.TFile(root_path, "RECREATE")
    g_alpha.Write()
    g_gamma.Write()
    out_file.Close()
    print(f"Saved {root_path}")


def _plot_avg_waveform_outliers(test_data,
                                output_name,
                                p0,
                                p1,
                                p2,
                                class_tag,
                                marker_color,
                                marker_style=20,
                                min_amplitude_frac=0.03,
                                residual_threshold=500.0):
    """Plot the class-average waveform (max-normalized) with max-normalized
    XGBoost SHAP overlaid. Samples whose quadratic-fit residual exceeds
    residual_threshold are highlighted as colored markers on the waveform,
    matching the labels on {class_tag}_quadratic_residual.pdf.

    Layout matches the feature-importance plots in parameter_study.py:
    gray average waveform, colored overlay, linear axes, time on x (ns).
    class_tag is the legend root symbol (e.g., "alpha" -> "#alpha").
    """
    mean_w = np.mean(test_data, axis=0)
    var_w = np.var(test_data, axis=0)

    thresh = min_amplitude_frac * float(mean_w.max())
    mask = (mean_w > thresh) & (var_w > 0)
    sample_indices_kept = np.where(mask)[0]
    mu = mean_w[mask].astype(np.float64)
    v = var_w[mask].astype(np.float64)
    fit_v = p0 + p1 * mu + p2 * mu**2
    residual = v - fit_v
    above = residual > residual_threshold
    outlier_orig_idx = sample_indices_kept[above]

    wf_max = float(np.max(mean_w))
    avg_wf_norm = mean_w / wf_max if wf_max > 0 else mean_w

    n_points = len(mean_w)
    x_values = np.arange(n_points, dtype=np.float64) * 2

    shap_mean = _load_xgb_shap_mean()

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas(False)
    if shap_mean is None:
        pad_top = None
        pad_bot = None
    else:
        pad_top = ROOT.TPad("pad_top", "", 0.0, 0.3, 1.0, 1.0)
        pad_top.SetBottomMargin(0.04)
        pad_top.SetTopMargin(0.12)
        pad_top.Draw()
        pad_bot = ROOT.TPad("pad_bot", "", 0.0, 0.0, 1.0, 0.3)
        pad_bot.SetTopMargin(0.04)
        pad_bot.SetBottomMargin(0.35)
        pad_bot.Draw()
        pad_top.cd()

    g_wf = ROOT.TGraph(n_points, x_values, avg_wf_norm.astype(np.float64))
    g_wf.SetLineColor(ROOT.kGray + 2)
    g_wf.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_wf.SetTitle("")
    if pad_top is None:
        g_wf.GetXaxis().SetTitle("Time [ns]")
    else:
        g_wf.GetXaxis().SetLabelSize(0)
        g_wf.GetXaxis().SetTitleSize(0)
        g_wf.GetYaxis().SetTitleOffset(0.9)
    g_wf.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    g_wf.GetXaxis().SetRangeUser(0, x_values[-1])
    g_wf.GetYaxis().SetRangeUser(-0.1, 1.1)
    g_wf.Draw("AL")

    g_shap = None
    if shap_mean is None:
        print(f"  SHAP cache not found in {SWEEP_CACHE_DIR}/, no SHAP overlay")
    else:
        shap_max = float(np.max(shap_mean))
        shap_norm = (shap_mean /
                     shap_max if shap_max > 0 else shap_mean).astype(
                         np.float64)
        n_shap = min(len(shap_norm), n_points)
        g_shap = ROOT.TGraph(n_shap, x_values[:n_shap], shap_norm[:n_shap])
        g_shap.SetLineColor(ROOT.kBlack)
        g_shap.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        g_shap.Draw("L SAME")

    g_outliers = None
    if len(outlier_orig_idx) > 0:
        outlier_t = (outlier_orig_idx.astype(np.float64) * 2)
        outlier_y = avg_wf_norm[outlier_orig_idx].astype(np.float64)
        g_outliers = ROOT.TGraph(len(outlier_orig_idx), outlier_t, outlier_y)
        g_outliers.SetMarkerStyle(marker_style)
        g_outliers.SetMarkerColor(marker_color)
        g_outliers.SetMarkerSize(1.4)
        g_outliers.Draw("P SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.5, 0.88, 0.5, 0.74)
    leg.AddEntry(g_wf, f"Average #{class_tag} Waveform", "l")
    if g_shap is not None:
        leg.AddEntry(g_shap, "SHAP / max(SHAP)", "l")
    if g_outliers is not None:
        leg.AddEntry(g_outliers,
                     f"Residual > {residual_threshold:g} ADC^{{2}}", "p")
    leg.SetMargin(0.1)
    leg.Draw()

    if pad_bot is not None:
        pad_top.SetTickx(0)
        pad_bot.cd()
        n_shap = min(len(shap_mean), n_points)
        cum = cumulative_shap(shap_mean[:n_shap])
        g_cum = ROOT.TGraph(n_shap, x_values[:n_shap], cum.astype(np.float64))
        g_cum.SetLineColor(ROOT.kBlack)
        g_cum.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        style_residual_pad_axes(g_cum, "Time [ns]")
        g_cum.GetXaxis().SetLimits(0, x_values[-1])
        g_cum.Draw("AL")

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _plot_quadratic_residual(test_data,
                             output_name,
                             p0,
                             p1,
                             p2,
                             class_tag,
                             marker_color,
                             marker_style=20,
                             min_amplitude_frac=0.03,
                             residual_threshold=500.0):
    """Subtract the interactively-fitted quadratic sigma^2_fit(mu) = p0 +
    p1*mu + p2*mu^2 from class_tag per-sample variance, then plot the
    residual against (dmu/dt)^2.

    Per-event timing jitter delta_t contributes Var(delta_t) * (dmu/dt)^2
    to sigma^2(t). If the quadratic absorbs only the floor, shot, and
    multiplicative terms, the residual vs (dmu/dt)^2 scatter should be
    linear through origin with slope = Var(delta_t).

    Samples whose residual exceeds residual_threshold are labeled with
    their time (sample_idx * 2 ns) on the main plot, and a companion plot
    of their residual vs max-normalized XGBoost SHAP is saved as
    "{output_name}_vs_shap.pdf".
    """
    mean_w = np.mean(test_data, axis=0)
    var_w = np.var(test_data, axis=0)

    thresh = min_amplitude_frac * float(mean_w.max())
    mask = (mean_w > thresh) & (var_w > 0)
    sample_indices_kept = np.where(mask)[0]
    mu = mean_w[mask].astype(np.float64)
    v = var_w[mask].astype(np.float64)

    fit_v = p0 + p1 * mu + p2 * mu**2
    residual = v - fit_v

    # dt = 2 ns/sample. Per-event jitter contributes Var(delta_t) *
    # (dmu/dt)^2 to sigma^2, so plotting vs (dmu/dt)^2 should be linear.
    dmu_dt = np.gradient(mean_w) / 2.0
    dx = (dmu_dt[mask]**2).astype(np.float64)

    above = residual > residual_threshold
    outlier_sample_indices = sample_indices_kept[above]

    print(f"  {class_tag.capitalize()} quadratic residual: "
          f"p0 = {p0}, p1 = {p1}, p2 = {p2}")
    print(f"  Residual range = [{float(residual.min()):.2f}, "
          f"{float(residual.max()):.2f}] ADC^2 over {len(dx)} samples")
    print(f"  {int(above.sum())} samples with residual > "
          f"{residual_threshold:g} ADC^2:")
    for i in np.where(above)[0]:
        s = int(sample_indices_kept[i])
        print(f"    sample {s:3d} (t = {s * 2:4d} ns): "
              f"residual = {float(residual[i]):.2f}, "
              f"(dmu/dt)^2 = {float(dx[i]):.3f}")

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas.SetLeftMargin(0.18)
    canvas.SetRightMargin(0.04)

    g = ROOT.TGraph(len(dx), dx, residual.astype(np.float64))
    g.SetName(f"{class_tag}_quadratic_residual_vs_dmu_dt_sq")
    g.SetTitle("")
    g.SetMarkerStyle(marker_style)
    g.SetMarkerSize(0.9)
    g.SetMarkerColor(marker_color)
    g.GetXaxis().SetTitle("(d#mu/dt)^{2} [(ADC/ns)^{2}]")
    g.GetYaxis().SetTitle(
        "#sigma^{2}(t) #minus (p_{0} + p_{1}#mu + p_{2}#mu^{2}) [ADC^{2}]")
    g.GetYaxis().SetTitleOffset(1.5)
    g.GetXaxis().SetNdivisions(505)
    g.Draw("AP")

    label_objs = []
    for i in np.where(above)[0]:
        s = int(sample_indices_kept[i])
        lbl = ROOT.TLatex(float(dx[i]), float(residual[i]), f" {s * 2}")
        lbl.SetTextSize(0.022)
        lbl.SetTextColor(marker_color)
        lbl.Draw()
        label_objs.append(lbl)

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")

    root_path = os.path.join(CACHE_DIR, f"{output_name}.root")
    out_file = ROOT.TFile(root_path, "RECREATE")
    g.Write()
    out_file.Close()
    print(f"Saved {root_path}")

    # Companion plot: residual vs max-normalized SHAP for the outlier samples.
    if not above.any():
        return
    shap_mean = _load_xgb_shap_mean()
    if shap_mean is None:
        print(f"  SHAP cache not found in {SWEEP_CACHE_DIR}/, skipping "
              f"residual-vs-SHAP companion plot")
        return
    shap_max = float(shap_mean.max())
    shap_norm = shap_mean / shap_max if shap_max > 0 else shap_mean

    outlier_residual = residual[above].astype(np.float64)
    outlier_shap = np.array([
        float(shap_norm[s]) if s < len(shap_norm) else float("nan")
        for s in outlier_sample_indices
    ],
                            dtype=np.float64)

    canvas2 = ROOT.PlottingUtils.GetConfiguredCanvas()
    canvas2.SetLeftMargin(0.19)
    canvas2.SetRightMargin(0.04)

    g_shap = ROOT.TGraph(len(outlier_sample_indices), outlier_shap,
                         outlier_residual)
    g_shap.SetName(f"{class_tag}_outlier_residual_vs_shap_norm")
    g_shap.SetTitle("")
    g_shap.SetMarkerStyle(marker_style)
    g_shap.SetMarkerSize(0.9)
    g_shap.SetMarkerColor(marker_color)
    g_shap.GetXaxis().SetTitle("SHAP / max(SHAP)")
    g_shap.GetYaxis().SetTitle(
        "#sigma^{2}(t) #minus (p_{0} + p_{1}#mu + p_{2}#mu^{2}) [ADC^{2}]")
    g_shap.GetYaxis().SetTitleOffset(1.4)
    g_shap.GetXaxis().SetLimits(0.0, 1.05)
    g_shap.Draw("AP")

    shap_label_objs = []
    for x, y, s in zip(outlier_shap, outlier_residual, outlier_sample_indices):
        lbl = ROOT.TLatex(float(x), float(y), f" {int(s) * 2}")
        lbl.SetTextSize(0.022)
        lbl.SetTextColor(marker_color)
        lbl.Draw()
        shap_label_objs.append(lbl)

    output_name_shap = f"{output_name}_vs_shap"
    ROOT.PlottingUtils.SaveFigure(canvas2, output_name_shap, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas2.Close()
    print(f"Saved {output_name_shap}")

    root_path_shap = os.path.join(CACHE_DIR, f"{output_name_shap}.root")
    out_file_shap = ROOT.TFile(root_path_shap, "RECREATE")
    g_shap.Write()
    out_file_shap.Close()
    print(f"Saved {root_path_shap}")


def _plot_shap_overlay(alpha_test, levels, shaps_per_level, level_name,
                       level_unit, output_name):
    """Overlay max-normalized SHAP curves -- one per noise level -- on top
    of the max-normalized alpha average waveform. Layout matches the
    feature-importance plots in parameter_study.py: gray waveform, colored
    overlays, linear axes, time on x (ns). Colors cycle through a fixed
    palette (cool -> warm) so low noise stays blue and high noise red.
    """
    mean_a = np.mean(alpha_test, axis=0)
    wf_max = float(np.max(mean_a))
    avg_wf_norm = mean_a / wf_max if wf_max > 0 else mean_a

    n_points = len(mean_a)
    x_values = np.arange(n_points, dtype=np.float64) * 2

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas(False)
    pad_top = ROOT.TPad("pad_top", "", 0.0, 0.3, 1.0, 1.0)
    pad_top.SetBottomMargin(0.04)
    pad_top.SetTopMargin(0.12)
    pad_top.Draw()
    pad_bot = ROOT.TPad("pad_bot", "", 0.0, 0.0, 1.0, 0.3)
    pad_bot.SetTopMargin(0.04)
    pad_bot.SetBottomMargin(0.35)
    pad_bot.Draw()
    pad_top.cd()

    g_wf = ROOT.TGraph(n_points, x_values, avg_wf_norm.astype(np.float64))
    g_wf.SetLineColor(ROOT.kGray + 2)
    g_wf.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_wf.SetTitle("")
    g_wf.GetXaxis().SetLabelSize(0)
    g_wf.GetXaxis().SetTitleSize(0)
    g_wf.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    g_wf.GetYaxis().SetTitleOffset(0.9)
    g_wf.GetXaxis().SetRangeUser(0, x_values[-1])
    g_wf.GetYaxis().SetRangeUser(-0.1, 1.1)
    g_wf.Draw("AL")

    palette = [
        ROOT.kAzure + 1, ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kOrange + 1,
        ROOT.kRed + 2, ROOT.kMagenta + 2, ROOT.kViolet + 2, ROOT.kPink + 2
    ]

    shap_graphs = []
    cum_graphs = []
    leg = ROOT.PlottingUtils.AddLegend(0.55, 0.88, 0.775, 0.325)
    leg.AddEntry(g_wf, "Average #alpha Waveform", "l")
    leg.SetMargin(0.15)
    for i, (level, shap_arr) in enumerate(zip(levels, shaps_per_level)):
        shap_max = float(np.max(shap_arr))
        shap_norm = (shap_arr / shap_max if shap_max > 0 else shap_arr)
        n_shap = min(len(shap_norm), n_points)
        color = palette[i % len(palette)]

        pad_top.cd()
        g_shap = ROOT.TGraph(n_shap, x_values[:n_shap],
                             shap_norm[:n_shap].astype(np.float64))
        g_shap.SetLineColor(color)
        g_shap.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        g_shap.Draw("L SAME")
        shap_graphs.append(g_shap)
        unit_str = f" {level_unit}" if level_unit else ""
        leg.AddEntry(g_shap, f"{level_name} = {level:g}{unit_str} SHAP", "l")

        pad_bot.cd()
        cum = cumulative_shap(shap_arr[:n_shap])
        g_cum = ROOT.TGraph(n_shap, x_values[:n_shap], cum.astype(np.float64))
        g_cum.SetLineColor(color)
        g_cum.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        if not cum_graphs:
            style_residual_pad_axes(g_cum, "Time [ns]")
            g_cum.GetXaxis().SetLimits(0, x_values[-1])
            g_cum.Draw("AL")
        else:
            g_cum.Draw("L SAME")
        cum_graphs.append(g_cum)

    pad_top.cd()
    leg.Draw()
    pad_top.SetTickx(0)

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _plot_two_curves(levels,
                     aucs_test_only,
                     aucs_test_only_err,
                     aucs_matched,
                     aucs_matched_err,
                     baseline_auc,
                     baseline_auc_err,
                     x_title,
                     output_name,
                     x_log,
                     rms_baseline=None):
    """AUC vs noise level with bootstrap error bars on every sweep point and
    +/-1 sigma translucent bands on the two baseline reference lines."""
    canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "", 1200,
                          600)

    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.72, 1.0)
    pad_plot.SetRightMargin(0.02)
    if x_log:
        pad_plot.SetLogx(True)
    pad_plot.Draw()

    pad_leg = ROOT.TPad("pad_leg", "", 0.72, 0.0, 1.0, 1.0)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_leg.Draw()

    pad_plot.cd()

    n = len(levels)
    x_arr = np.array(levels, dtype=np.float64)
    ex_arr = np.zeros(n, dtype=np.float64)

    g_test = ROOT.TGraphErrors(n, x_arr,
                               np.array(aucs_test_only, dtype=np.float64),
                               ex_arr,
                               np.array(aucs_test_only_err, dtype=np.float64))
    g_test.SetLineColor(ROOT.kBlue + 2)
    g_test.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_test.SetMarkerColor(ROOT.kBlue + 2)
    g_test.SetMarkerStyle(20)
    g_test.SetMarkerSize(1.2)

    g_match = ROOT.TGraphErrors(n, x_arr,
                                np.array(aucs_matched, dtype=np.float64),
                                ex_arr,
                                np.array(aucs_matched_err, dtype=np.float64))
    g_match.SetLineColor(ROOT.kRed + 2)
    g_match.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_match.SetMarkerColor(ROOT.kRed + 2)
    g_match.SetMarkerStyle(21)
    g_match.SetMarkerSize(1.2)

    all_y_lo = ([a - e for a, e in zip(aucs_test_only, aucs_test_only_err)] +
                [a - e for a, e in zip(aucs_matched, aucs_matched_err)] + [
                    baseline_auc - baseline_auc_err,
                ])
    all_y_hi = ([a + e for a, e in zip(aucs_test_only, aucs_test_only_err)] +
                [a + e for a, e in zip(aucs_matched, aucs_matched_err)] + [
                    baseline_auc + baseline_auc_err,
                ])
    y_min = min(all_y_lo)
    y_max = max(all_y_hi)
    pad = max(0.005, 0.2 * (y_max - y_min))

    g_test.SetTitle("")
    g_test.GetXaxis().SetTitle(x_title)
    g_test.GetYaxis().SetTitle("ROC AUC")
    g_test.GetYaxis().SetRangeUser(y_min - pad, y_max + pad)
    g_test.Draw("APL")
    g_match.Draw("PL SAME")

    # Baseline +/-1 sigma bands: 2-point TGraphErrors with constant y and ey
    # rendered with "3" (filled) draw option for the band, "L" for the line.
    band_x = np.array([x_arr[0], x_arr[-1]], dtype=np.float64)
    band_ex = np.zeros(2, dtype=np.float64)

    base_band = ROOT.TGraphErrors(
        2, band_x, np.array([baseline_auc, baseline_auc], dtype=np.float64),
        band_ex,
        np.array([baseline_auc_err, baseline_auc_err], dtype=np.float64))
    base_band.SetFillColorAlpha(ROOT.kGray + 2, 0.3)
    base_band.SetLineColor(ROOT.kGray + 2)
    base_band.SetLineStyle(2)
    base_band.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    base_band.Draw("L3 SAME")

    rms_line = None
    if rms_baseline is not None:
        rms_line = ROOT.TLine(rms_baseline, y_min - pad, rms_baseline,
                              y_max + pad)
        rms_line.SetLineColor(ROOT.kMagenta + 2)
        rms_line.SetLineStyle(3)
        rms_line.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        rms_line.Draw()

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.0, 0.95, 0.25, 0.85)
    leg.SetMargin(0.15)
    leg.AddEntry(g_test, "#splitline{Trained unmodified}{/ tested noisy}",
                 "lpe")
    leg.AddEntry(g_match, "#splitline{Trained noisy}{/ tested noisy}", "lpe")
    leg.AddEntry(base_band, "#splitline{Unmodified}{baseline}", "lf")
    if rms_line is not None:
        leg.AddEntry(rms_line,
                     "#splitline{Baseline RMS}{= %g ADC}" % rms_baseline, "l")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs("plots", exist_ok=True)
    results_cache = os.path.join(CACHE_DIR, f"results_{LO_TAG}.pkl")

    # Raw split is always loaded so sample-waveform plots can run regardless
    # of whether AUC results are cached. Loading from raw_split.npz is fast.
    (a_train, g_train, a_test, g_test, a_test_lo,
     g_test_lo) = _load_or_build_split()
    print(f"Loaded (LO {LO_LOWER['alpha']}-{LO_UPPER['alpha']}): "
          f"{len(a_train)} alpha + {len(g_train)} gamma train  |  "
          f"{len(a_test)} alpha + {len(g_test)} gamma test")

    a_train = _balanced_subsample(a_train, N_TRAIN_PER_CLASS, SUBSAMPLE_SEED)
    g_train = _balanced_subsample(g_train, N_TRAIN_PER_CLASS,
                                  SUBSAMPLE_SEED + 1)
    a_test, a_test_lo = _balanced_subsample_pair(a_test, a_test_lo,
                                                 N_TEST_PER_CLASS,
                                                 SUBSAMPLE_SEED + 2)
    g_test, g_test_lo = _balanced_subsample_pair(g_test, g_test_lo,
                                                 N_TEST_PER_CLASS,
                                                 SUBSAMPLE_SEED + 3)
    print(f"Balanced: {len(a_train)} alpha + {len(g_train)} gamma train  |  "
          f"{len(a_test)} alpha + {len(g_test)} gamma test")

    # Recompute if the cache predates the bootstrap-error refactor: the
    # cached dict will be missing the new *_err keys.
    required_keys = {
        "baseline_auc_err",
        "white_aucs_test_only_err",
        "white_aucs_matched_err",
        "shot_aucs_test_only_err",
        "shot_aucs_matched_err",
        "white_shaps_per_level",
        "shot_shaps_per_level",
        "white_shaps_test_only_per_level",
        "shot_shaps_test_only_per_level",
    }
    results = None
    if os.path.exists(results_cache):
        with open(results_cache, "rb") as fh:
            cached = pickle.load(fh)
        if required_keys.issubset(cached.keys()):
            print(f"Loading cached results from {results_cache}")
            results = cached
        else:
            print(f"Cache {results_cache} is missing bootstrap-error keys; "
                  f"recomputing")

    if results is None:
        if not os.path.exists(MODEL_PATH):
            raise FileNotFoundError(
                f"{MODEL_PATH} not found. Run analysis.py first to train "
                f"and cache the XGBoost regressor.")
        with open(MODEL_PATH, "rb") as fh:
            clean_model = pickle.load(fh)

        baseline_auc, baseline_auc_err = _score_clean_trained(
            clean_model, a_test, g_test, a_test_lo, g_test_lo)
        print(f"Baseline (unmodified, clean-trained) AUC = "
              f"{baseline_auc:.4f} +/- {baseline_auc_err:.4f}")

        print("White noise sweeps (sigma in ADC counts / LSB):")
        (white_test_only, white_test_only_err,
         white_test_only_shaps) = _sweep_test_only(clean_model,
                                                   a_test,
                                                   g_test,
                                                   a_test_lo,
                                                   g_test_lo,
                                                   WHITE_SIGMAS_ADC,
                                                   _add_white_noise,
                                                   "sigma_LSB",
                                                   compute_shap=True)
        white_matched, white_matched_err, white_shaps = _sweep_matched(
            a_train, g_train, a_test, g_test, a_test_lo, g_test_lo,
            WHITE_SIGMAS_ADC, _add_white_noise, "sigma_LSB")

        print("Shot noise sweeps (sigma = m * sqrt(|signal|)):")
        (shot_test_only, shot_test_only_err,
         shot_test_only_shaps) = _sweep_test_only(clean_model,
                                                  a_test,
                                                  g_test,
                                                  a_test_lo,
                                                  g_test_lo,
                                                  SHOT_MULTIPLIERS,
                                                  _add_shot_noise,
                                                  "multiplier",
                                                  compute_shap=True)
        shot_matched, shot_matched_err, shot_shaps = _sweep_matched(
            a_train, g_train, a_test, g_test, a_test_lo, g_test_lo,
            SHOT_MULTIPLIERS, _add_shot_noise, "multiplier")

        results = {
            "baseline_auc": baseline_auc,
            "baseline_auc_err": baseline_auc_err,
            "white_sigmas": WHITE_SIGMAS_ADC,
            "white_aucs_test_only": white_test_only,
            "white_aucs_test_only_err": white_test_only_err,
            "white_shaps_test_only_per_level": white_test_only_shaps,
            "white_aucs_matched": white_matched,
            "white_aucs_matched_err": white_matched_err,
            "white_shaps_per_level": white_shaps,
            "shot_multipliers": SHOT_MULTIPLIERS,
            "shot_aucs_test_only": shot_test_only,
            "shot_aucs_test_only_err": shot_test_only_err,
            "shot_shaps_test_only_per_level": shot_test_only_shaps,
            "shot_aucs_matched": shot_matched,
            "shot_aucs_matched_err": shot_matched_err,
            "shot_shaps_per_level": shot_shaps,
        }
        with open(results_cache, "wb") as fh:
            pickle.dump(results, fh)
        print(f"Results cached to {results_cache}")

    _plot_two_curves(results["white_sigmas"],
                     results["white_aucs_test_only"],
                     results["white_aucs_test_only_err"],
                     results["white_aucs_matched"],
                     results["white_aucs_matched_err"],
                     results["baseline_auc"],
                     results["baseline_auc_err"],
                     "White noise #sigma [ADC]",
                     "noise_white_auc",
                     x_log=False,
                     rms_baseline=BASELINE_RMS_ADC)

    _plot_two_curves(
        results["shot_multipliers"],
        results["shot_aucs_test_only"],
        results["shot_aucs_test_only_err"],
        results["shot_aucs_matched"],
        results["shot_aucs_matched_err"],
        results["baseline_auc"],
        results["baseline_auc_err"],
        "Shot noise multiplier #it{m} (#sigma = #it{m}#sqrt{|f(t)|})",
        "noise_shot_auc",
        x_log=False)

    print("Plotting SHAP overlay across white-noise levels (matched)...")
    _plot_shap_overlay(a_test, results["white_sigmas"],
                       results["white_shaps_per_level"], "#sigma", "ADC",
                       "shap_overlay_white_matched")

    print("Plotting SHAP overlay across shot-noise levels (matched)...")
    _plot_shap_overlay(a_test, results["shot_multipliers"],
                       results["shot_shaps_per_level"], "#it{m}", "",
                       "shap_overlay_shot_matched")

    print("Plotting SHAP overlay across white-noise levels (test-only)...")
    _plot_shap_overlay(a_test, results["white_sigmas"],
                       results["white_shaps_test_only_per_level"], "#sigma",
                       "ADC", "shap_overlay_white_testonly")

    print("Plotting SHAP overlay across shot-noise levels (test-only)...")
    _plot_shap_overlay(a_test, results["shot_multipliers"],
                       results["shot_shaps_test_only_per_level"], "#it{m}", "",
                       "shap_overlay_shot_testonly")

    print("Plotting sample waveforms at each white-noise level...")
    _plot_noise_sample_waveforms(a_test, g_test, results["white_sigmas"],
                                 _add_white_noise, "sigma", "LSB",
                                 "noise_waveforms_white")

    print("Plotting sample waveforms at each shot-noise level...")
    _plot_noise_sample_waveforms(a_test, g_test, results["shot_multipliers"],
                                 _add_shot_noise, "m", "",
                                 "noise_waveforms_shot")

    print("Filtered train+test AUC: Butterworth LPF sweep...")
    lpf_aucs, lpf_errs = _sweep_filter(a_train, g_train, a_test, g_test,
                                       a_test_lo, g_test_lo,
                                       LPF_CUTOFFS_MHZ, _apply_butter_lpf,
                                       "butter")
    _plot_filter_auc(
        LPF_CUTOFFS_MHZ, lpf_aucs, lpf_errs, results["baseline_auc"],
        results["baseline_auc_err"], "Butterworth LPF cutoff [MHz]",
        f"#splitline{{Butterworth LPF}}"
        f"{{(order {LPF_ORDER}, train+test)}}", ROOT.kBlue + 2, 20,
        "filter_auc_butter")

    print("Filtered train+test AUC: Gaussian smoothing sweep...")
    gauss_aucs, gauss_errs = _sweep_filter(a_train, g_train, a_test, g_test,
                                           a_test_lo, g_test_lo,
                                           GAUSSIAN_SIGMAS_SAMPLES,
                                           _apply_gaussian, "gauss")
    _plot_filter_auc(GAUSSIAN_SIGMAS_SAMPLES, gauss_aucs, gauss_errs,
                     results["baseline_auc"], results["baseline_auc_err"],
                     "Gaussian smoothing #sigma [samples]",
                     "#splitline{Gaussian smoothing}{(train+test)}",
                     ROOT.kRed + 2, 21, "filter_auc_gauss")

    print(f"Notch experiment at {NOTCH_FREQS_MHZ} MHz (Q={NOTCH_Q})...")
    notch_results = _run_notch_experiment(a_train, g_train, a_test, g_test,
                                          a_test_lo, g_test_lo)
    _plot_notch_results(notch_results, "notch_auc")

    print("Plotting LO spectra at each white-noise level...")
    _plot_noise_light_output(a_test, g_test, results["white_sigmas"],
                             _add_white_noise, "sigma", "LSB",
                             "noise_lo_white")

    print("Plotting LO spectra at each shot-noise level...")
    _plot_noise_light_output(a_test, g_test, results["shot_multipliers"],
                             _add_shot_noise, "m", "", "noise_lo_shot")

    print("Plotting per-sample variance vs mean (shot-noise diagnostic)...")
    _plot_variance_vs_mean(a_test, g_test, "variance_vs_mean")

    print("Plotting normalized FFT of average alpha/gamma waveforms...")
    _plot_avg_waveform_fft(a_test, g_test, "avg_waveform_fft")
    _plot_avg_waveform_fft(a_test,
                           g_test,
                           "avg_waveform_fft_power",
                           power=True)
    _plot_avg_waveform_fft(a_test,
                           g_test,
                           "avg_waveform_fft_power_zoom",
                           power=True,
                           f_min_mhz=50.0)

    print("Plotting ensemble-averaged per-event residual noise PSD...")
    _plot_residual_psd(a_test, g_test, "residual_psd")
    _plot_residual_psd(a_test, g_test, "residual_psd_zoom", f_min_mhz=50.0)

    print("Plotting per-class sample-value distributions at fixed times...")
    _plot_sample_value_distributions(a_test, g_test,
                                     "sample_value_distributions")

    print("Plotting variance vs |dmu/dt| (direct jitter diagnostic)...")
    _plot_variance_vs_dmu_dt(a_test, g_test, "variance_vs_dmu_dt")

    print("Plotting alpha quadratic-fit residual vs (dmu/dt)^2...")
    _plot_quadratic_residual(a_test,
                             "alpha_quadratic_residual",
                             ALPHA_VAR_QUAD_P0,
                             ALPHA_VAR_QUAD_P1,
                             ALPHA_VAR_QUAD_P2,
                             "alpha",
                             ROOT.kRed + 2,
                             marker_style=20)

    print("Plotting gamma quadratic-fit residual vs (dmu/dt)^2...")
    _plot_quadratic_residual(g_test,
                             "gamma_quadratic_residual",
                             GAMMA_VAR_QUAD_P0,
                             GAMMA_VAR_QUAD_P1,
                             GAMMA_VAR_QUAD_P2,
                             "gamma",
                             ROOT.kBlue + 2,
                             marker_style=21)

    print("Plotting alpha average waveform with residual outliers + SHAP...")
    _plot_avg_waveform_outliers(a_test,
                                "alpha_avg_waveform_outliers",
                                ALPHA_VAR_QUAD_P0,
                                ALPHA_VAR_QUAD_P1,
                                ALPHA_VAR_QUAD_P2,
                                "alpha",
                                ROOT.kRed + 2,
                                marker_style=20)

    print("Plotting gamma average waveform with residual outliers + SHAP...")
    _plot_avg_waveform_outliers(g_test,
                                "gamma_avg_waveform_outliers",
                                GAMMA_VAR_QUAD_P0,
                                GAMMA_VAR_QUAD_P1,
                                GAMMA_VAR_QUAD_P2,
                                "gamma",
                                ROOT.kBlue + 2,
                                marker_style=21)


if __name__ == "__main__":
    main()
