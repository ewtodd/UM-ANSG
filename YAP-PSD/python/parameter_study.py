import argparse
import os
import pickle
import numpy as np
import ROOT
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.metrics import make_scorer, roc_auc_score
from scipy.stats import randint, uniform, loguniform
from joblib import parallel_backend
from xgboost import XGBRegressor
from analysis_utilities.io import load_tree_data
from psd_utils import (process_waveforms, ROOT_FILES_DIR, cumulative_shap,
                       style_residual_pad_axes, compute_mean_abs_shap,
                       per_lo_auc_from_scores, error_weighted_auc)
from torch_models import TorchMLPRegressor
import analysis_utilities

analysis_utilities.load_cpp_library()
ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPDF)

CACHE_DIR = "sweep_cache"

N_SEEDS = 3
#SEEDS = [42, 123, 256, 789, 1024][:N_SEEDS]
SEEDS = [42, 123, 256][:N_SEEDS]
DEFAULT_TRAIN_PER_CLASS = 10000

RANDOM_SEARCH_N_ITER = 50
RANDOM_SEARCH_CV_FOLDS = 3

RF_CONFIG = dict(
    name="Random Forest",
    prefix="rf",
    model_class=RandomForestRegressor,
    explainer="tree",
    color=ROOT.kRed + 2,
    default_params=dict(
        n_estimators=100,
        max_depth=None,
        max_samples=None,
        max_features=1.0,
        n_jobs=-1,
    ),
    sweeps=[
        dict(
            sweep_name="n_estimators",
            values=[5, 10, 25, 50, 100, 150, 200, 300],
            x_title="Number of Trees",
            param_key="n_estimators",
        ),
        dict(
            sweep_name="n_training_samples",
            values=[500, 1000, 10000, 25000, 50000, 100000],
            x_title="Training Samples per Class",
            param_key=None,  # special handling
        ),
        dict(
            sweep_name="max_depth",
            values=[3, 5, 10, 20, 30, None],
            x_title="Max Depth (50 = None)",
            param_key="max_depth",
        ),
        dict(
            sweep_name="max_samples",
            values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.632, 0.8, 1.0],
            x_title="Max Samples (Bootstrap Fraction)",
            param_key="max_samples",
        ),
    ],
)

GB_CONFIG = dict(
    name="Gradient Boosting",
    prefix="gb",
    model_class=GradientBoostingRegressor,
    explainer="tree",
    color=ROOT.kBlue + 2,
    default_params=dict(
        n_estimators=200,
        max_depth=5,
        learning_rate=0.12256,
        verbose=1,
    ),
    sweeps=[],
)

MLP_CONFIG = dict(
    name="MLP",
    prefix="mlp",
    model_class=TorchMLPRegressor,
    explainer="kernel",
    color=ROOT.kOrange + 2,
    default_params=dict(
        hidden_layer_sizes=(128, 64),
        max_iter=500,
    ),
    sweeps=[],
)

# Mapping from sweep prefixes to the regressor names registered in
# regressors.py; used by --skip-randomized-search to pull tuned params from
# the single source of truth instead of duplicating them here.
PREFIX_TO_REGRESSOR_NAME = {"rf": "Random Forest", "xgb": "XGBoost"}


def _best_params_from_regressors(config):
    """Pull the tuned parameter values for this config from regressors.py.

    Returns only the params that this config actually sweeps (per
    config["sweeps"]), so the result is shape-compatible with
    RandomizedSearchCV.best_params_.
    """
    from regressors import get_default_regressors
    target = PREFIX_TO_REGRESSOR_NAME[config["prefix"]]
    for reg in get_default_regressors():
        if reg["name"] == target:
            full_params = reg["model"].get_params()
            tunable = [
                s["param_key"] for s in config["sweeps"]
                if s["param_key"] is not None
            ]
            return {k: full_params[k] for k in tunable if k in full_params}
    raise KeyError(
        f"No regressor named {target!r} in get_default_regressors()")


XGB_CONFIG = dict(
    name="XGBoost",
    prefix="xgb",
    model_class=XGBRegressor,
    explainer="tree",
    color=ROOT.kGreen + 2,
    default_params=dict(
        n_estimators=100,
        max_depth=6,
        learning_rate=0.3,
        n_jobs=-1,
        verbosity=1,
    ),
    sweeps=[
        dict(
            sweep_name="n_estimators",
            values=[1, 5, 10, 25, 50, 100, 150, 200, 300, 500],
            x_title="Number of Boosting Rounds",
            param_key="n_estimators",
        ),
        dict(
            sweep_name="n_training_samples",
            values=[500, 1000, 10000, 25000, 50000, 100000],
            x_title="Training Samples per Class",
            param_key=None,
        ),
        dict(
            sweep_name="max_depth",
            values=[1, 2, 3, 5, 7, 10, 15],
            x_title="Max Depth",
            param_key="max_depth",
        ),
        dict(
            sweep_name="learning_rate",
            values=[0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
            x_title="Learning Rate",
            param_key="learning_rate",
        ),
    ],
)


def _prepare_data(alpha_waveforms,
                  gamma_waveforms,
                  alpha_features,
                  gamma_features,
                  n_train_per_class=10000,
                  random_state=42,
                  return_lo=False):
    """Filter by light output, balance classes, split train/test, and process waveforms.

    When return_lo=True, also returns test_lo: a 1-D array of length len(x_test)
    holding the light-output of each test event (aligned with x_test rows).
    """
    lo_lower = {"alpha": 375, "gamma": 0}
    lo_upper = {"alpha": 1575, "gamma": 1750}

    alpha_mask_idx = np.where(
        (alpha_features["light_output"] >= lo_lower["alpha"])
        & (alpha_features["light_output"] <= lo_upper["alpha"]))[0]
    gamma_mask_idx = np.where(
        (gamma_features["light_output"] >= lo_lower["gamma"])
        & (gamma_features["light_output"] <= lo_upper["gamma"]))[0]

    alpha_masked = alpha_waveforms[alpha_mask_idx]
    gamma_masked = gamma_waveforms[gamma_mask_idx]

    min_samples = min(n_train_per_class, len(alpha_masked), len(gamma_masked))
    print(f"Balanced training pool size per class: {min_samples}")

    rng = np.random.RandomState(random_state)
    alpha_train_sel = rng.choice(len(alpha_masked),
                                 size=min_samples,
                                 replace=False)
    gamma_train_sel = rng.choice(len(gamma_masked),
                                 size=min_samples,
                                 replace=False)

    train_alpha_wf = process_waveforms(alpha_masked[alpha_train_sel])
    train_gamma_wf = process_waveforms(gamma_masked[gamma_train_sel])

    x_train = np.vstack([train_alpha_wf, train_gamma_wf])
    y_train = np.array([0] * len(train_alpha_wf) + [1] * len(train_gamma_wf))

    alpha_train_set = set(alpha_train_sel.tolist())
    gamma_train_set = set(gamma_train_sel.tolist())

    alpha_test_sel = [
        i for i in range(len(alpha_masked)) if i not in alpha_train_set
    ]
    gamma_test_sel = [
        i for i in range(len(gamma_masked)) if i not in gamma_train_set
    ]

    test_alpha_wf = process_waveforms(alpha_masked[alpha_test_sel])
    test_gamma_wf = process_waveforms(gamma_masked[gamma_test_sel])

    x_test = np.vstack([test_alpha_wf, test_gamma_wf])
    y_test = np.array([0] * len(test_alpha_wf) + [1] * len(test_gamma_wf))

    print(
        f"Train: {len(x_train)} ({len(train_alpha_wf)} alpha + {len(train_gamma_wf)} gamma)"
    )
    print(
        f"Test:  {len(x_test)} ({len(test_alpha_wf)} alpha + {len(test_gamma_wf)} gamma)"
    )

    if return_lo:
        test_alpha_lo = alpha_features["light_output"].values[
            alpha_mask_idx[alpha_test_sel]]
        test_gamma_lo = gamma_features["light_output"].values[
            gamma_mask_idx[gamma_test_sel]]
        test_lo = np.concatenate([test_alpha_lo, test_gamma_lo])
        return x_train, y_train, x_test, y_test, test_lo

    return x_train, y_train, x_test, y_test


def _lo_averaged_auc(y_test, scores, test_lo):
    """Bin scores by light output and return (mean, err) via inverse-variance
    weighting. Replaces single-pool bootstrap_auc so reported AUC is not
    inflated by the gamma/alpha population imbalance."""
    alpha_mask = (y_test == 0)
    gamma_mask = (y_test == 1)
    aucs, errs = per_lo_auc_from_scores(scores[alpha_mask],
                                        scores[gamma_mask],
                                        test_lo[alpha_mask],
                                        test_lo[gamma_mask])
    return error_weighted_auc(aucs, errs)


def _plot_sweep(x_values, y_values, y_errors, x_title, prefix, output_name,
                color):
    """Plot a single sweep with colored line + black markers/error bars."""
    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    x_arr = np.array(x_values, dtype=np.float64)
    y_arr = np.array(y_values, dtype=np.float64)
    ex_arr = np.zeros(len(x_arr), dtype=np.float64)
    ey_arr = np.array(y_errors, dtype=np.float64)

    graph_line = ROOT.TGraph(len(x_arr), x_arr, y_arr)
    graph_line.SetLineColor(color)
    graph_line.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

    graph_err = ROOT.TGraphErrors(len(x_arr), x_arr, y_arr, ex_arr, ey_arr)
    graph_err.SetLineColor(ROOT.kBlack)
    graph_err.SetMarkerColor(ROOT.kBlack)
    graph_err.SetMarkerStyle(20)
    graph_err.SetMarkerSize(1.2)

    graph_line.SetTitle("")
    if x_title == "Training Samples per Class":
        canvas.SetLogx(True)
    if "Max Samples" in x_title:
        graph_line.GetYaxis().SetTitleOffset(1.7)
        canvas.SetLeftMargin(0.2)
    if x_title == "Max Depth":
        graph_line.GetYaxis().SetTitleOffset(1.4)

    graph_line.GetXaxis().SetTitle(x_title)
    graph_line.GetYaxis().SetTitle("ROC AUC")

    y_min = min(y_values)
    y_max = max(y_values)
    y_range = y_max - y_min if y_max > y_min else 0.01
    graph_line.GetYaxis().SetRangeUser(y_min - 0.5 * y_range,
                                       y_max + 0.2 * y_range)

    graph_line.Draw("AL")
    graph_err.Draw("P SAME")

    if (prefix == "rf"):
        _ = ROOT.PlottingUtils.AddText("Random Forest", 0.85, 0.25)
    else:
        _ = ROOT.PlottingUtils.AddText("XGBoost", 0.85, 0.25)

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _run_sweep(config, sweep_name, values, param_key, x_train, y_train, x_test,
               y_test, test_lo):
    """Run a parameter sweep, caching results to disk."""
    prefix = config["prefix"]
    model_class = config["model_class"]
    default_params = config["default_params"]
    n_per_class = len(x_train) // 2

    cache_file = os.path.join(CACHE_DIR, f"{prefix}_sweep_{sweep_name}.npz")
    if os.path.exists(cache_file):
        print(f"  Loading cached results: {cache_file}")
        data = np.load(cache_file)
        means = data["means"].tolist()
        stds = data["stds"].tolist()
        for val, m, s in zip(values, means, stds):
            print(f"  {sweep_name}={val}    AUC = {m:.4f} +/- {s:.4f}")
        return means, stds

    # For regular sweeps, subsample to DEFAULT_TRAIN_PER_CLASS so we don't
    # train with the entire (potentially huge) pool meant for the
    # n_training_samples sweep.
    if param_key is not None and n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(42)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        default_idx = np.concatenate([alpha_sel, gamma_sel])
        x_train_sweep = x_train[default_idx]
        y_train_sweep = y_train[default_idx]
        print(f"  Subsampled to {len(x_train_sweep)} training samples "
              f"({DEFAULT_TRAIN_PER_CLASS} per class)")
    else:
        x_train_sweep = x_train
        y_train_sweep = y_train
        print(f"  Using full training pool: {len(x_train_sweep)} samples")

    means, stds = [], []
    for val in values:
        print(f"  {sweep_name}={val}")
        params = dict(default_params, random_state=42)

        if param_key is not None:
            params[param_key] = val
            model = model_class(**params)
            model.fit(x_train_sweep, y_train_sweep)
        else:
            # n_training_samples: subsample from full training pool
            rng = np.random.RandomState(42)
            alpha_sel = rng.choice(n_per_class, size=val, replace=False)
            gamma_sel = rng.choice(n_per_class, size=val,
                                   replace=False) + n_per_class
            idx = np.concatenate([alpha_sel, gamma_sel])
            model = model_class(**params)
            model.fit(x_train[idx], y_train[idx])

        scores = model.predict(x_test)
        auc_mean, auc_std = _lo_averaged_auc(y_test, scores, test_lo)
        means.append(auc_mean)
        stds.append(auc_std)
        print(f"    AUC = {auc_mean:.4f} +/- {auc_std:.4f}")

    np.savez(cache_file, means=np.array(means), stds=np.array(stds))
    return means, stds


def _plot_feature_importance(model_class, default_params, prefix, color, name,
                             x_train, y_train, avg_waveform):
    """Train or load N_SEEDS models and plot averaged feature importance."""
    # Subsample to DEFAULT_TRAIN_PER_CLASS if pool is larger
    n_per_class = len(x_train) // 2
    if n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(0)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        idx = np.concatenate([alpha_sel, gamma_sel])
        x_train = x_train[idx]
        y_train = y_train[idx]

    all_importances = []
    for seed in SEEDS:
        model_file = os.path.join(CACHE_DIR,
                                  f"{prefix}_importance_seed{seed}.pkl")
        if os.path.exists(model_file):
            print(f"  Loading cached model: {model_file}")
            with open(model_file, "rb") as fh:
                model = pickle.load(fh)
        else:
            print(f"  Training seed {seed}...")
            params = dict(default_params, random_state=seed)
            model = model_class(**params)
            model.fit(x_train, y_train)
            with open(model_file, "wb") as fh:
                pickle.dump(model, fh)
        all_importances.append(model.feature_importances_)

    importances = np.array(all_importances)
    mean_imp = np.mean(importances, axis=0)
    std_imp = np.std(importances, axis=0)

    wf_max = np.max(avg_waveform)
    avg_wf_norm = avg_waveform / wf_max if wf_max > 0 else avg_waveform

    imp_max = np.max(mean_imp)
    mean_imp_norm = mean_imp / imp_max if imp_max > 0 else mean_imp
    std_imp_norm = std_imp / imp_max if imp_max > 0 else std_imp

    n_points = len(avg_waveform)
    x_values = np.arange(n_points, dtype=np.float64) * 2

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    graph_waveform = ROOT.TGraph(n_points, x_values,
                                 avg_wf_norm.astype(np.float64))
    graph_waveform.SetLineColor(ROOT.kGray + 2)
    graph_waveform.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_waveform.SetTitle("")
    graph_waveform.GetXaxis().SetTitle("Time [ns]")
    graph_waveform.GetYaxis().SetTitleOffset(1.0)
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(-0.1, 1.1)
    graph_waveform.Draw("AL")

    ex = np.zeros(n_points, dtype=np.float64)
    graph_importance = ROOT.TGraphErrors(n_points, x_values,
                                         mean_imp_norm.astype(np.float64), ex,
                                         std_imp_norm.astype(np.float64))
    graph_importance.SetLineColor(color)
    graph_importance.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_importance.SetFillColorAlpha(color, 0.2)
    graph_importance.Draw("L3 SAME")

    if (name == "Random Forest" or name == "Gradient Boosting"):
        leg = ROOT.PlottingUtils.AddLegend(0.42, 0.88, 0.6, 0.85)
    else:
        leg = ROOT.PlottingUtils.AddLegend(0.5, 0.85, 0.6, 0.85)
    leg.AddEntry(graph_waveform, "Average #alpha Waveform", "l")
    leg.AddEntry(graph_importance, f"{name} Feature Importance", "lf")
    leg.SetMargin(0.1)
    leg.Draw()

    output_file = f"feature_importance_{prefix}"
    ROOT.PlottingUtils.SaveFigure(canvas, output_file, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_file}")


def _plot_shap_importance(model_class, default_params, prefix, color, name,
                          x_train, y_train, avg_waveform, explainer_type):
    """Train or load N_SEEDS models and plot SHAP-based feature importance.

    explainer_type="tree" uses interventional TreeSHAP (exact for tree
    models). explainer_type="kernel" uses KernelSHAP (model-agnostic
    sampling). Both are run with the same background and explain sets so
    tree-model and MLP importances are estimating the same interventional
    Shapley values.
    """
    n_per_class = len(x_train) // 2
    if n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(0)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        idx = np.concatenate([alpha_sel, gamma_sel])
        x_train = x_train[idx]
        y_train = y_train[idx]

    all_mean_abs_shap = []
    for seed in SEEDS:
        shap_cache = os.path.join(
            CACHE_DIR, f"{prefix}_{explainer_type}_shap_values_seed{seed}.npy")

        if os.path.exists(shap_cache):
            print(f"  Loading cached SHAP values: {shap_cache}")
            mean_abs_shap = np.load(shap_cache)
        else:
            model_file = os.path.join(CACHE_DIR,
                                      f"{prefix}_shap_seed{seed}.pkl")
            if os.path.exists(model_file):
                print(f"  Loading cached model: {model_file}")
                with open(model_file, "rb") as fh:
                    model = pickle.load(fh)
            else:
                print(f"  Training {name} seed {seed}...")
                params = dict(default_params, random_state=seed)
                model = model_class(**params)
                model.fit(x_train, y_train)
                with open(model_file, "wb") as fh:
                    pickle.dump(model, fh)

            print(f"  Computing {explainer_type} SHAP values "
                  f"for seed {seed}...")
            mean_abs_shap = compute_mean_abs_shap(model, x_train, seed,
                                                  explainer_type)
            np.save(shap_cache, mean_abs_shap)

        all_mean_abs_shap.append(mean_abs_shap)

    importances = np.array(all_mean_abs_shap)
    mean_imp = np.mean(importances, axis=0)
    std_imp = np.std(importances, axis=0)

    wf_max = np.max(avg_waveform)
    avg_wf_norm = avg_waveform / wf_max if wf_max > 0 else avg_waveform

    imp_max = np.max(mean_imp)
    mean_imp_norm = mean_imp / imp_max if imp_max > 0 else mean_imp
    std_imp_norm = std_imp / imp_max if imp_max > 0 else std_imp

    n_points = len(avg_waveform)
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

    graph_waveform = ROOT.TGraph(n_points, x_values,
                                 avg_wf_norm.astype(np.float64))
    graph_waveform.SetLineColor(ROOT.kGray + 2)
    graph_waveform.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_waveform.SetTitle("")
    graph_waveform.GetXaxis().SetLabelSize(0)
    graph_waveform.GetXaxis().SetTitleSize(0)
    graph_waveform.GetYaxis().SetTitleOffset(0.9)
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(-0.1, 1.1)
    graph_waveform.Draw("AL")

    ex = np.zeros(n_points, dtype=np.float64)
    graph_importance = ROOT.TGraphErrors(n_points, x_values,
                                         mean_imp_norm.astype(np.float64), ex,
                                         std_imp_norm.astype(np.float64))
    graph_importance.SetLineColor(color)
    graph_importance.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_importance.SetFillColorAlpha(color, 0.2)
    graph_importance.Draw("L3 SAME")

    if (name == "Random Forest" or name == "Gradient Boosting"):
        leg = ROOT.PlottingUtils.AddLegend(0.56, 0.85, 0.6, 0.85)
    else:
        leg = ROOT.PlottingUtils.AddLegend(0.59, 0.85, 0.6, 0.85)
    leg.AddEntry(graph_waveform, "Average #alpha Waveform", "l")
    leg.AddEntry(graph_importance, f"{name} SHAP", "lf")
    leg.SetMargin(0.1)
    leg.Draw()
    pad_top.SetTickx(0)

    pad_bot.cd()
    cum = cumulative_shap(mean_imp)
    graph_cum = ROOT.TGraph(n_points, x_values, cum.astype(np.float64))
    graph_cum.SetLineColor(color)
    graph_cum.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    style_residual_pad_axes(graph_cum, "Time [ns]")
    graph_cum.GetXaxis().SetLimits(0, x_values[-1])
    graph_cum.Draw("AL")

    output_file = f"feature_importance_{prefix}_shap"
    ROOT.PlottingUtils.SaveFigure(canvas, output_file, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_file}")


def _recommend_value(values, means, default_value=None):
    """Find the recommended parameter value from a sweep.

    If the AUC curve has a clear peak, return the peak value.
    Otherwise, find the plateau onset: the first value where AUC >= 0.99.

    Then compare against the library default: if the default gives higher AUC
    than the plateau onset, prefer the default.
    """
    best_idx = int(np.argmax(means))
    best_auc = means[best_idx]

    # Check if it's a peak (AUC drops on both sides) vs a plateau
    is_peak = (best_idx > 0 and best_idx < len(means) - 1)

    if is_peak:
        rec_val, rec_auc = values[best_idx], best_auc
    else:
        # Plateau: first value that crosses 0.99 AUC
        rec_val, rec_auc = values[best_idx], best_auc
        for i, m in enumerate(means):
            if m >= 0.99:
                rec_val, rec_auc = values[i], m
                break

    # Compare against library default if available
    if default_value is not None and default_value in values:
        default_idx = values.index(default_value)
        default_auc = means[default_idx]
        if default_auc > rec_auc:
            print(
                f"  >> Plateau onset at {rec_val} (AUC={rec_auc:.4f}), "
                f"but default {default_value} is better (AUC={default_auc:.4f})"
            )
            rec_val, rec_auc = default_value, default_auc

    return rec_val, rec_auc


def _run_all_sweeps(config, x_train, y_train, x_test, y_test, test_lo):
    """Run all sweeps for a given model configuration.

    Returns a dict of recommended parameter values (param_key -> value),
    excluding n_training_samples.
    """
    color = config["color"]
    n_per_class = len(x_train) // 2
    recommended = {}

    for sweep in config["sweeps"]:
        sweep_name = sweep["sweep_name"]
        values = sweep["values"]
        x_title = sweep["x_title"]
        param_key = sweep["param_key"]
        output_name = f"{config['prefix']}_auc_vs_{sweep_name}"

        if sweep_name == "n_training_samples":
            values = [v for v in values if v <= n_per_class]

        print(f"{config['name']}: AUC vs {sweep_name}")

        means, stds = _run_sweep(config, sweep_name, values, param_key,
                                 x_train, y_train, x_test, y_test, test_lo)

        default_value = config["default_params"].get(
            param_key) if param_key else None
        rec_val, rec_auc = _recommend_value(values, means, default_value)
        print(
            f"  >> Recommended {sweep_name} = {rec_val} (AUC = {rec_auc:.4f})")

        if param_key is not None:
            recommended[param_key] = rec_val

        plot_x = values
        if sweep_name == "max_depth" and None in values:
            plot_x = [d if d is not None else 50 for d in values]

        _plot_sweep(plot_x, means, stds, x_title, config['prefix'],
                    output_name, color)

    return recommended


def _build_search_space(config, recommended):
    """Build RandomizedSearchCV parameter distributions around recommended values.

    For each tunable parameter, creates a distribution centered on the
    recommended value from the OAT sweeps, covering a reasonable neighborhood
    to explore interactions.
    """
    distributions = {}

    for sweep in config["sweeps"]:
        param_key = sweep["param_key"]
        if param_key is None:  # skip n_training_samples
            continue

        rec_val = recommended.get(param_key)
        if rec_val is None and param_key == "max_depth":
            # RF recommended unlimited depth; include None + some large values
            distributions[param_key] = [20, 30, 40, 50, None]
            continue
        if rec_val is None:
            continue

        values = [v for v in sweep["values"] if v is not None]

        if param_key == "n_estimators":
            low = max(10, int(rec_val * 0.5))
            high = int(rec_val * 2.0)
            distributions[param_key] = randint(low, high + 1)

        elif param_key == "max_depth":
            low = max(1, rec_val - 3)
            high = rec_val + 5
            distributions[param_key] = randint(low, high + 1)

        elif param_key == "learning_rate":
            low = max(1e-4, rec_val / 5.0)
            high = min(1.0, rec_val * 5.0)
            distributions[param_key] = loguniform(low, high)

        elif param_key == "max_samples":
            low = max(0.05, rec_val - 0.2)
            high = min(1.0, rec_val + 0.2)
            distributions[param_key] = uniform(low, high - low)

        elif param_key == "max_features":
            low = max(0.1, rec_val - 0.3)
            high = min(1.0, rec_val + 0.3)
            distributions[param_key] = uniform(low, high - low)

    return distributions


def _run_randomized_search(config, recommended, x_train, y_train, x_test,
                           y_test, test_lo):
    """Run RandomizedSearchCV around the OAT-recommended values.

    Uses the training set for cross-validated search, then evaluates the
    best parameters on the held-out test set with LO-bin-averaged AUC for a
    direct comparison against the OAT results.
    """
    prefix = config["prefix"]
    model_class = config["model_class"]
    default_params = config["default_params"]

    cache_file = os.path.join(CACHE_DIR, f"{prefix}_randomized_search.pkl")
    if os.path.exists(cache_file):
        print(f"  Loading cached RandomizedSearchCV: {cache_file}")
        with open(cache_file, "rb") as fh:
            search = pickle.load(fh)
        print(f"  Best CV score (ROC AUC): {search.best_score_:.4f}")
        print(f"  Best params: {search.best_params_}")

        # Evaluate on test set
        scores = search.best_estimator_.predict(x_test)
        auc_mean, auc_std = _lo_averaged_auc(y_test, scores, test_lo)
        print(f"  Test AUC: {auc_mean:.4f} +/- {auc_std:.4f}")
        return search.best_params_, search.best_score_

    # Subsample training data to DEFAULT_TRAIN_PER_CLASS per class
    n_per_class = len(x_train) // 2
    if n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(42)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        idx = np.concatenate([alpha_sel, gamma_sel])
        x_train_search = x_train[idx]
        y_train_search = y_train[idx]
        print(f"  Subsampled to {len(x_train_search)} training samples "
              f"({DEFAULT_TRAIN_PER_CLASS} per class)")
    else:
        x_train_search = x_train
        y_train_search = y_train

    # Build search space around recommended values
    distributions = _build_search_space(config, recommended)
    print(f"  Search distributions:")
    for k, v in distributions.items():
        print(f"    {k}: {v}")

    # Fixed params that aren't being searched
    fixed_params = {
        k: v
        for k, v in default_params.items()
        if k not in distributions and k != "random_state"
    }

    # Scorer: use regression predictions directly for ROC AUC
    auc_scorer = make_scorer(roc_auc_score)

    base_model = model_class(random_state=42, **fixed_params)

    cv_strategy = StratifiedKFold(n_splits=RANDOM_SEARCH_CV_FOLDS,
                                  shuffle=True,
                                  random_state=42)

    search = RandomizedSearchCV(
        base_model,
        param_distributions=distributions,
        n_iter=RANDOM_SEARCH_N_ITER,
        scoring=auc_scorer,
        cv=cv_strategy,
        random_state=42,
        n_jobs=-1,
        verbose=1,
        refit=True,
    )

    print(f"  Running RandomizedSearchCV ({RANDOM_SEARCH_N_ITER} iterations, "
          f"{RANDOM_SEARCH_CV_FOLDS}-fold CV)...")
    with parallel_backend("threading"):
        search.fit(x_train_search, y_train_search)

    with open(cache_file, "wb") as fh:
        pickle.dump(search, fh)

    print(f"  Best CV score (ROC AUC): {search.best_score_:.4f}")
    print(f"  Best params: {search.best_params_}")

    # Compare: OAT recommended params on test set
    oat_params = dict(default_params, **recommended, random_state=42)
    oat_model = model_class(**oat_params)
    oat_model.fit(x_train_search, y_train_search)
    oat_scores = oat_model.predict(x_test)
    oat_auc_mean, oat_auc_std = _lo_averaged_auc(y_test, oat_scores, test_lo)
    print(
        f"  OAT-recommended test AUC: {oat_auc_mean:.4f} +/- {oat_auc_std:.4f}"
    )

    # Best from randomized search on test set
    best_scores = search.best_estimator_.predict(x_test)
    best_auc_mean, best_auc_std = _lo_averaged_auc(y_test, best_scores,
                                                   test_lo)
    print(
        f"  RandomizedSearchCV test AUC: {best_auc_mean:.4f} +/- {best_auc_std:.4f}"
    )

    delta = best_auc_mean - oat_auc_mean
    print(f"  Delta: {delta:+.4f}")

    return search.best_params_, search.best_score_


def main():
    parser = argparse.ArgumentParser(
        description="Hyperparameter sweeps, randomized search, and "
        "feature-importance plots for the YAP-PSD regressors.")
    parser.add_argument(
        "--skip-randomized-search",
        action="store_true",
        help="Skip RandomizedSearchCV and use the hardcoded best params in "
        "SKIP_SEARCH_BEST_PARAMS. Use this when re-running just to "
        "regenerate plots after the search results are already known.")
    args = parser.parse_args()

    os.makedirs("plots", exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    data_cache = os.path.join(CACHE_DIR, "prepared_data.npz")
    wf_cache = os.path.join(CACHE_DIR, "avg_waveform.npy")
    gamma_wf_cache = os.path.join(CACHE_DIR, "avg_gamma_waveform.npy")

    cache_has_test_lo = (os.path.exists(data_cache)
                         and "test_lo" in np.load(data_cache).files)
    if os.path.exists(data_cache) and not cache_has_test_lo:
        print(f"  Cache {data_cache} predates test_lo; rebuilding.")
    if (os.path.exists(data_cache) and os.path.exists(wf_cache)
            and cache_has_test_lo):
        print("Loading cached prepared data...")
        data = np.load(data_cache)
        x_train, y_train = data["x_train"], data["y_train"]
        x_test, y_test = data["x_test"], data["y_test"]
        test_lo = data["test_lo"]
        avg_waveform = np.load(wf_cache)
        if os.path.exists(gamma_wf_cache):
            avg_gamma_waveform = np.load(gamma_wf_cache)
        else:
            f = ROOT.TFile.Open(ROOT_FILES_DIR + "Na22.root")
            avg_wf_graph = f.Get("average_waveform")
            avg_gamma_waveform = np.array([
                avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())
            ])
            f.Close()
            np.save(gamma_wf_cache, avg_gamma_waveform)
        print(f"  Train: {len(x_train)}, Test: {len(x_test)}")
    else:
        print("Loading alpha data (Am-241)...")
        alpha_features, alpha_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Am241.root",
            array_branch="Samples",
        )
        print(
            f"Alpha events: {len(alpha_features)}, waveform shape: {alpha_waveforms.shape}"
        )

        print("Loading gamma data (Na-22)...")
        gamma_features, gamma_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Na22.root",
            array_branch="Samples",
        )
        print(
            f"Gamma events: {len(gamma_features)}, waveform shape: {gamma_waveforms.shape}"
        )

        # Find the max n_training_samples across all configs so _prepare_data
        # allocates enough training data for the full sweep.
        max_train = DEFAULT_TRAIN_PER_CLASS
        for config in [RF_CONFIG, XGB_CONFIG]:
            for sweep in config["sweeps"]:
                if sweep["sweep_name"] == "n_training_samples":
                    max_train = max(max_train, max(sweep["values"]))

        print("Preparing train/test data")
        x_train, y_train, x_test, y_test, test_lo = _prepare_data(
            alpha_waveforms,
            gamma_waveforms,
            alpha_features,
            gamma_features,
            n_train_per_class=max_train,
            return_lo=True)

        f = ROOT.TFile.Open(ROOT_FILES_DIR + "Am241.root")
        avg_wf_graph = f.Get("average_waveform")
        avg_waveform = np.array(
            [avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())])
        f.Close()

        f = ROOT.TFile.Open(ROOT_FILES_DIR + "Na22.root")
        avg_wf_graph = f.Get("average_waveform")
        avg_gamma_waveform = np.array(
            [avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())])
        f.Close()

        np.savez(data_cache,
                 x_train=x_train,
                 y_train=y_train,
                 x_test=x_test,
                 y_test=y_test,
                 test_lo=test_lo)
        np.save(wf_cache, avg_waveform)
        np.save(gamma_wf_cache, avg_gamma_waveform)
        print(f"Prepared data cached to {CACHE_DIR}/")

    for config in [RF_CONFIG, XGB_CONFIG]:
        print(f"{config['name']} Hyperparameter Study")
        print(f"One-at-a-time sweeps")
        recommended = _run_all_sweeps(config, x_train, y_train, x_test,
                                      y_test, test_lo)
        optimized_params = dict(config["default_params"])
        optimized_params.update(recommended)
        print(f"{config['name']}: OAT recommended params: {recommended}")
        if args.skip_randomized_search:
            best_params = _best_params_from_regressors(config)
            print(f"Skipping randomized search; using params from "
                  f"regressors.py: {best_params}")
        else:
            print(f"RandomizedSearchCV around recommended values")
            best_params, _ = _run_randomized_search(config, recommended,
                                                    x_train, y_train, x_test,
                                                    y_test, test_lo)
        # Use the better params (from randomized search) for feature importance
        final_params = dict(config["default_params"])
        final_params.update(best_params)
        print(f"Feature importance (averaged over seeds)")
        _plot_shap_importance(config["model_class"], final_params,
                              config["prefix"], config["color"],
                              config["name"], x_train, y_train, avg_waveform,
                              config["explainer"])
        _plot_feature_importance(config["model_class"], final_params,
                                 config["prefix"], config["color"],
                                 config["name"], x_train, y_train,
                                 avg_waveform)

    for config in [GB_CONFIG]:
        print(f"{config['name']} Feature Importance")
        _plot_shap_importance(config["model_class"], config["default_params"],
                              config["prefix"], config["color"],
                              config["name"], x_train, y_train, avg_waveform,
                              config["explainer"])
        _plot_feature_importance(config["model_class"],
                                 config["default_params"], config["prefix"],
                                 config["color"], config["name"], x_train,
                                 y_train, avg_waveform)

    for config in [MLP_CONFIG]:
        print(f"{config['name']} SHAP Feature Importance")
        _plot_shap_importance(config["model_class"], config["default_params"],
                              config["prefix"], config["color"],
                              config["name"], x_train, y_train, avg_waveform,
                              config["explainer"])

    print("Done. All plots saved.")


if __name__ == "__main__":
    main()
