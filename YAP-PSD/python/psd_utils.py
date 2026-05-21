import os
import time
from pathlib import Path
import ROOT
import numpy as np
import pandas as pd
import pickle
from sklearn.metrics import roc_curve, auc
from joblib import Parallel, delayed, parallel_backend
from analysis_utilities.init_utils import set_root_preferences

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ROOT_FILES_DIR = str(PROJECT_ROOT / "root_files") + "/"

set_root_preferences(plots_dir=PROJECT_ROOT / "plots",
                     root_files_dir=PROJECT_ROOT / "root_files")

N_JOBS = 32
N_BOOTSTRAP = 250
ANALYSIS_CACHE_DIR = "analysis_cache"


def _bootstrap_chunk(y_true, scores, indices):
    """Compute AUC for a chunk of bootstrap resamples."""
    aucs = []
    for idx in indices:
        y_b = y_true[idx]
        s_b = scores[idx]
        if len(np.unique(y_b)) < 2:
            continue
        fpr_b, tpr_b, _ = roc_curve(y_b, s_b)
        aucs.append(auc(fpr_b, tpr_b))
    return aucs


def bootstrap_auc(y_true,
                  scores,
                  n_bootstrap=N_BOOTSTRAP,
                  random_state=42,
                  n_jobs=N_JOBS):
    """Estimate AUC and its uncertainty via bootstrap resampling.

    Returns
    -------
    auc_mean : float
    auc_std : float
    """
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    all_idx = rng.randint(0, n, size=(n_bootstrap, n))
    chunks = np.array_split(all_idx, n_jobs)
    with parallel_backend('threading', n_jobs=n_jobs):
        results = Parallel()(delayed(_bootstrap_chunk)(y_true, scores, chunk)
                             for chunk in chunks)
    aucs = np.array([a for chunk_aucs in results for a in chunk_aucs])
    return float(np.mean(aucs)), float(np.std(aucs))


def balance_test_set(alpha_features,
                     gamma_features,
                     random_state=42):
    """Subsample the majority class so |alpha| == |gamma|.

    AUC on a heavily imbalanced test set rewards majority-class accuracy;
    forcing balance gives a class-symmetric reading.
    """
    n_balanced = min(len(alpha_features), len(gamma_features))
    if len(alpha_features) > n_balanced:
        alpha_features = alpha_features.sample(
            n=n_balanced, random_state=random_state).reset_index(drop=True)
    if len(gamma_features) > n_balanced:
        gamma_features = gamma_features.sample(
            n=n_balanced, random_state=random_state).reset_index(drop=True)
    return alpha_features, gamma_features


def process_waveforms(waveforms, n_jobs=N_JOBS):
    """Normalize each waveform by its maximum value.

    Parameters
    ----------
    waveforms : numpy.ndarray
        2-D array of shape (n_waveforms, n_samples).
    n_jobs : int
        Unused, kept for interface compatibility.

    Returns
    -------
    numpy.ndarray
        Normalized waveforms.
    """
    maxvals = np.max(waveforms, axis=1, keepdims=True)
    maxvals[maxvals == 0] = 1.0
    return waveforms / maxvals


def column_name(regressor_name):
    """Convert a regressor display name to a DataFrame column name."""
    return regressor_name.replace(" ", "_") + "_Output"


def cumulative_shap(shap):
    """Cumulative SHAP, normalized so the final value is 1.

    Mean(|SHAP|) is non-negative so the cumulative is monotonically
    non-decreasing; the normalized curve runs from 0 to 1.
    """
    cum = np.cumsum(np.asarray(shap, dtype=np.float64))
    final = cum[-1]
    if final > 0:
        cum = cum / final
    return cum


def compute_mean_abs_shap(model,
                          X,
                          seed,
                          explainer,
                          n_bg=100,
                          n_explain=500):
    """Compute mean(|SHAP|) per input feature for a fitted model.

    explainer="tree" uses interventional TreeSHAP (exact for tree models:
    sklearn forests/GBMs, XGBoost, LightGBM). explainer="kernel" uses
    KernelSHAP (model-agnostic sampling). Both use a background sample of
    size n_bg drawn from X and explain n_explain rows from X, so tree and
    non-tree models can be compared under the same interventional Shapley
    formulation with the same background distribution.
    """
    import shap
    rng = np.random.RandomState(seed)
    bg_idx = rng.choice(len(X), size=min(n_bg, len(X)), replace=False)
    explain_idx = rng.choice(len(X),
                             size=min(n_explain, len(X)),
                             replace=False)
    if explainer == "tree":
        expl = shap.TreeExplainer(model,
                                  data=X[bg_idx],
                                  feature_perturbation="interventional")
        shap_vals = expl.shap_values(X[explain_idx], check_additivity=False)
    elif explainer == "kernel":
        expl = shap.KernelExplainer(model.predict, X[bg_idx])
        shap_vals = expl.shap_values(X[explain_idx])
    else:
        raise ValueError(f"Unknown explainer {explainer!r}")
    return np.mean(np.abs(shap_vals), axis=0)


def style_residual_pad_axes(graph,
                            x_title,
                            y_title="#splitline{Cumulative}{SHAP [a.u.]}"):
    """Apply residual-pad styling matching PlottingUtils::PlotFitWithResiduals.

    Mirrors the bottom-pad style of the C++ fit-with-residuals layout: larger
    NDC title/label sizes to compensate for the 30%-tall pad, dense y-axis
    title centered, fixed [0, 1] cumulative range with a small margin.

    Default y-title uses #splitline so that the long label fits in the
    short pad: under the 90 deg CCW rotation, the second argument lands
    closer to the y-axis (so "SHAP" reads next to the axis, "Cumulative"
    one column out).
    """
    graph.SetTitle("")
    graph.GetXaxis().SetTitle(x_title)
    graph.GetYaxis().SetTitle(y_title)
    graph.GetXaxis().SetTitleSize(0.13)
    graph.GetYaxis().SetTitleSize(0.13)
    graph.GetXaxis().SetLabelSize(0.12)
    graph.GetYaxis().SetLabelSize(0.12)
    graph.GetXaxis().SetTitleOffset(1.0)
    graph.GetYaxis().SetTitleOffset(0.4)
    graph.GetYaxis().SetNdivisions(505)
    graph.GetXaxis().SetNdivisions(510)
    graph.GetYaxis().CenterTitle(True)
    graph.GetYaxis().SetRangeUser(-0.05, 1.05)


def _train_or_load(regressor_cfg, x_train, y_train):
    """Train a regressor or load from cache. Returns (model, train_time_s or None)."""
    model_file = regressor_cfg["file"]
    name = regressor_cfg["name"]

    if os.path.exists(model_file):
        print(f"Loading existing model: {name} ({model_file})")
        with open(model_file, "rb") as f:
            return pickle.load(f), None

    print(f"Training new model: {name}...")
    model = regressor_cfg["model"]
    t0 = time.perf_counter()
    model.fit(x_train, y_train)
    train_time = time.perf_counter() - t0
    with open(model_file, "wb") as f:
        pickle.dump(model, f)
    return model, train_time


def _get_training_indices(alpha_features, gamma_features, random_state=42):
    """Reconstruct which event indices were selected for training.

    Mirrors the selection logic in regress_waveforms so that downstream code
    (e.g. proof_of_concept) can exclude training events deterministically.

    Returns
    -------
    alpha_train_idx, gamma_train_idx : numpy arrays of original-index positions
    """
    lo_lower = {"alpha": 375, "gamma": 0}
    lo_upper = {"alpha": 1575, "gamma": 1750}

    lo_mask_alpha = ((alpha_features["light_output"] <= lo_upper["alpha"])
                     & (alpha_features["light_output"] >= lo_lower["alpha"]))
    lo_mask_gamma = ((gamma_features["light_output"] <= lo_upper["gamma"])
                     & (gamma_features["light_output"] >= lo_lower["gamma"]))

    alpha_mask_idx = np.where(lo_mask_alpha.values)[0]
    gamma_mask_idx = np.where(lo_mask_gamma.values)[0]

    min_samples = 10000
    rng = np.random.RandomState(random_state)
    alpha_sel = rng.choice(len(alpha_mask_idx),
                           size=min_samples,
                           replace=False)
    gamma_sel = rng.choice(len(gamma_mask_idx),
                           size=min_samples,
                           replace=False)

    return alpha_mask_idx[alpha_sel], gamma_mask_idx[gamma_sel]


def _save_timing_table(timing_records, n_cores, cache_dir=None):
    """Save training and inference timing data to a LaTeX table."""
    if cache_dir is None:
        cache_dir = ANALYSIS_CACHE_DIR
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"    \centering")
    lines.append(r"    \begin{tabular}{lcc}")
    lines.append(r"        \hline \hline")
    lines.append(r"        \textbf{Method} & \textbf{Train [s]} & "
                 r"\textbf{Infer [$\mu$s/event]} \\")
    lines.append(r"        \hline")
    for r in timing_records:
        train = (f"{r['train_time']:.3f}"
                 if r.get("train_time") is not None else "--")
        infer = (f"{r['infer_per_event_us']:.1f}"
                 if r.get("infer_per_event_us") is not None else "--")
        lines.append(f"        {r['name']} & {train} & {infer} " + r"\\")
    lines.append(r"        \hline \hline")
    lines.append(r"    \end{tabular}")
    lines.append(r"    \caption{Training and inference time for each "
                 f"regressor ({n_cores} threads)" + r"}")
    lines.append(r"    \label{tab:timing}")
    lines.append(r"\end{table}")

    table_str = "\n".join(lines)
    output_path = os.path.join(cache_dir, "timing_table.txt")
    with open(output_path, "w") as f:
        f.write(table_str + "\n")
    print(f"\nLaTeX timing table saved to {output_path}")
    print(table_str)


HYPERPARAM_SPECS = [
    {
        "name": "Random Forest",
        "label": "tab:rf",
        "params": ["n_estimators", "max_depth", "max_samples", "max_features"],
    },
    {
        "name": "Gradient Boosting",
        "label": "tab:gb",
        "params": ["n_estimators", "max_depth", "learning_rate"],
    },
    {
        "name": "XGBoost",
        "label": "tab:xgb",
        "params": ["n_estimators", "max_depth", "learning_rate"],
    },
    {
        "name": "MLP",
        "label": "tab:mlp",
        "params": ["hidden_layer_sizes", "max_iter"],
    },
]

HYPERPARAM_PAIRS = [(0, 1), (2, 3)]


def _fmt_hyperparam_value(value):
    """Format a hyperparameter value for the LaTeX table."""
    if isinstance(value, str):
        return rf"\say{{{value}}}"
    if isinstance(value, tuple):
        return "(" + ", ".join(str(x) for x in value) + ")"
    return str(value)


def _build_hyperparam_subtable(name, label, rows, phantom_rows):
    """Build the LaTeX lines for a single hyperparameter subtable."""
    indent = "    "
    out = []
    out.append(rf"{indent}\begin{{subtable}}[t]{{0.45\textwidth}}")
    out.append(rf"{indent}    \centering")
    out.append(rf"{indent}    \begin{{tabular}}{{lc}}")
    out.append(rf"{indent}        \hline \hline")
    out.append(rf"{indent}        \textbf{{Parameter}} & \textbf{{Value}} \\")
    out.append(rf"{indent}        \hline")
    for param, value in rows:
        out.append(rf"{indent}        {param} & {value} \\")
    for param, value in phantom_rows:
        out.append(rf"{indent}        \phantom{{{param}}} & "
                   rf"\phantom{{{value}}} \\")
    out.append(rf"{indent}        \hline \hline")
    out.append(rf"{indent}    \end{{tabular}}")
    out.append(rf"{indent}    \caption{{{name}}}")
    out.append(rf"{indent}    \label{{{label}}}")
    out.append(rf"{indent}\end{{subtable}}")
    return out


def _save_hyperparameter_table(cache_dir=None):
    """Save a LaTeX 2x2 subtable showing each regressor's hyperparameters."""
    from regressors import get_default_regressors

    if cache_dir is None:
        cache_dir = ANALYSIS_CACHE_DIR

    models_by_name = {r["name"]: r["model"] for r in get_default_regressors()}

    rows_per_spec = []
    for spec in HYPERPARAM_SPECS:
        model = models_by_name[spec["name"]]
        rows = [(key.replace("_",
                             r"\_"), _fmt_hyperparam_value(getattr(model,
                                                                   key)))
                for key in spec["params"]]
        rows_per_spec.append(rows)

    phantoms_per_spec = [[] for _ in HYPERPARAM_SPECS]
    for i, j in HYPERPARAM_PAIRS:
        ni, nj = len(rows_per_spec[i]), len(rows_per_spec[j])
        if ni > nj:
            phantoms_per_spec[j] = rows_per_spec[i][nj:]
        elif nj > ni:
            phantoms_per_spec[i] = rows_per_spec[j][ni:]

    lines = []
    lines.append(r"\begin{table}[h!]")
    lines.append(r"    \centering")
    lines.append(r"    %")
    for pair_idx, (i, j) in enumerate(HYPERPARAM_PAIRS):
        if pair_idx > 0:
            lines.append(r"    %")
            lines.append(r"    \vspace{1em}")
            lines.append(r"    %")
        for k, spec_idx in enumerate((i, j)):
            spec = HYPERPARAM_SPECS[spec_idx]
            lines.extend(
                _build_hyperparam_subtable(spec["name"], spec["label"],
                                           rows_per_spec[spec_idx],
                                           phantoms_per_spec[spec_idx]))
            if k == 0:
                lines.append(r"    \hfill")
    lines.append(r"    %")
    lines.append(r"    \caption{Parameters used for each of the models}")
    lines.append(r"    \label{tab:ml}")
    lines.append(r"\end{table}")

    table_str = "\n".join(lines)
    output_path = os.path.join(cache_dir, "hyperparameter_table.txt")
    os.makedirs(cache_dir, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(table_str + "\n")
    print(f"\nLaTeX hyperparameter table saved to {output_path}")
    print(table_str)


def _run_plots(test_alpha_features,
               test_gamma_features,
               test_alpha_wf,
               test_gamma_wf,
               regressor_names,
               lo_lower_dict,
               lo_upper_dict,
               plot_prefix="",
               include_legacy_psd=True,
               auc_lo_subplot_label=None,
               auc_lo_y_range=(0.6, 1.01)):
    """Run all plotting and analysis on cached/computed test data."""
    alpha_lo_mask = (
        test_alpha_features["light_output"] <= lo_upper_dict["alpha"]) & (
            test_alpha_features["light_output"] >= lo_lower_dict["alpha"])

    gamma_lo_mask = (
        test_gamma_features["light_output"] <= lo_upper_dict["gamma"]) & (
            test_gamma_features["light_output"] >= lo_lower_dict["gamma"])

    test_alpha_features_filtered = test_alpha_features[alpha_lo_mask]
    test_gamma_features_filtered = test_gamma_features[gamma_lo_mask]

    if (plot_prefix == "" or plot_prefix == "full"):
        for name in regressor_names:
            col = column_name(name)
            safe_name = name.replace(" ", "_").lower()
            _plot_score_histogram(
                test_alpha_features_filtered[col].values,
                test_gamma_features_filtered[col].values,
                f"Test Set Scores ({name})",
                f"{plot_prefix}test_score_histogram_{safe_name}",
                regressor_name=name,
            )

        alpha_lo_mask_900_1200 = (
            test_alpha_features["light_output"]
            <= 1200) & (test_alpha_features["light_output"] >= 900)

        gamma_lo_mask_900_1200 = (
            test_gamma_features["light_output"]
            <= 1200) & (test_gamma_features["light_output"] >= 900)

        test_alpha_features_900_1200 = test_alpha_features[
            alpha_lo_mask_900_1200]
        test_gamma_features_900_1200 = test_gamma_features[
            gamma_lo_mask_900_1200]

        for name in regressor_names:
            col = column_name(name)
            safe_name = name.replace(" ", "_").lower()
            _plot_score_histogram(
                test_alpha_features_900_1200[col].values,
                test_gamma_features_900_1200[col].values,
                f"Test Set Scores ({name})",
                f"{plot_prefix}test_score_histogram_900_1200_{safe_name}",
                regressor_name=name,
            )

    _analyze_all_methods(test_alpha_features_filtered,
                         test_gamma_features_filtered,
                         regressor_names,
                         plot_prefix=plot_prefix,
                         include_legacy_psd=include_legacy_psd)

    _plot_auc_vs_light_output(test_alpha_features,
                              test_gamma_features,
                              regressor_names,
                              plot_prefix=plot_prefix,
                              include_legacy_psd=include_legacy_psd,
                              subplot_label=auc_lo_subplot_label,
                              y_range=auc_lo_y_range)

    return (
        (test_alpha_wf, test_gamma_wf),
        (test_alpha_features, test_gamma_features),
    )


def regress_waveforms(waveforms,
                      features,
                      regressors,
                      process_func=process_waveforms,
                      random_state=42,
                      cache_dir=None,
                      plot_prefix="",
                      include_legacy_psd=True,
                      auc_lo_subplot_label=None,
                      auc_lo_y_range=(0.6, 1.01),
                      skip_plots=False):
    """Train and evaluate ML models using waveforms from ROOT files.

    Parameters
    ----------
    waveforms : tuple of (alpha_waveforms_ndarray, gamma_waveforms_ndarray)
    features : tuple of (alpha_features_df, gamma_features_df)
    regressors : list of dict
        Each dict has keys "name", "model" (unfitted sklearn regressor),
        and "file" (path for caching the trained model).
    cache_dir : str or None
        Directory for caching analysis results. Defaults to ANALYSIS_CACHE_DIR.
    plot_prefix : str
        Prefix prepended to all plot output filenames.
    """
    if cache_dir is None:
        cache_dir = ANALYSIS_CACHE_DIR

    lo_lower_dict = {"alpha": 375, "gamma": 0}
    lo_upper_dict = {"alpha": 1575, "gamma": 1750}

    # Check for cached analysis results
    cache_files = {
        "alpha_feat": os.path.join(cache_dir, "test_alpha_features.pkl"),
        "gamma_feat": os.path.join(cache_dir, "test_gamma_features.pkl"),
        "waveforms": os.path.join(cache_dir, "test_waveforms.npz"),
        "names": os.path.join(cache_dir, "regressor_names.pkl"),
    }
    if all(os.path.exists(f) for f in cache_files.values()):
        print("Loading cached analysis results...")
        test_alpha_features = pd.read_pickle(cache_files["alpha_feat"])
        test_gamma_features = pd.read_pickle(cache_files["gamma_feat"])
        wf_data = np.load(cache_files["waveforms"])
        test_alpha_wf = wf_data["alpha"]
        test_gamma_wf = wf_data["gamma"]
        with open(cache_files["names"], "rb") as f:
            regressor_names = pickle.load(f)
        print(f"  Loaded {len(regressor_names)} regressors: "
              f"{', '.join(regressor_names)}")
        if skip_plots:
            return ((test_alpha_wf, test_gamma_wf),
                    (test_alpha_features, test_gamma_features))
        # Skip to plotting
        return _run_plots(test_alpha_features,
                          test_gamma_features,
                          test_alpha_wf,
                          test_gamma_wf,
                          regressor_names,
                          lo_lower_dict,
                          lo_upper_dict,
                          plot_prefix=plot_prefix,
                          include_legacy_psd=include_legacy_psd,
                          auc_lo_subplot_label=auc_lo_subplot_label,
                          auc_lo_y_range=auc_lo_y_range)

    os.makedirs(cache_dir, exist_ok=True)

    alpha_waveforms, gamma_waveforms = waveforms
    alpha_features, gamma_features = features

    alpha_train_original_indices, gamma_train_original_indices = \
        _get_training_indices(alpha_features, gamma_features, random_state)

    train_alpha_wf = alpha_waveforms[alpha_train_original_indices]
    train_gamma_wf = gamma_waveforms[gamma_train_original_indices]

    min_samples = len(alpha_train_original_indices)
    print(f"Samples for balanced training: {min_samples}")

    # Process training waveforms
    with parallel_backend('threading', n_jobs=2):
        train_results = Parallel()(
            delayed(process_func)(group)
            for group in [train_alpha_wf, train_gamma_wf])

    _plot_sample_waveforms(train_results,
                           n_samples=5,
                           random_state=random_state,
                           plot_prefix=plot_prefix)

    x_train = np.vstack(train_results)
    y_train = np.array([0] * len(train_results[0]) +
                       [1] * len(train_results[1]))

    # Train or load each regressor
    n_cores = 32
    trained_models = {}
    timing_records = []
    for cfg in regressors:
        model, train_time = _train_or_load(cfg, x_train, y_train)
        trained_models[cfg["name"]] = model
        if train_time is not None:
            print(f"  Training time: {train_time:.3f} s "
                  f"({train_time * n_cores:.3f} s on single-core)")
            timing_records.append({
                "name": cfg["name"],
                "train_time": train_time,
                "train_time_per_core": train_time * n_cores,
            })

    # Create test data by dropping training samples
    all_alpha_idx = np.arange(len(alpha_waveforms))
    test_alpha_idx = np.setdiff1d(all_alpha_idx, alpha_train_original_indices)
    test_alpha_wf = alpha_waveforms[test_alpha_idx]
    test_alpha_features = alpha_features.iloc[test_alpha_idx].reset_index(
        drop=True)

    all_gamma_idx = np.arange(len(gamma_waveforms))
    test_gamma_idx = np.setdiff1d(all_gamma_idx, gamma_train_original_indices)
    test_gamma_wf = gamma_waveforms[test_gamma_idx]
    test_gamma_features = gamma_features.iloc[test_gamma_idx].reset_index(
        drop=True)

    # Process test waveforms for prediction
    with parallel_backend('threading', n_jobs=2):
        test_results = Parallel()(delayed(process_func)(group)
                                  for group in [test_alpha_wf, test_gamma_wf])

    X_test = np.vstack(test_results)

    # Predict with each regressor and add output columns
    test_alpha_features = test_alpha_features.copy()
    test_gamma_features = test_gamma_features.copy()

    n_test = len(X_test)
    regressor_names = []
    for name, model in trained_models.items():
        col = column_name(name)
        t0 = time.perf_counter()
        y_pred = model.predict(X_test)
        infer_time = time.perf_counter() - t0
        infer_per_event = infer_time / n_test
        print(f"{name}: inference {infer_time:.3f} s total, "
              f"{infer_per_event * 1e6:.1f} us/event, "
              f"{infer_per_event * n_cores * 1e6:.1f} us/event on single core")
        # Update timing record or create one (if model was loaded from cache)
        record = next((r for r in timing_records if r["name"] == name), None)
        if record is None:
            timing_records.append({
                "name":
                name,
                "train_time":
                None,
                "train_time_per_core":
                None,
                "infer_time":
                infer_time,
                "infer_per_event_us":
                infer_per_event * 1e6,
                "infer_per_event_per_core_us":
                infer_per_event * n_cores * 1e6,
            })
        else:
            record["infer_time"] = infer_time
            record["infer_per_event_us"] = infer_per_event * 1e6
            record[
                "infer_per_event_per_core_us"] = infer_per_event * n_cores * 1e6
        test_alpha_features[col] = y_pred[:len(test_alpha_wf)]
        test_gamma_features[col] = y_pred[len(test_alpha_wf):]
        regressor_names.append(name)

    # Only save timing table if we actually trained (otherwise we'd overwrite
    # the table with missing training times)
    if any(r.get("train_time") is not None for r in timing_records):
        _save_timing_table(timing_records, n_cores, cache_dir=cache_dir)

    _save_hyperparameter_table(cache_dir=cache_dir)

    # Save analysis cache
    test_alpha_features.to_pickle(cache_files["alpha_feat"])
    test_gamma_features.to_pickle(cache_files["gamma_feat"])
    np.savez(cache_files["waveforms"],
             alpha=test_alpha_wf,
             gamma=test_gamma_wf)
    with open(cache_files["names"], "wb") as f:
        pickle.dump(regressor_names, f)
    print(f"Analysis cache saved to {cache_dir}/")

    if skip_plots:
        return ((test_alpha_wf, test_gamma_wf),
                (test_alpha_features, test_gamma_features))
    return _run_plots(test_alpha_features,
                      test_gamma_features,
                      test_alpha_wf,
                      test_gamma_wf,
                      regressor_names,
                      lo_lower_dict,
                      lo_upper_dict,
                      plot_prefix=plot_prefix,
                      include_legacy_psd=include_legacy_psd,
                      auc_lo_subplot_label=auc_lo_subplot_label,
                      auc_lo_y_range=auc_lo_y_range)


def _analyze_all_methods(test_alpha_features,
                         test_gamma_features,
                         regressor_names,
                         plot_prefix="",
                         include_legacy_psd=True):
    """ROC analysis comparing ML regressors, Charge Comparison, and Shape Indicator PSD.

    The test set is class-balanced (majority class subsampled) so AUC and the
    5%-FPR thresholds are not inflated by the gamma/alpha population imbalance.
    Legacy PSD methods (CC/SI) operate on the full waveform and are not
    meaningful comparators for waveform-slice studies; set
    include_legacy_psd=False to omit them.
    """
    test_alpha_features, test_gamma_features = balance_test_set(
        test_alpha_features, test_gamma_features)
    print(f"  Balanced ROC test set: {len(test_alpha_features)} alpha, "
          f"{len(test_gamma_features)} gamma")

    test_features = pd.concat([test_alpha_features,
                               test_gamma_features]).reset_index(drop=True)

    y_true = np.array([0] * len(test_alpha_features) +
                      [1] * len(test_gamma_features))

    all_methods = [column_name(n) for n in regressor_names]
    all_method_names = list(regressor_names)

    if include_legacy_psd:
        all_methods += ["charge_comparison", "clean_shape_indicator"]
        all_method_names += ["Charge Comparison", "Shape Indicator"]

    _plot_roc_curves(test_features, y_true, all_methods, all_method_names,
                     f"{plot_prefix}roc_curves")


def compute_lo_bin_edges(alpha_lo,
                         gamma_lo,
                         lo_min=375,
                         lo_max=1575,
                         n_bins=14):
    """Return LO bin edges that target equal minority-class counts per bin."""
    fine_edges = np.arange(lo_min, lo_max + 25, 25)
    alpha_hist, _ = np.histogram(alpha_lo, bins=fine_edges)
    gamma_hist, _ = np.histogram(gamma_lo, bins=fine_edges)
    minority_hist = np.minimum(alpha_hist, gamma_hist)
    cum_minority = np.cumsum(minority_hist)
    total_minority = cum_minority[-1]
    target_per_bin = total_minority / n_bins

    bin_edges = [lo_min]
    running = 0
    for i in range(len(minority_hist)):
        running += minority_hist[i]
        if running >= target_per_bin and len(bin_edges) < n_bins:
            bin_edges.append(fine_edges[i + 1])
            running = 0
    bin_edges.append(lo_max)
    return bin_edges


def per_lo_auc_from_scores(alpha_scores,
                           gamma_scores,
                           alpha_lo,
                           gamma_lo,
                           random_state=42,
                           min_per_bin=5):
    """Compute balanced per-LO-bin AUC + bootstrap error from raw scores.

    Bins are chosen to equalize minority-class counts per bin. Within each
    bin the majority class is subsampled to the minority size so AUC is
    insensitive to class imbalance.

    Returns (aucs, errs) — lists aligned with the LO bins. Bins with too
    few events get NaN.
    """
    alpha_scores = np.asarray(alpha_scores)
    gamma_scores = np.asarray(gamma_scores)
    alpha_lo = np.asarray(alpha_lo)
    gamma_lo = np.asarray(gamma_lo)

    bin_edges = compute_lo_bin_edges(alpha_lo, gamma_lo)
    aucs = []
    errs = []
    rng = np.random.RandomState(random_state)
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        a_mask = (alpha_lo >= lo) & (alpha_lo < hi)
        g_mask = (gamma_lo >= lo) & (gamma_lo < hi)
        a_s = alpha_scores[a_mask]
        g_s = gamma_scores[g_mask]
        n_balanced = min(len(a_s), len(g_s))
        if n_balanced < min_per_bin:
            aucs.append(float("nan"))
            errs.append(float("nan"))
            continue
        if len(a_s) > n_balanced:
            a_s = a_s[rng.choice(len(a_s), n_balanced, replace=False)]
        if len(g_s) > n_balanced:
            g_s = g_s[rng.choice(len(g_s), n_balanced, replace=False)]
        y_true = np.array([0] * n_balanced + [1] * n_balanced)
        scores = np.concatenate([a_s, g_s])
        m, s = bootstrap_auc(y_true, scores)
        aucs.append(m)
        errs.append(s)
    return aucs, errs


def _compute_per_lo_auc(test_alpha_features,
                        test_gamma_features,
                        regressor_names,
                        include_legacy_psd=True):
    """Compute per-LO-bin balanced AUC for every method.

    Returns
    -------
    bin_centers : list of float
    lo_bins : list of (lo, hi)
    auc_results : dict[method_name -> list of float]
    auc_errors : dict[method_name -> list of float]
    all_method_names : list of str (display names)
    """
    alpha_lo = test_alpha_features["light_output"].values
    gamma_lo = test_gamma_features["light_output"].values
    bin_edges = compute_lo_bin_edges(alpha_lo, gamma_lo)

    lo_bins = [(bin_edges[i], bin_edges[i + 1])
               for i in range(len(bin_edges) - 1)]
    bin_centers = [0.5 * (lo + hi) for lo, hi in lo_bins]

    all_methods = [column_name(n) for n in regressor_names]
    all_method_names = list(regressor_names)
    if include_legacy_psd:
        all_methods += ["charge_comparison", "clean_shape_indicator"]
        all_method_names += ["Charge Comparison", "Shape Indicator"]

    auc_results = {name: [] for name in all_method_names}
    auc_errors = {name: [] for name in all_method_names}

    for lo, hi in lo_bins:
        alpha_mask = (test_alpha_features["light_output"]
                      >= lo) & (test_alpha_features["light_output"] < hi)
        gamma_mask = (test_gamma_features["light_output"]
                      >= lo) & (test_gamma_features["light_output"] < hi)

        alpha_bin = test_alpha_features[alpha_mask]
        gamma_bin = test_gamma_features[gamma_mask]

        n_alpha = len(alpha_bin)
        n_gamma = len(gamma_bin)
        n_balanced = min(n_alpha, n_gamma)
        print(f"Light output bin [{lo}, {hi}) keVee: {n_alpha} alpha, "
              f"{n_gamma} gamma  -> balanced to {n_balanced} each")

        if n_balanced < 5:
            for name in all_method_names:
                auc_results[name].append(np.nan)
                auc_errors[name].append(np.nan)
            continue

        if n_alpha > n_balanced:
            alpha_bin = alpha_bin.sample(
                n=n_balanced, random_state=42).reset_index(drop=True)
        if n_gamma > n_balanced:
            gamma_bin = gamma_bin.sample(
                n=n_balanced, random_state=42).reset_index(drop=True)

        combined = pd.concat([alpha_bin, gamma_bin]).reset_index(drop=True)
        y_true = np.array([0] * n_balanced + [1] * n_balanced)

        for method, name in zip(all_methods, all_method_names):
            scores = combined[method].values
            if method in ("raw_shape_indicator", "clean_shape_indicator"):
                scores = -scores
            auc_mean, auc_std = bootstrap_auc(y_true, scores)
            auc_results[name].append(auc_mean)
            auc_errors[name].append(auc_std)

    return bin_centers, lo_bins, auc_results, auc_errors, all_method_names


def error_weighted_auc(aucs, errs):
    """Inverse-variance-weighted mean and combined uncertainty.

    NaNs (from bins with insufficient statistics) are dropped. Returns
    (nan, nan) if no valid bins remain.
    """
    aucs = np.asarray(aucs, dtype=np.float64)
    errs = np.asarray(errs, dtype=np.float64)
    valid = ~(np.isnan(aucs) | np.isnan(errs) | (errs <= 0))
    if not np.any(valid):
        return float("nan"), float("nan")
    w = 1.0 / (errs[valid] ** 2)
    mean = float(np.sum(w * aucs[valid]) / np.sum(w))
    err = float(1.0 / np.sqrt(np.sum(w)))
    return mean, err


def _plot_auc_vs_light_output(test_alpha_features,
                              test_gamma_features,
                              regressor_names,
                              plot_prefix="",
                              include_legacy_psd=True,
                              subplot_label=None,
                              y_range=(0.6, 1.01)):
    """Plot ROC AUC as a function of light output for all classifiers."""
    bin_centers, lo_bins, auc_results, auc_errors, all_method_names = (
        _compute_per_lo_auc(test_alpha_features,
                            test_gamma_features,
                            regressor_names,
                            include_legacy_psd=include_legacy_psd))

    # Plot
    colors = list(ROOT.PlottingUtils.GetDefaultColors())
    canvas = ROOT.TCanvas("auc_lo", "", 1200, 600)

    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.72, 1.0)
    pad_plot.SetRightMargin(0.02)
    pad_plot.Draw()

    pad_leg = ROOT.TPad("pad_leg", "", 0.72, 0.0, 1.0, 1.0)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_leg.Draw()

    pad_plot.cd()

    x_arr = np.array(bin_centers, dtype=np.float64)
    graphs = []
    graph_names = []

    ex_arr = np.array([0.5 * (hi - lo) for lo, hi in lo_bins],
                      dtype=np.float64)

    for i, name in enumerate(all_method_names):
        y_arr = np.array(auc_results[name], dtype=np.float64)
        ey_arr = np.array(auc_errors[name], dtype=np.float64)
        valid = ~np.isnan(y_arr)
        if not np.any(valid):
            continue

        graph = ROOT.TGraphErrors(int(np.sum(valid)), x_arr[valid],
                                  y_arr[valid], ex_arr[valid], ey_arr[valid])
        graph.SetLineColor(colors[i])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        graph.SetMarkerColor(colors[i])
        graph.SetMarkerStyle(20 + i)
        graph.SetMarkerSize(1.2)
        graphs.append(graph)
        graph_names.append(name)

        if len(graphs) == 1:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Light Output [keVee]")
            graph.GetYaxis().SetTitle("ROC AUC")
            graph.GetYaxis().SetTitleOffset(1)
            graph.GetXaxis().SetRangeUser(300, 2000)
            graph.GetYaxis().SetRangeUser(*y_range)
            graph.Draw("APE")
        else:
            graph.Draw("P SAME")

    if subplot_label is not None:
        _subplot_label_text = ROOT.PlottingUtils.AddText(
            subplot_label, 0.92, 0.84)

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.0, 0.95, 0.2, 0.85)
    leg.SetMargin(0.15)
    for graph, name in zip(graphs, graph_names):
        leg.AddEntry(graph, name, "lp")

    leg.Draw()
    ROOT.PlottingUtils.SaveFigure(canvas, f"{plot_prefix}auc_vs_light_output",
                                  "", ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()

    print(
        f"AUC vs light output plot saved to {plot_prefix}auc_vs_light_output")


def _plot_score_histogram(alpha_scores,
                          gamma_scores,
                          title,
                          output_path,
                          regressor_name="Regressor"):
    """Plot score histogram using ROOT"""
    canvas = ROOT.PlottingUtils.GetConfiguredCanvas(ROOT.kTRUE)

    all_scores = np.concatenate([alpha_scores, gamma_scores])
    score_min = np.min(all_scores)
    score_max = np.max(all_scores)

    h_alpha = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", 125,
                        score_min, score_max)
    h_gamma = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", 125,
                        score_min, score_max)

    for val in alpha_scores:
        h_alpha.Fill(val)
    for val in gamma_scores:
        h_gamma.Fill(val)

    ROOT.PlottingUtils.ConfigureHistogram(h_alpha, ROOT.kRed + 2)
    ROOT.PlottingUtils.ConfigureHistogram(h_gamma, ROOT.kBlue + 2)

    h_alpha.GetXaxis().SetTitle(f"{regressor_name} Score")
    h_alpha.GetYaxis().SetTitle("Counts")
    h_alpha.SetTitle("")

    max_val = max(h_alpha.GetMaximum(), h_gamma.GetMaximum())
    h_alpha.SetMaximum(max_val * 1.2)
    h_alpha.Draw("HIST")
    h_gamma.Draw("HIST SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.4, 0.6, 0.7, 0.85)
    leg.AddEntry(h_alpha, f"Am-241 (#alpha)", "f")
    leg.AddEntry(h_gamma, f"Na-22 (#gamma)", "f")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(
        canvas,
        output_path,
        "",
        ROOT.PlotSaveOptions.kLOG,
    )
    canvas.Close()
    h_alpha.Delete()
    h_gamma.Delete()


def _plot_roc_curves(test_features, y_true, methods, method_names,
                     output_name):
    """Plot ROC curves for all methods using ROOT"""
    colors = list(ROOT.PlottingUtils.GetDefaultColors())

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    roc_graphs = []
    leg = ROOT.PlottingUtils.AddLegend(0.42, 0.84, 0.2, 0.6)
    leg.SetMargin(0.1)

    for i, (method, name) in enumerate(zip(methods, method_names)):
        scores = test_features[method].values

        if method == "raw_shape_indicator" or method == "clean_shape_indicator":
            scores_to_use = -scores
        else:
            scores_to_use = scores

        fpr, tpr, thresholds = roc_curve(y_true, scores_to_use)

        index = np.argmin(np.abs(fpr - 0.05))
        threshold_at_5pct_fpr = thresholds[index]
        tpr_at_5pct_fpr = tpr[index]
        actual_fpr = fpr[index]

        print(f"{name}:")
        if method == "raw_shape_indicator" or method == "clean_shape_indicator":
            original_threshold = -threshold_at_5pct_fpr
            print(
                f"  Inverted threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}"
            )
            print(
                f"  Original threshold (lower is better): {original_threshold:.6f}"
            )
        else:
            print(
                f"  Threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}"
            )
        print(f"  TPR at {actual_fpr:.3f} FPR: {tpr_at_5pct_fpr:.3f}")

        auc_mean, auc_std = bootstrap_auc(y_true, scores_to_use)

        roc_graph = ROOT.TGraph(len(fpr), fpr, tpr)
        roc_graph.SetLineColor(colors[i])
        roc_graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        roc_graphs.append(roc_graph)

        if i == 0:
            roc_graph.SetTitle("")
            roc_graph.GetXaxis().SetTitle(
                "False Positive Rate (1 - Specificity)")
            roc_graph.GetYaxis().SetTitle("True Positive Rate (Sensitivity)")
            roc_graph.GetXaxis().SetRangeUser(0, 1)
            roc_graph.GetYaxis().SetRangeUser(0, 1)
            roc_graph.Draw("AL")
        else:
            roc_graph.Draw("L SAME")

        leg.AddEntry(roc_graph, f"{name} (AUC = {auc_mean:.3f})", "l")

    diagonal = ROOT.TLine(0, 0, 1, 1)
    diagonal.SetLineColor(ROOT.kBlack)
    diagonal.SetLineStyle(2)
    diagonal.Draw()

    leg.Draw()
    ROOT.PlottingUtils.SaveFigure(
        canvas,
        output_name,
        "",
        ROOT.PlotSaveOptions.kLINEAR,
    )


def _plot_sample_waveforms(waveforms_tuple,
                           n_samples=5,
                           random_state=42,
                           plot_prefix=""):
    """Plot sample waveforms for alpha and gamma classes after normalization.

    Parameters
    ----------
    waveforms_tuple : tuple of (alpha_waveforms, gamma_waveforms)
        Raw waveform arrays before normalization.
    n_samples : int
        Number of sample waveforms to plot per class.
    random_state : int
        Random seed for reproducibility.
    """
    alpha_waveforms, gamma_waveforms = waveforms_tuple

    rng = np.random.RandomState(random_state)
    alpha_indices = rng.choice(len(alpha_waveforms),
                               size=min(n_samples, len(alpha_waveforms)),
                               replace=False)
    gamma_indices = rng.choice(len(gamma_waveforms),
                               size=min(n_samples, len(gamma_waveforms)),
                               replace=False)

    alpha_samples = alpha_waveforms[alpha_indices]
    gamma_samples = gamma_waveforms[gamma_indices]

    n_points = alpha_samples.shape[1]
    x_values = np.arange(n_points) * 2

    canvas_alpha = ROOT.PlottingUtils.GetConfiguredCanvas()

    graphs_alpha = []
    leg_alpha = ROOT.PlottingUtils.AddLegend(0.65, 0.88, 0.65, 0.85)

    colors = [
        ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kOrange + 2,
        ROOT.kMagenta + 2
    ]

    for i, waveform in enumerate(alpha_samples):
        graph = ROOT.TGraph(n_points, x_values.astype(np.float64),
                            waveform.astype(np.float64))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

        if i == 0:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Time [ns]")
            graph.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
            graph.GetXaxis().SetRangeUser(0, x_values[-1])
            graph.GetYaxis().SetRangeUser(-0.1, 1.1)
            graph.Draw("AL")
        else:
            graph.Draw("L SAME")

        graphs_alpha.append(graph)
        leg_alpha.AddEntry(graph, f"Am-241 Waveform {i+1}", "l")

    text = ROOT.PlottingUtils.AddText("(a) Am-241", 0.78, 0.78)
    text.SetTextSize(35)

    ROOT.PlottingUtils.SaveFigure(canvas_alpha,
                                  f"{plot_prefix}sample_waveforms_alpha", "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas_alpha.Close()

    canvas_gamma = ROOT.PlottingUtils.GetConfiguredCanvas()

    graphs_gamma = []
    leg_gamma = ROOT.PlottingUtils.AddLegend(0.65, 0.88, 0.65, 0.85)

    for i, waveform in enumerate(gamma_samples):
        graph = ROOT.TGraph(n_points, x_values.astype(np.float64),
                            waveform.astype(np.float64))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

        if i == 0:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Time [ns]")
            graph.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
            graph.GetXaxis().SetRangeUser(0, x_values[-1])
            graph.GetYaxis().SetRangeUser(-0.1, 1.1)
            graph.Draw("AL")
        else:
            graph.Draw("L SAME")

        graphs_gamma.append(graph)
        leg_gamma.AddEntry(graph, f"Na-22 Waveform {i+1}", "l")

    text = ROOT.PlottingUtils.AddText("(b) Na-22", 0.78, 0.78)
    text.SetTextSize(35)

    ROOT.PlottingUtils.SaveFigure(canvas_gamma,
                                  f"{plot_prefix}sample_waveforms_gamma", "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas_gamma.Close()

    print(f"Sample waveform plots saved: {plot_prefix}sample_waveforms_alpha, "
          f"{plot_prefix}sample_waveforms_gamma")
