"""Shuffle study: train CNN on input-permuted waveforms.

Tests whether the CNN's convolutional inductive bias actually exploits
adjacency in the time samples. If shuffled-input AUC matches raw-input AUC
across multiple fixed permutations (with same permutation applied to train
and test), the conv layers aren't extracting extra signal from local
correlations -- which closes the rebuttal that trees/MLP miss performance
because they're order-blind.
"""
import os
import pickle
import numpy as np
import ROOT

from analysis_utilities.io import load_tree_data
from psd_utils import (ROOT_FILES_DIR, cumulative_shap,
                       style_residual_pad_axes, compute_mean_abs_shap,
                       per_lo_auc_from_scores, error_weighted_auc)
from torch_models import TorchCNNRegressor
import analysis_utilities

analysis_utilities.load_cpp_library()
ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPDF)

CACHE_DIR = "shuffle_study_cache"
SWEEP_CACHE_DIR = "sweep_cache"
DATA_CACHE = os.path.join(SWEEP_CACHE_DIR, "prepared_data.npz")
AVG_WF_CACHE = os.path.join(SWEEP_CACHE_DIR, "avg_waveform.npy")

N_PERMUTATIONS = 3
N_SEEDS_PER_PERM = 3
BASELINE_SEEDS = [42, 123, 256]
PERM_BASE_SEED = 1000
TRAIN_PER_CLASS = 10000

# KernelSHAP settings -- matched to parameter_study so the CNN curve is
# directly comparable to the MLP curve there.
SHAP_BG_SIZE = 100
SHAP_EXPLAIN_SIZE = 500


def _load_or_prepare_data():
    """Returns (x_train, y_train, x_test, y_test, test_lo).

    test_lo is a 1-D array aligned with x_test rows. Older caches that
    lack the test_lo array are rebuilt automatically so per-LO-bin AUC
    aggregation has the metadata it needs.
    """
    if os.path.exists(DATA_CACHE):
        d = np.load(DATA_CACHE)
        if "test_lo" in d.files:
            print(f"Loading prepared data from {DATA_CACHE}")
            return (d["x_train"], d["y_train"], d["x_test"], d["y_test"],
                    d["test_lo"])
        print(f"Cache {DATA_CACHE} predates test_lo; rebuilding.")
    else:
        print(f"{DATA_CACHE} not found, preparing from scratch")
    from parameter_study import _prepare_data
    print("Loading alpha data (Am-241)...")
    alpha_features, alpha_waveforms = load_tree_data(ROOT_FILES_DIR +
                                                     "Am241.root",
                                                     array_branch="Samples")
    print("Loading gamma data (Na-22)...")
    gamma_features, gamma_waveforms = load_tree_data(ROOT_FILES_DIR +
                                                     "Na22.root",
                                                     array_branch="Samples")
    return _prepare_data(alpha_waveforms,
                         gamma_waveforms,
                         alpha_features,
                         gamma_features,
                         n_train_per_class=TRAIN_PER_CLASS,
                         return_lo=True)


def _subsample_balanced(x, y, n_per_class, seed=42):
    rng = np.random.RandomState(seed)
    idx0 = np.where(y == 0)[0]
    idx1 = np.where(y == 1)[0]
    take0 = rng.choice(idx0, size=min(n_per_class, len(idx0)), replace=False)
    take1 = rng.choice(idx1, size=min(n_per_class, len(idx1)), replace=False)
    keep = np.concatenate([take0, take1])
    return x[keep], y[keep]


def _make_cnn(seed):
    return TorchCNNRegressor(conv_channels=(16, 32, 64),
                             kernel_size=5,
                             fc_sizes=(64, ),
                             dropout=0.2,
                             max_iter=50,
                             random_state=seed,
                             verbose=False)


def _lo_averaged_auc(y_test, scores, test_lo):
    """Bin scores by LO and return (mean, err) via inverse-variance weighting."""
    alpha_mask = (y_test == 0)
    gamma_mask = (y_test == 1)
    aucs, errs = per_lo_auc_from_scores(scores[alpha_mask], scores[gamma_mask],
                                        test_lo[alpha_mask],
                                        test_lo[gamma_mask])
    return error_weighted_auc(aucs, errs)


def _train_and_score(x_train, y_train, x_test, y_test, test_lo, seed):
    model = _make_cnn(seed)
    model.fit(x_train, y_train)
    auc, err = _lo_averaged_auc(y_test, model.predict(x_test), test_lo)
    return auc, err


def _compute_cnn_shap(model, x_train_view, seed):
    """KernelSHAP mean(|phi|) per *input column* of x_train_view.

    Caller is responsible for un-permuting back to original-feature order.
    """
    return compute_mean_abs_shap(model,
                                 x_train_view,
                                 seed,
                                 explainer="kernel",
                                 n_bg=SHAP_BG_SIZE,
                                 n_explain=SHAP_EXPLAIN_SIZE)


def _train_and_shap(x_train, y_train, x_test, y_test, test_lo, perm, seed):
    """Fit CNN on x_train[:, perm], score LO-averaged test AUC, compute SHAP,
    un-permute.

    Returns (auc, auc_err, shap_orig) where shap_orig is mean(|phi|) indexed
    by the *original* feature axis (i.e. shap_orig[k] is the importance of
    original sample k regardless of which permutation was used during
    training).
    """
    xs_train = x_train[:, perm]
    xs_test = x_test[:, perm]
    model = _make_cnn(seed)
    model.fit(xs_train, y_train)
    auc, err = _lo_averaged_auc(y_test, model.predict(xs_test), test_lo)
    shap_perm = _compute_cnn_shap(model, xs_train, seed)
    shap_orig = np.empty_like(shap_perm)
    shap_orig[perm] = shap_perm
    return auc, err, shap_orig


def _shap_cache_path(tag, seed):
    return os.path.join(CACHE_DIR, f"cnn_shap_{tag}_seed{seed}.npy")


def _summarize(label, aucs):
    a = np.array(aucs)
    return (f"  {label:>18}  mean={a.mean():.4f}  "
            f"std={a.std():.4f}  n={len(a)}")


def _plot_results(results, output_name):
    n_perms = len(results["shuffled_groups"])
    n_points = 1 + n_perms

    baseline = np.array(results["baseline"])
    b_mean = float(baseline.mean())
    b_std = float(baseline.std())

    means = [b_mean] + [float(np.mean(g)) for g in results["shuffled_groups"]]
    stds = [b_std] + [float(np.std(g)) for g in results["shuffled_groups"]]

    y_lo = min(means[i] - stds[i] for i in range(n_points))
    y_hi = max(means[i] + stds[i] for i in range(n_points))
    pad = max(0.005, 0.25 * (y_hi - y_lo))

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    frame = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", n_points,
                      -0.5, n_points - 0.5)
    frame.SetStats(0)
    frame.GetXaxis().SetBinLabel(1, "Raw")
    for i in range(n_perms):
        frame.GetXaxis().SetBinLabel(i + 2, f"Perm {i + 1}")
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitle("ROC AUC")
    frame.GetYaxis().SetRangeUser(y_lo - pad, y_hi + pad)
    frame.Draw()

    band = ROOT.TBox(-0.5, b_mean - b_std, n_points - 0.5, b_mean + b_std)
    band.SetFillColorAlpha(ROOT.kGray, 0.4)
    band.SetLineWidth(0)
    band.Draw("SAME")

    line = ROOT.TLine(-0.5, b_mean, n_points - 0.5, b_mean)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineStyle(2)
    line.Draw()

    x_arr = np.arange(n_points, dtype=np.float64)
    y_arr = np.array(means, dtype=np.float64)
    ex_arr = np.zeros(n_points, dtype=np.float64)
    ey_arr = np.array(stds, dtype=np.float64)

    graph = ROOT.TGraphErrors(n_points, x_arr, y_arr, ex_arr, ey_arr)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.4)
    graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph.Draw("PE SAME")

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _load_avg_waveform():
    if os.path.exists(AVG_WF_CACHE):
        return np.load(AVG_WF_CACHE)
    f = ROOT.TFile.Open(ROOT_FILES_DIR + "Am241.root")
    g = f.Get("average_waveform")
    wf = np.array([g.GetPointY(i) for i in range(g.GetN())])
    f.Close()
    return wf


def _plot_shap_comparison(baseline_shaps, shuffled_shaps, avg_waveform,
                          output_name):
    """Overlay raw-CNN SHAP vs shuffled-pooled CNN SHAP on the avg alpha
    waveform. Both SHAP arrays are in *original-feature* order (the shuffled
    set has already been un-permuted)."""
    baseline = np.array(baseline_shaps)
    shuffled = np.array(shuffled_shaps)

    base_mean = baseline.mean(axis=0)
    base_std = baseline.std(axis=0)
    shuf_mean = shuffled.mean(axis=0)
    shuf_std = shuffled.std(axis=0)

    wf_max = float(np.max(avg_waveform))
    avg_wf_norm = avg_waveform / wf_max if wf_max > 0 else avg_waveform

    # Normalize each curve to its own max so shape (not absolute scale) is
    # what's compared, matching parameter_study._plot_kernel_shap_importance.
    base_max = float(np.max(base_mean))
    shuf_max = float(np.max(shuf_mean))
    base_mean_n = base_mean / base_max if base_max > 0 else base_mean
    base_std_n = base_std / base_max if base_max > 0 else base_std
    shuf_mean_n = shuf_mean / shuf_max if shuf_max > 0 else shuf_mean
    shuf_std_n = shuf_std / shuf_max if shuf_max > 0 else shuf_std

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
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetYaxis().SetTitleOffset(1)
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(-0.1, 1.1)
    graph_waveform.Draw("AL")

    ex = np.zeros(n_points, dtype=np.float64)

    raw_color = ROOT.kAzure + 2
    graph_raw = ROOT.TGraphErrors(n_points, x_values,
                                  base_mean_n.astype(np.float64), ex,
                                  base_std_n.astype(np.float64))
    graph_raw.SetLineColor(raw_color)
    graph_raw.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_raw.SetFillColorAlpha(raw_color, 0.2)
    graph_raw.Draw("L3 SAME")

    shuf_color = ROOT.kOrange + 7
    graph_shuf = ROOT.TGraphErrors(n_points, x_values,
                                   shuf_mean_n.astype(np.float64), ex,
                                   shuf_std_n.astype(np.float64))
    graph_shuf.SetLineColor(shuf_color)
    graph_shuf.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_shuf.SetFillColorAlpha(shuf_color, 0.2)
    graph_shuf.Draw("L3 SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.45, 0.875, 0.6, 0.85)
    leg.AddEntry(graph_waveform, "Average #alpha Waveform", "l")
    leg.AddEntry(graph_raw, "CNN SHAP (raw)", "lf")
    leg.AddEntry(graph_shuf, "CNN SHAP (shuffled, un-permuted)", "lf")
    leg.SetMargin(0.1)
    leg.Draw()
    pad_top.SetTickx(0)

    pad_bot.cd()
    cum_raw = cumulative_shap(base_mean)
    g_cum_raw = ROOT.TGraph(n_points, x_values, cum_raw.astype(np.float64))
    g_cum_raw.SetLineColor(raw_color)
    g_cum_raw.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    style_residual_pad_axes(g_cum_raw, "Time [ns]")
    g_cum_raw.GetXaxis().SetLimits(0, x_values[-1])
    g_cum_raw.Draw("AL")

    cum_shuf = cumulative_shap(shuf_mean)
    g_cum_shuf = ROOT.TGraph(n_points, x_values, cum_shuf.astype(np.float64))
    g_cum_shuf.SetLineColor(shuf_color)
    g_cum_shuf.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    g_cum_shuf.Draw("L SAME")

    ROOT.PlottingUtils.SaveFigure(canvas, output_name, "",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def _save_latex_table(results, output_path):
    baseline = np.array(results["baseline"])
    pooled = np.concatenate([np.array(g) for g in results["shuffled_groups"]])

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{CNN test ROC AUC on raw vs. permuted "
                 r"waveform inputs. Each permutation is a fixed shuffle "
                 r"of the sample axis applied identically to train and "
                 r"test data. The `Raw' row gives the seed-only variance.}")
    lines.append(r"  \label{tab:shuffle}")
    lines.append(r"  \begin{tabular}{lccc}")
    lines.append(r"    \toprule")
    lines.append(r"    Input & Mean AUC & Std & Seeds \\")
    lines.append(r"    \midrule")
    lines.append(f"    Raw & {baseline.mean():.4f} & {baseline.std():.4f} & "
                 f"{len(baseline)} \\\\")
    for i, g in enumerate(results["shuffled_groups"]):
        arr = np.array(g)
        lines.append(f"    Perm {i + 1} & {arr.mean():.4f} & "
                     f"{arr.std():.4f} & {len(arr)} \\\\")
    lines.append(r"    \midrule")
    lines.append(f"    Shuffled (pooled) & {pooled.mean():.4f} & "
                 f"{pooled.std():.4f} & {len(pooled)} \\\\")
    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")
    table = "\n".join(lines)
    with open(output_path, "w") as f:
        f.write(table + "\n")
    print(f"LaTeX table written to {output_path}")
    print(table)


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs("plots", exist_ok=True)
    results_cache = os.path.join(CACHE_DIR, "results.pkl")

    if os.path.exists(results_cache):
        print(f"Loading cached results from {results_cache}")
        with open(results_cache, "rb") as fh:
            results = pickle.load(fh)
    else:
        x_train, y_train, x_test, y_test, test_lo = _load_or_prepare_data()
        if (y_train == 0).sum() > TRAIN_PER_CLASS:
            x_train, y_train = _subsample_balanced(x_train, y_train,
                                                   TRAIN_PER_CLASS)
        n_features = x_train.shape[1]
        print(f"x_train: {x_train.shape}  x_test: {x_test.shape}  "
              f"n_features: {n_features}")

        results = {"baseline": [], "shuffled_groups": []}

        print("Baseline (raw waveforms, varying seed):")
        for seed in BASELINE_SEEDS:
            auc, err = _train_and_score(x_train, y_train, x_test, y_test,
                                        test_lo, seed)
            print(f"  seed={seed}  AUC = {auc:.4f} +/- {err:.4f}")
            results["baseline"].append(auc)

        print("Shuffled inputs:")
        for p in range(N_PERMUTATIONS):
            perm_rng = np.random.RandomState(PERM_BASE_SEED + p)
            perm = perm_rng.permutation(n_features)
            xs_train = x_train[:, perm]
            xs_test = x_test[:, perm]
            group = []
            for s in range(N_SEEDS_PER_PERM):
                seed = BASELINE_SEEDS[s]
                auc, err = _train_and_score(xs_train, y_train, xs_test, y_test,
                                            test_lo, seed)
                print(f"  perm={p + 1} seed={seed}  "
                      f"AUC = {auc:.4f} +/- {err:.4f}")
                group.append(auc)
            results["shuffled_groups"].append(group)

        with open(results_cache, "wb") as fh:
            pickle.dump(results, fh)
        print(f"Results cached to {results_cache}")

    print()
    print(_summarize("Raw", results["baseline"]))
    for i, g in enumerate(results["shuffled_groups"]):
        print(_summarize(f"Perm {i + 1}", g))
    pooled = [a for g in results["shuffled_groups"] for a in g]
    print(_summarize("Shuffled pooled", pooled))

    _plot_results(results, "shuffle_study_auc")
    _save_latex_table(results, os.path.join(CACHE_DIR, "shuffle_table.txt"))

    _run_shap_pass()


def _run_shap_pass():
    """Compute KernelSHAP for baseline (raw) and each permutation, store
    everything in original-feature order, then plot the comparison."""
    print()
    print("=== SHAP pass ===")

    # Decide up front what's missing so we only load the dataset if we need to.
    baseline_paths = [(seed, _shap_cache_path("baseline", seed))
                      for seed in BASELINE_SEEDS]
    perm_paths = []
    for p in range(N_PERMUTATIONS):
        for s in range(N_SEEDS_PER_PERM):
            seed = BASELINE_SEEDS[s]
            perm_paths.append((p, seed, _shap_cache_path(f"perm{p}", seed)))

    need_baseline = [t for t in baseline_paths if not os.path.exists(t[1])]
    need_perm = [t for t in perm_paths if not os.path.exists(t[2])]

    if need_baseline or need_perm:
        x_train, y_train, x_test, y_test, test_lo = _load_or_prepare_data()
        if (y_train == 0).sum() > TRAIN_PER_CLASS:
            x_train, y_train = _subsample_balanced(x_train, y_train,
                                                   TRAIN_PER_CLASS)
        n_features = x_train.shape[1]
        identity = np.arange(n_features)

        for seed, path in need_baseline:
            print(f"  baseline SHAP seed={seed}...")
            _, _, shap_orig = _train_and_shap(x_train, y_train, x_test, y_test,
                                              test_lo, identity, seed)
            np.save(path, shap_orig)
            print(f"    saved {path}")

        # Recompute the same permutations the AUC loop used so they line up.
        perms = []
        for p in range(N_PERMUTATIONS):
            perm_rng = np.random.RandomState(PERM_BASE_SEED + p)
            perms.append(perm_rng.permutation(n_features))

        for p, seed, path in need_perm:
            print(f"  perm={p + 1} SHAP seed={seed}...")
            _, _, shap_orig = _train_and_shap(x_train, y_train, x_test, y_test,
                                              test_lo, perms[p], seed)
            np.save(path, shap_orig)
            print(f"    saved {path}")
    else:
        print("  all SHAP caches present, skipping compute")

    baseline_shaps = [np.load(path) for _, path in baseline_paths]
    shuffled_shaps = [np.load(path) for _, _, path in perm_paths]

    avg_waveform = _load_avg_waveform()
    _plot_shap_comparison(baseline_shaps, shuffled_shaps, avg_waveform,
                          "shuffle_study_shap")


if __name__ == "__main__":
    main()
