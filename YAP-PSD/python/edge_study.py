import numpy as np
import os
import pandas as pd
from analysis_utilities.io import load_tree_data
from psd_utils import (regress_waveforms, process_waveforms,
                       _compute_per_lo_auc, error_weighted_auc, column_name,
                       ANALYSIS_CACHE_DIR, ROOT_FILES_DIR)
from regressors import get_default_regressors
import analysis_utilities

analysis_utilities.load_cpp_library()
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPDF)

EDGE_STUDY_CACHE_DIR = "edge_study_cache"

SWEEP_CACHE_DIR = "sweep_cache"


def _load_peak_sample():
    """Find the peak sample index from the average alpha waveform."""
    wf_cache = os.path.join(SWEEP_CACHE_DIR, "avg_waveform.npy")
    if os.path.exists(wf_cache):
        avg_waveform = np.load(wf_cache)
    else:
        f = ROOT.TFile.Open(ROOT_FILES_DIR + "Am241.root")
        avg_wf_graph = f.Get("average_waveform")
        avg_waveform = np.array(
            [avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())])
        f.Close()
    peak = int(np.argmax(avg_waveform))
    print(f"Peak sample from average waveform: {peak}")
    return peak


FIRST_20NS_SAMPLES = 10  # 2 ns/sample digitizer -> 20 ns = 10 samples


def _make_edge_processors(peak_sample):
    """Return rising/falling/first-20ns edge processing functions.

    Each truncates first, then normalizes by the max of the *truncated*
    window. Normalizing the full waveform and then slicing leaves pre-peak
    samples divided by the full-pulse peak, which leaks light-output (and
    therefore class) information into windows that should carry no shape.
    """

    def process_rising_edge(waveforms):
        """Keep the rising edge (up to peak), then normalize."""
        return process_waveforms(waveforms[:, :peak_sample + 1])

    def process_falling_edge(waveforms):
        """Keep the falling edge (from peak), then normalize."""
        return process_waveforms(waveforms[:, peak_sample:])

    def process_first_20ns(waveforms):
        """Keep the first 20 ns, then normalize."""
        return process_waveforms(waveforms[:, :FIRST_20NS_SAMPLES])

    return process_rising_edge, process_falling_edge, process_first_20ns


def _variant_cache_dir(variant_name):
    """Cache directory for a variant. 'full' reuses analysis.py's cache."""
    if variant_name == "full":
        return ANALYSIS_CACHE_DIR
    return os.path.join(EDGE_STUDY_CACHE_DIR, variant_name)


def _get_regressors_for_variant(variant_name):
    """Return regressor configs with cache paths specific to a variant."""
    regressors = get_default_regressors()
    cache_dir = _variant_cache_dir(variant_name)
    os.makedirs(cache_dir, exist_ok=True)
    for r in regressors:
        r["file"] = os.path.join(cache_dir, os.path.basename(r["file"]))
    return regressors


def main():
    os.makedirs("plots", exist_ok=True)

    peak_sample = _load_peak_sample()
    process_rising_edge, process_falling_edge, process_first_20ns = (
        _make_edge_processors(peak_sample))

    variants = {
        "full": process_waveforms,
        "rising": process_rising_edge,
        "falling": process_falling_edge,
        "first_20ns": process_first_20ns,
    }

    all_cached = all(
        os.path.isdir(_variant_cache_dir(v)) and all(
            os.path.exists(os.path.join(_variant_cache_dir(v), f)) for f in [
                "test_alpha_features.pkl", "test_gamma_features.pkl",
                "test_waveforms.npz", "regressor_names.pkl"
            ]) for v in variants)

    if all_cached:
        print("All edge study caches found — skipping ROOT file loading.")
        alpha_waveforms = gamma_waveforms = None
        alpha_features = gamma_features = None
    else:
        print("Loading alpha data (Am-241)...")
        alpha_features, alpha_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Am241.root",
            array_branch="Samples",
        )
        print(f"Alpha events: {len(alpha_features)}, "
              f"waveform shape: {alpha_waveforms.shape}")

        print("Loading gamma data (Na-22)...")
        gamma_features, gamma_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Na22.root",
            array_branch="Samples",
        )
        print(f"Gamma events: {len(gamma_features)}, "
              f"waveform shape: {gamma_waveforms.shape}")

    for name, process_func in variants.items():
        print(f"  Running variant: {name}")

        cache_dir = _variant_cache_dir(name)
        regressors = _get_regressors_for_variant(name)

        regress_waveforms(
            (alpha_waveforms, gamma_waveforms),
            (alpha_features, gamma_features),
            regressors=regressors,
            process_func=process_func,
            cache_dir=cache_dir,
            plot_prefix=f"{name}_",
            skip_plots=True,
        )

    _write_auc_table(list(variants.keys()))

    print("Edge study complete")


VARIANT_LABELS = {
    "full": "Full Waveform",
    "rising": "Rising Edge",
    "falling": "Falling Edge",
    "first_20ns": "First 20 ns",
}

LO_LOWER = {"alpha": 375, "gamma": 0}
LO_UPPER = {"alpha": 1575, "gamma": 1750}


def _write_auc_table(variant_names):
    """Compute AUC per model per variant and emit a LaTeX table.

    Each cell is the inverse-variance-weighted mean of the per-LO-bin AUC
    values, matching the per-bin balancing used in the AUC-vs-LO plot.
    """
    model_names = [r["name"] for r in get_default_regressors()]

    auc_table = {}
    err_table = {}
    for variant in variant_names:
        cache_dir = _variant_cache_dir(variant)
        alpha = pd.read_pickle(
            os.path.join(cache_dir, "test_alpha_features.pkl"))
        gamma = pd.read_pickle(
            os.path.join(cache_dir, "test_gamma_features.pkl"))

        alpha = alpha[(alpha["light_output"] >= LO_LOWER["alpha"])
                      & (alpha["light_output"] <= LO_UPPER["alpha"])]
        gamma = gamma[(gamma["light_output"] >= LO_LOWER["gamma"])
                      & (gamma["light_output"] <= LO_UPPER["gamma"])]

        print(f"  Variant {variant}: computing per-LO-bin AUCs")
        _, _, per_bin_auc, per_bin_err, _ = _compute_per_lo_auc(
            alpha, gamma, model_names, include_legacy_psd=False)

        auc_table[variant] = {}
        err_table[variant] = {}
        for name in model_names:
            mean, err = error_weighted_auc(per_bin_auc[name],
                                           per_bin_err[name])
            auc_table[variant][name] = mean
            err_table[variant][name] = err

    lines = []
    lines.append(r"\begin{table}[h!]")
    lines.append(r"    \centering")
    col_spec = "l" + "c" * len(variant_names)
    lines.append(r"    \begin{tabular}{" + col_spec + "}")
    lines.append(r"        \hline \hline")
    header = (r"        \textbf{Model} & " +
              " & ".join(rf"\textbf{{{VARIANT_LABELS[v]}}}"
                         for v in variant_names) + r" \\")
    lines.append(header)
    lines.append(r"        \hline")
    for name in model_names:
        row_vals = " & ".join(
            rf"{auc_table[v][name]:.3f} $\pm$ {err_table[v][name]:.3f}"
            for v in variant_names)
        lines.append(f"        {name} & {row_vals} " + r"\\")
    lines.append(r"        \hline \hline")
    lines.append(r"    \end{tabular}")
    lines.append(r"    \caption{AUC values for each model across different "
                 r"waveform regions}")
    lines.append(r"    \label{tab:auc}")
    lines.append(r"\end{table}")

    table_str = "\n".join(lines)
    output_path = os.path.join(EDGE_STUDY_CACHE_DIR, "auc_table.txt")
    os.makedirs(EDGE_STUDY_CACHE_DIR, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(table_str + "\n")
    print(f"\nLaTeX AUC table saved to {output_path}")
    print(table_str)


if __name__ == "__main__":
    main()
