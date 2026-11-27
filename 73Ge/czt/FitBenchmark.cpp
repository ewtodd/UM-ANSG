// How long production's two kinds of CZT fit take on the RooFit backend of
// the build they run in, for the dissertation's appendix on the CUDA
// photopeak kernels of Analysis-Utilities. Each fit runs as production runs
// it:
//   direct  Method 2's direct fit of the 01/15 no-shield signal run: one
//           channel, a cold start, the Ge tail held (GainTransferFit.hpp's
//           FitDirectGeOnce).
//   insitu  Method 1's in-situ fit of the 01/13 Cu-shield pair: the signal
//           and background runs together, started from production's saved
//           post-cal state (RunBkgGeSim).
// Each fit runs once untimed, where a CUDA build creates its context and
// loads its kernels, then N_TIMED times. A time is the whole fit call:
// building the model, minimizing, the errors and drawing the fit's figure,
// the same work on every backend. The CPU time over the wall time is how many
// cores the fit kept busy. Ge must come out the same in every configuration.
//
// The backend is the build's: flake.nix's use73GeCuda switch gives the
// 73Ge shell au.cuda (RooFit's CUDA backend) when true and au.default (its
// CPU backend, Analysis-Utilities' own lineshape terms spread over OpenMP
// threads) when false. The cores are those the process may run on, which
// also sizes OpenMP's thread pool:
//   root -b -q 'FitBenchmark.cpp++'                 each build
//   taskset -c 0 root -b -q 'FitBenchmark.cpp++'    the CPU build, one core
// ("++" rebuilds the macro for the build it runs in.) Each run writes
// results/fit_benchmark_<cuda|cpu>_<cores>.result; FitBenchmarkTable.cpp
// writes the dissertation's table from them. Nothing else should use the GPU
// or the cores meanwhile.
#include "GainTransferFit.hpp"
#include "InSituFitState.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"

#include <TF1.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>

#include <sched.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

const Int_t N_TIMED = 3;

const TString DIRECT_DATE = "20260115";
const TString INSITU_PAIR = "CuShield10";

TString BackendName() {
#if defined(AU_ROOFIT_BACKEND_CUDA) && AU_ROOFIT_BACKEND_CUDA
  return "cuda";
#else
  return "cpu";
#endif
}

// The cores this process may run on, which taskset restricts.
Int_t AllowedCores() {
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) != 0)
    return 0;
  return CPU_COUNT(&set);
}

TString CpuModel() {
  std::ifstream in("/proc/cpuinfo");
  std::string line;
  while (std::getline(in, line)) {
    if (line.compare(0, 10, "model name") != 0)
      continue;
    size_t colon = line.find(':');
    if (colon == std::string::npos)
      continue;
    TString model = line.substr(colon + 1).c_str();
    model = model.Strip(TString::kBoth);
    return model;
  }
  return "unknown";
}

TString GpuModel() {
  if (BackendName() != "cuda")
    return "none";
  TString name = gSystem->GetFromPipe(
      "nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null");
  name = name.Strip(TString::kBoth);
  return name.Length() > 0 ? name : TString("unknown");
}

// A calibration's keV per raw unit: a day's (20260115) or a 01/13 pair's
// (20260113_<pair>).
Double_t CalibrationSlope(const TString &calibration_tag) {
  TFile *file = IO::OpenForReading(OUT_SUBDIR + "/calibration_function_low_" +
                                   calibration_tag + ".root");
  Double_t slope = 0;
  if (file && !file->IsZombie()) {
    TF1 *cal = static_cast<TF1 *>(file->Get("calibration"));
    if (cal)
      slope = cal->GetParameter(1);
    file->Close();
  }
  delete file;
  return slope;
}

Long64_t CountInWindow(const std::vector<Double_t> &events, Double_t lo,
                       Double_t hi) {
  Long64_t n = 0;
  for (size_t i = 0; i < events.size(); i++)
    if (events[i] >= lo && events[i] < hi)
      n++;
  return n;
}

Bool_t CopyFile(const TString &source, const TString &target) {
  std::ifstream in(source.Data(), std::ios::binary);
  std::ofstream out(target.Data(), std::ios::binary);
  if (!in.is_open() || !out.is_open())
    return kFALSE;
  out << in.rdbuf();
  return kTRUE;
}

// One fit's timed runs and the Ge it gives.
struct Timing {
  std::vector<Double_t> wall, cpu;
  Double_t ge_kev = 0, ge_err_kev = 0;
  Bool_t valid = kTRUE;
};

Double_t Sum(const std::vector<Double_t> &v) {
  Double_t s = 0;
  for (size_t i = 0; i < v.size(); i++)
    s += v[i];
  return s;
}

Double_t Smallest(const std::vector<Double_t> &v) {
  Double_t m = v.empty() ? 0 : v[0];
  for (size_t i = 1; i < v.size(); i++)
    if (v[i] < m)
      m = v[i];
  return m;
}

Double_t Largest(const std::vector<Double_t> &v) {
  Double_t m = v.empty() ? 0 : v[0];
  for (size_t i = 1; i < v.size(); i++)
    if (v[i] > m)
      m = v[i];
  return m;
}

// TStopwatch's readers are not const.
void Record(Timing &timing, Int_t run, TStopwatch &watch, Bool_t converged,
            Double_t ge, Double_t ge_err, const TString &name) {
  std::cout << "BENCH " << name << " run " << run << ": wall " << std::fixed
            << std::setprecision(2) << watch.RealTime() << " s, cpu "
            << watch.CpuTime() << " s, Ge " << std::setprecision(5) << ge
            << " +/- " << ge_err << " keV" << (converged ? "" : " [FAILED]")
            << (run == 0 ? " (untimed)" : "") << std::endl;
  if (!converged)
    timing.valid = kFALSE;
  if (run == 0)
    return;
  timing.wall.push_back(watch.RealTime());
  timing.cpu.push_back(watch.CpuTime());
  timing.ge_kev = ge;
  timing.ge_err_kev = ge_err;
}

// Method 2 as production runs it (CalibrateAndFit.cpp's FitDirectGe): the
// run's 95.902 keV line measured once for its escapes, then the held fit.
Timing TimeDirectFit(const TString &project_root, Long64_t &n_window) {
  Timing timing;
  const TString sig_run = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;
  std::vector<Double_t> all = LoadCalibrated(sig_run);
  std::vector<Double_t> events;
  for (size_t i = 0; i < all.size(); i++)
    if (InDirectFitRange(all[i]))
      events.push_back(all[i]);
  std::vector<Double_t>().swap(all);
  n_window = CountInWindow(events, DIRECT_WINDOW_LO, DIRECT_WINDOW_HI);
  const Double_t bin_kev = COMB_PERIOD_RAW_KEV * CalibrationSlope(DIRECT_DATE);
  RateSubShape shape;
  const TString tail_path = DirectTailPath(
      project_root + "/plots/czt/chain/" + RESULT_TAG, DIRECT_DATE);
  if (events.empty() || !(bin_kev > 0) || !ReadDirectTail(tail_path, shape)) {
    std::cerr << "BENCH-FATAL: no events, calibration or held tail ("
              << tail_path << ") for the direct fit" << std::endl;
    timing.valid = kFALSE;
    return timing;
  }
  shape.escapes = EscapeLinesRequested() ? &ESCAPES_USED : nullptr;
  std::vector<Double_t> escape_mus, escape_yields;
  if (shape.escapes) {
    CaptureLine line =
        MeasureDirectCaptureLine(events, bin_kev, shape, "benchmark");
    if (!line.valid) {
      std::cerr << "BENCH-FATAL: no 95.902 keV line" << std::endl;
      timing.valid = kFALSE;
      return timing;
    }
    EscapeLinesOf(line, *shape.escapes, 1.0, escape_mus, escape_yields);
  }
  for (Int_t run = 0; run <= N_TIMED; run++) {
    TStopwatch watch;
    watch.Start();
    FitResult fit =
        FitDirectGeOnce(events, bin_kev, shape, escape_mus, escape_yields,
                        kTRUE, nullptr, sig_run, "Benchmark_direct");
    watch.Stop();
    const Bool_t converged = fit.valid && !fit.peaks.empty();
    Record(timing, run, watch, converged, converged ? fit.peaks[0].mu : 0.0,
           converged ? fit.peaks[0].mu_error : 0.0, "direct");
  }
  return timing;
}

// Method 1's post-cal fit of the Cu-shield pair, as production calls it
// (CalibrateAndFit.cpp), started each time from production's saved state,
// which the fitter finds beside this macro's own fits.
Timing TimeInSituFit(const TString &project_root, const TString &bench_dir,
                     Long64_t &n_sig_window, Long64_t &n_bkg_window) {
  Timing timing;
  const TString bkg_run = Constants::CUSHIELDBACKGROUND_10PERCENT_20260113;
  const TString sig_run = Constants::CUSHIELDSIGNAL_10PERCENT_20260113;
  const TString bkg_label = bkg_run + "_postcal_" + INSITU_PAIR;
  const TString sig_label = sig_run + "_postcal_" + INSITU_PAIR;
  const TString fit_label = "PostCal_" + INSITU_PAIR;
  const TString production_state = project_root + "/plots/czt/chain/" +
                                   RESULT_TAG + "/fits/" + fit_label + "_" +
                                   sig_label + ".simroofits";
  const TString start_state =
      bench_dir + "/fits/" + fit_label + "_" + sig_label + ".simroofits";
  FitState state = ReadFitState(production_state);
  std::vector<Double_t> bkg_events = LoadCalibrated(bkg_run);
  std::vector<Double_t> sig_events = LoadCalibrated(sig_run);
  const Double_t bin_kev =
      COMB_PERIOD_RAW_KEV * CalibrationSlope("20260113_" + INSITU_PAIR);
  if (!state.ok || bkg_events.empty() || sig_events.empty() || !(bin_kev > 0)) {
    std::cerr << "BENCH-FATAL: no state (" << production_state
              << "), events or calibration for the in-situ fit" << std::endl;
    timing.valid = kFALSE;
    return timing;
  }
  n_sig_window = CountInWindow(sig_events, state.range_lo, state.range_hi);
  n_bkg_window = CountInWindow(bkg_events, state.range_lo, state.range_hi);
  const FitResult bkg_seed = ChannelSeed(state, "bkg");
  const FitResult sig_seed = ChannelSeed(state, "sig");
  InSituEscapes escapes;
  escapes.fractions = EscapeLinesRequested() ? &ESCAPES_USED : nullptr;
  escapes.units_per_kev = 1.0;
  escapes.comb_period = bin_kev;
  std::vector<Double_t> pb_mus;
  pb_mus.push_back((Double_t)E_PB_KA2);
  pb_mus.push_back((Double_t)E_PB_KA1);
  for (Int_t run = 0; run <= N_TIMED; run++) {
    // The fit saves its converged state where it found the start; each run
    // starts from production's again.
    if (!CopyFile(production_state, start_state)) {
      std::cerr << "BENCH-FATAL: cannot copy the state to " << start_state
                << std::endl;
      timing.valid = kFALSE;
      return timing;
    }
    TStopwatch watch;
    watch.Start();
    BkgGeSimResult fit =
        RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_label, sig_label,
                    fit_label, bin_kev, &bkg_seed, &sig_seed, escapes);
    watch.Stop();
    const Bool_t converged = fit.valid && !fit.sig_channel.peaks.empty();
    Record(timing, run, watch, converged,
           converged ? fit.sig_channel.peaks.back().mu : 0.0,
           converged ? fit.sig_channel.peaks.back().mu_error : 0.0, "insitu");
  }
  return timing;
}

void WriteTiming(std::ofstream &out, const TString &name,
                 const Timing &timing) {
  const Double_t wall = Sum(timing.wall);
  out << name << "_converged " << (timing.valid ? 1 : 0) << std::endl;
  out << name << "_wall_mean " << wall / timing.wall.size() << std::endl;
  out << name << "_wall_min " << Smallest(timing.wall) << std::endl;
  out << name << "_wall_max " << Largest(timing.wall) << std::endl;
  out << name << "_cores_busy " << (wall > 0 ? Sum(timing.cpu) / wall : 0)
      << std::endl;
  out << name << "_ge_kev " << std::setprecision(8) << timing.ge_kev
      << std::endl;
  out << name << "_ge_err_kev " << timing.ge_err_kev << std::endl;
  out << std::setprecision(4);
}

} // namespace

void FitBenchmark() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  const TString backend = BackendName();
  const Int_t cores = AllowedCores();
  const TString config = TString::Format("%s_%d", backend.Data(), cores);
  const TString bench_dir =
      project_root + "/plots/czt/studies/fit_benchmark/" + config;
  gSystem->mkdir(bench_dir + "/fits", kTRUE);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, bench_dir,
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  std::cout << "BENCH-CONFIG " << config << ", result tag " << RESULT_TAG
            << ", CPU " << CpuModel() << ", GPU " << GpuModel() << std::endl;

  Long64_t n_direct = 0, n_sig = 0, n_bkg = 0;
  const Timing direct = TimeDirectFit(project_root, n_direct);
  const Timing insitu = TimeInSituFit(project_root, bench_dir, n_sig, n_bkg);
  if (direct.wall.empty() || insitu.wall.empty()) {
    std::cerr << "BENCH-FATAL: a fit did not run; no result written"
              << std::endl;
    return;
  }

  const TString path =
      project_root + "/results/fit_benchmark_" + config + ".result";
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    std::cerr << "BENCH-FATAL: cannot write " << path << std::endl;
    return;
  }
  out << "# FitBenchmark.cpp: production's CZT fits on the " << backend
      << " backend with " << cores << " cores allowed; times in s over "
      << N_TIMED << " runs after one untimed" << std::endl;
  out << "backend " << backend << std::endl;
  out << "cores " << cores << std::endl;
  out << "cpu_model " << CpuModel() << std::endl;
  out << "gpu_model " << GpuModel() << std::endl;
  out << "direct_events " << n_direct << std::endl;
  out << "insitu_sig_events " << n_sig << std::endl;
  out << "insitu_bkg_events " << n_bkg << std::endl;
  out << std::fixed << std::setprecision(4);
  WriteTiming(out, "direct", direct);
  WriteTiming(out, "insitu", insitu);
  std::cout << "BENCH wrote " << path << std::endl;
}
