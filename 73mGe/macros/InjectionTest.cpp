#include "GainTransferFit.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitPhotopeakPdfs.hpp"

#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooMsgService.h>
#include <RooPolynomial.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RtypesCore.h>
#include <TCollection.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TList.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Edit these and re-run. No command-line arguments by design.
const TString INJ_SHAPE = "both";  // "model" | "empirical" | "both"
const Int_t INJ_PAIR_INDEX = 0;    // -1 = every pair
const Int_t INJ_ENERGY_INDEX = -1; // -1 = every injection energy
const Int_t INJ_N_TRIALS = 5;
const Int_t INJ_TRIAL_OFFSET = 0;
const TString INJ_RUN_TAG = "extended"; // results go to plots/injection_<tag>/
const TString INJ_SOURCE = "halves";    // "halves" | "real"
const Bool_t INJ_SYNTHETIC = kFALSE;
// Low linear tail in the injected/synthetic lineshape. Note the shared
// RunBkgGeSim does not take this flag (production fixes it), so it only affects
// how empirical/synthetic events are generated, not the fit.
const Bool_t INJ_USE_LOW_LIN = kFALSE;
// Keep the per-trial fit/pull PNGs. Off by default: at 72 cells x 20 trials it
// is thousands of files nobody reads, and painting them is a real cost. Turn on
// to inspect a single cell (set INJ_PAIR_INDEX/INJ_ENERGY_INDEX/INJ_N_TRIALS
// down too).
const Bool_t INJ_KEEP_TRIAL_PLOTS = kTRUE;
// Run the simultaneous fit with the hard settings (strategy 2, final HESSE,
// 200k calls) through the shared RunBkgGeSim FIT_HARD hook, so there is one
// implementation of those settings. Off reproduces the production defaults.
const Bool_t INJ_FIT_HARD = kTRUE;
// Rate-sub fit variants. The rate-sub fit is the cheap half of a trial, so
// each trial runs it against EVERY (range, tail) combination below on the SAME
// injected data -- the expensive simultaneous fit is done once per trial and
// shared. Each variant is written as its own CSV row, tagged in the method
// column, and the summarizer groups them separately. Production is 62.0-71.5
// with no tail; the scan spans the production range, the old 62.0-76.0 window,
// and a wider 60.0-78.0, each with and without the low-energy tail.
struct RsVariant {
  Double_t lo;
  Double_t hi;
  Bool_t tail;
  TString method;
};
const std::vector<RsVariant> INJ_RS_VARIANTS = {
    {62.0, 71.5, kFALSE, "rs_62_0_71_5_notail"},
    {62.0, 71.5, kTRUE, "rs_62_0_71_5_tail"},
    {62.0, 76.0, kFALSE, "rs_62_0_76_0_notail"},
    {62.0, 76.0, kTRUE, "rs_62_0_76_0_tail"},
    {60.0, 78.0, kFALSE, "rs_60_0_78_0_notail"},
    {60.0, 78.0, kTRUE, "rs_60_0_78_0_tail"},
};

namespace {

const Double_t kInjEnergies[9] = {68.650, 68.700, 68.725, 68.750, 68.775,
                                  68.800, 68.850, 66.5,   70.5};

struct InjPair {
  TString tag;   // state-file tag, e.g. CuShield10
  TString key;   // unique result key, e.g. CuShield10_0114
  TString label; // plot label
  TString bkg_run;
  TString sig_run;
  TString am_run;
};

std::vector<InjPair> MakePairs() {
  std::vector<InjPair> pairs;
  pairs.push_back({"CuShield10", "CuShield10", "Cu 10% (01/13)",
                   Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                   Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113});
  pairs.push_back({"CdShield10", "CdShield10", "Cd 10% (01/13)",
                   Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                   Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113});
  pairs.push_back({"CdShield25", "CdShield25", "Cd 25% (01/13)",
                   Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                   Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113});
  pairs.push_back({"CuShield10", "CuShield10_0114", "Cu 10% (01/14)",
                   Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                   Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                   Constants::POSTREACTOR_AM241_20260113});
  return pairs;
}

TString ProjectRoot() {
  return TString(gSystem->Getenv("HOME")) + "/UM-ANSG/73mGe";
}

struct SavedState {
  Bool_t ok = kFALSE;
  Double_t range_lo = 0, range_hi = 0;
  std::map<std::string, Double_t> values;
  Double_t Get(const std::string &key, Double_t fallback = -1e30) const {
    std::map<std::string, Double_t>::const_iterator it = values.find(key);
    return it == values.end() ? fallback : it->second;
  }
};

SavedState ReadSavedState(const TString &path) {
  SavedState s;
  std::ifstream in(path.Data());
  if (!in.is_open())
    return s;
  std::string token;
  if (!(in >> token) || token != "RANGE")
    return s;
  if (!(in >> s.range_lo >> s.range_hi))
    return s;
  std::string name;
  Double_t value, error;
  Int_t fixed;
  while (in >> name >> value >> error >> fixed)
    s.values[name] = value;
  s.ok = s.values.count("sig:Mu3") && s.values.count("sig:Sigma3") &&
         s.values.count("sig:GausAmplitude3");
  return s;
}

struct GeShape {
  Bool_t ok = kFALSE;
  Double_t sigma = 0;
  Double_t tau_ratio = 0;
  Double_t lin_slope = 0;
  // The low exponential tail amplitude as the saved state holds it: the
  // full-range ratio R when the file carries Analysis-Utilities'
  // ExpTailAmplitudeConvention marker, the old in-window ratio r otherwise.
  // InWindowExpRatio() gives a particular peak's in-window ratio.
  Double_t exp_amplitude = 0;
  Bool_t exp_is_full_range = kFALSE;
  Double_t r_lin = 0;
  Double_t range_lo = 0, range_hi = 0;
  Double_t gaus_yield = 0;
  Double_t area = 0;
};

// The tail's yield inside the saved window over the Gaussian's, for a peak at
// mu with width sigma. An old in-window ratio is returned as saved: the old
// fits shared that number between peaks as it is.
Double_t InWindowExpRatio(const GeShape &g, Double_t mu, Double_t sigma) {
  if (!g.exp_is_full_range)
    return g.exp_amplitude;
  return g.exp_amplitude *
         LowExpTailWindowFactor(g.range_lo, g.range_hi, mu, sigma, g.tau_ratio);
}

GeShape ExtractGeShape(const SavedState &s) {
  GeShape g;
  if (!s.ok)
    return g;
  g.sigma = s.Get("sig:Sigma3");
  g.gaus_yield = s.Get("sig:GausAmplitude3");
  g.exp_amplitude =
      s.Get("bkg:LowExpTailAmplitude1", s.Get("bkg:LowExpTailAmplitude2", 0.0));
  g.exp_is_full_range = s.Get("ExpTailAmplitudeConvention", 0.0) >= 1.0;
  g.r_lin =
      s.Get("bkg:LowLinTailAmplitude1", s.Get("bkg:LowLinTailAmplitude2", 0.0));
  g.tau_ratio =
      s.Get("bkg:LowExpTailRatio1", s.Get("bkg:LowExpTailRatio2", 0.0));
  g.lin_slope =
      s.Get("bkg:LowLinTailSlope1", s.Get("bkg:LowLinTailSlope2", 0.0));
  g.range_lo = s.range_lo;
  g.range_hi = s.range_hi;
  g.area = g.gaus_yield *
           (1.0 + InWindowExpRatio(g, s.Get("sig:Mu3"), g.sigma) + g.r_lin);
  g.ok = g.sigma > 0 && g.gaus_yield > 0 && g.tau_ratio > 0 &&
         g.range_hi > g.range_lo;
  return g;
}

FitResult MakeBkgSeed(const SavedState &s) {
  FitResult seed;
  seed.peaks.resize(2);
  for (Int_t i = 0; i < 2; i++) {
    TString suffix = TString::Format("%d", i + 1);
    seed.peaks[i].mu = s.Get(("bkg:Mu" + suffix).Data());
    seed.peaks[i].sigma =
        s.Get(("bkg:Sigma" + suffix).Data(), s.Get("bkg:Sigma1"));
    seed.peaks[i].gaus_amplitude = s.Get(("bkg:GausAmplitude" + suffix).Data());
  }
  seed.bkg_constant = s.Get("bkg:BkgConstant", 0.0);
  seed.lin_bkg_slope = s.Get("bkg:BkgSlope", 0.0);
  seed.valid = kTRUE;
  return seed;
}

// Signal-channel seed: the two Pb peaks plus the Ge peak moved to e_inj.
// Production's postcal fit seeds BOTH channels; passing only the bkg seed
// leaves the sig channel to start from whatever the state load gives it, so the
// linked parameters (LinkPeakShape links Mu, Sigma and the tails) begin
// inconsistent and the joint fit has to reconcile them inside its call budget.
FitResult MakeSigSeed(const SavedState &s, Double_t e_inj) {
  FitResult seed;
  seed.peaks.resize(3);
  for (Int_t i = 0; i < 3; i++) {
    TString suffix = TString::Format("%d", i + 1);
    seed.peaks[i].mu = s.Get(("sig:Mu" + suffix).Data());
    seed.peaks[i].sigma =
        s.Get(("sig:Sigma" + suffix).Data(), s.Get("sig:Sigma1"));
    seed.peaks[i].gaus_amplitude = s.Get(("sig:GausAmplitude" + suffix).Data());
  }
  seed.peaks[2].mu = e_inj;
  seed.bkg_constant = s.Get("sig:BkgConstant", 0.0);
  seed.lin_bkg_slope = s.Get("sig:BkgSlope", 0.0);
  seed.valid = kTRUE;
  return seed;
}

// Copy the real post-cal state to a trial path, replacing the Ge start value
// with E_inj + jitter. Never touches the source file.
Bool_t WriteTrialState(const TString &source, const TString &target,
                       Double_t ge_start_kev) {
  std::ifstream in(source.Data());
  std::ofstream out(target.Data());
  if (!in.is_open() || !out.is_open()) {
    std::cerr << "ERROR: cannot copy saved state " << source << " -> " << target
              << std::endl;
    return kFALSE;
  }
  std::string line;
  Bool_t replaced = kFALSE;
  while (std::getline(in, line)) {
    if (line.compare(0, 7, "sig:Mu3") == 0 &&
        (line.size() == 7 || line[7] == ' ' || line[7] == '\t')) {
      std::istringstream ss(line);
      std::string name;
      Double_t value, error;
      Int_t fixed;
      ss >> name >> value >> error >> fixed;
      out << name << " " << std::setprecision(17) << ge_start_kev << " "
          << error << " " << fixed << std::endl;
      replaced = kTRUE;
    } else {
      out << line << std::endl;
    }
  }
  return replaced;
}

Double_t RunLiveTime(const TString &run) {
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  if (!file || file->IsZombie()) {
    delete file;
    return 0;
  }
  Double_t total = 0;
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TParameter<Float_t> *p = static_cast<TParameter<Float_t> *>(
        file->Get(Form("LiveTime_Filtered_Crystal%d_s", c)));
    if (p)
      total += p->GetVal();
  }
  file->Close();
  delete file;
  return total;
}

std::vector<Double_t> ParityHalf(const std::vector<Double_t> &events,
                                 Int_t parity) {
  std::vector<Double_t> out;
  out.reserve(events.size() / 2 + 1);
  for (size_t i = (size_t)parity; i < events.size(); i += 2)
    out.push_back(events[i]);
  return out;
}

std::vector<Double_t> PoissonResample(const std::vector<Double_t> &events,
                                      TRandom3 &rng) {
  std::vector<Double_t> out;
  out.reserve(events.size() * 2);
  for (size_t i = 0; i < events.size(); i++) {
    Int_t copies = rng.Poisson(1.0);
    for (Int_t k = 0; k < copies; k++)
      out.push_back(events[i]);
  }
  return out;
}

// Fitted Ge lineshape: Gaussian plus the shared low exponential and linear
// tails. The saved tail amplitudes are ratios over the Gaussian yield.
std::vector<Double_t> SampleModel(Int_t n, Double_t e_inj,
                                  const GeShape &shape) {
  std::vector<Double_t> out;
  if (n <= 0)
    return out;
  RooRealVar x("x", "x", shape.range_lo, shape.range_hi);
  RooRealVar mu("mu", "mu", e_inj);
  RooRealVar sigma("sigma", "sigma", shape.sigma);
  RooRealVar tau("tau", "tau", shape.tau_ratio);
  RooRealVar slope("slope", "slope", shape.lin_slope);
  mu.setConstant(kTRUE);
  sigma.setConstant(kTRUE);
  tau.setConstant(kTRUE);
  slope.setConstant(kTRUE);
  RooGaussian gauss("gauss", "gauss", x, mu, sigma);
  RooLowExpTail lowexp("lowexp", "lowexp", x, mu, sigma, tau);
  RooLowLinTail lowlin("lowlin", "lowlin", x, mu, sigma, slope);
  Double_t r_exp = InWindowExpRatio(shape, e_inj, shape.sigma);
  Double_t denom = 1.0 + r_exp + shape.r_lin;
  RooConstVar f_gaus("f_gaus", "f_gaus", 1.0 / denom);
  RooConstVar f_exp("f_exp", "f_exp", r_exp / denom);
  RooConstVar f_lin("f_lin", "f_lin", shape.r_lin / denom);
  RooAddPdf total("total", "total", RooArgList(gauss, lowexp, lowlin),
                  RooArgList(f_gaus, f_exp, f_lin));
  RooDataSet *data = total.generate(RooArgSet(x), n);
  out.reserve(n);
  for (Int_t i = 0; i < data->numEntries(); i++)
    out.push_back(data->get(i)->getRealValue("x"));
  delete data;
  return out;
}

// Empirical Am-241 lineshape. The peak window is +-3 keV around the 59.54 keV
// line; the two +-3..+-6 keV sidebands give the background per bin. The
// sideband-subtracted histogram is sampled and each event mapped to
// E_inj + (E - mu_Am) * sigma_Ge / sigma_Am.
std::vector<Double_t> SampleEmpirical(Int_t n, Double_t e_inj,
                                      Double_t sigma_ge,
                                      const std::vector<Double_t> &am_events) {
  std::vector<Double_t> out;
  if (n <= 0 || am_events.empty())
    return out;
  const Double_t mu_am = E_AM241;
  const Double_t half_window = 3.0;
  const Double_t lo = mu_am - half_window;
  const Double_t hi = mu_am + half_window;
  const Double_t slo2 = mu_am + half_window;
  const Double_t shi2 = mu_am + 2.0 * half_window;
  const Int_t nbins = (Int_t)std::lround((hi - lo) / 0.05);
  const Int_t nbins_side = (Int_t)std::lround((shi2 - slo2) / 0.05);
  TH1D peak("am_peak", "Am peak", nbins, lo, hi);
  TH1D side_hi("am_side_hi", "Am high sideband", nbins_side, slo2, shi2);
  for (size_t i = 0; i < am_events.size(); i++) {
    Double_t e = am_events[i];
    if (e >= lo && e < hi)
      peak.Fill(e);
    else if (e >= slo2 && e < shi2)
      side_hi.Fill(e);
  }
  // Background from the HIGH sideband only: the low sideband sits on the
  // Am peak's low-energy tail and subtracting it would remove part of the
  // signal, biasing the injected core high. The background is flat over this
  // few-keV window, so scale the high-side rate to the peak bin width.
  Double_t bkg_per_kev = side_hi.Integral() / (shi2 - slo2);
  Double_t bkg_per_bin = bkg_per_kev * peak.GetBinWidth(1);
  TH1D subtracted("am_subtracted", "Am sideband-subtracted", nbins, lo, hi);
  for (Int_t b = 1; b <= nbins; b++) {
    Double_t value = peak.GetBinContent(b) - bkg_per_bin;
    subtracted.SetBinContent(b, value > 0.0 ? value : 0.0);
  }
  if (subtracted.Integral() <= 0) {
    std::cerr << "ERROR: Am sideband subtraction left no entries" << std::endl;
    return out;
  }

  TF1 gaus("am_gaus", "gaus", lo, hi);
  gaus.SetParameters(subtracted.GetMaximum(),
                     subtracted.GetBinCenter(subtracted.GetMaximumBin()), 0.8);
  TFitResultPtr fit = subtracted.Fit(&gaus, "Q R N S");
  Double_t mu_gaus = gaus.GetParameter(1);
  Double_t sigma_fit = gaus.GetParameter(2);
  if (!fit.Get() || !fit->IsValid() || !(sigma_fit > 0.1 && sigma_fit < 2.0)) {
    sigma_fit = subtracted.GetRMS();
    std::cerr << "WARNING: Am Gaussian fit failed; using histogram RMS"
              << std::endl;
  }
  if (!(sigma_fit > 0.05)) {
    std::cerr << "ERROR: Am lineshape has no usable width" << std::endl;
    return out;
  }
  // The calibration anchors the Am fit's core to the tabulated line energy,
  // so align 59.5409 to E_inj. A plain Gaussian fit to the empirical shape is
  // pulled down by its low tail and would bias the injected peak high.
  const Double_t mu_anchor = E_AM241;
  std::cout << "AM-SHAPE sigma_gaus " << std::setprecision(5) << sigma_fit
            << " mu_gaus " << mu_gaus << " mode "
            << subtracted.GetBinCenter(subtracted.GetMaximumBin()) << " anchor "
            << mu_anchor << std::endl;

  out.reserve(n);
  for (Int_t i = 0; i < n; i++) {
    Double_t sampled = subtracted.GetRandom();
    out.push_back(e_inj + (sampled - mu_anchor) * (sigma_ge / sigma_fit));
  }
  return out;
}

// Production rate subtraction, routed through the shared core in
// GainTransferFit.hpp so the injection test and the measurement fit with the
// same code. Two deliberate differences from the production call site: the
// background start comes from the window edge farther from the injected peak,
// and the mu limit is relaxed to +-0.2 keV so the 70.5 keV injection is not
// clamped at the boundary.
RateSubResult FitRateSubtractedEvents(const std::vector<Double_t> &sig_events,
                                      const std::vector<Double_t> &bkg_events,
                                      Double_t live_sig, Double_t live_bkg,
                                      Double_t e_inj, const TString &tag,
                                      Double_t rs_lo, Double_t rs_hi,
                                      Bool_t rs_tail) {
  const Double_t lo = rs_lo, hi = rs_hi;
  // Build over a fixed wide window so the fit range can be scanned without
  // rebuilding the histograms (and so the background start has room).
  const Double_t hist_lo = 60.0, hist_hi = 80.0;
  TH1F *sig = RateHistFromEvents(sig_events, live_sig, TString("rs_sig_") + tag,
                                 hist_lo, hist_hi);
  TH1F *bkg = RateHistFromEvents(bkg_events, live_bkg, TString("rs_bkg_") + tag,
                                 hist_lo, hist_hi);
  TH1F *res =
      static_cast<TH1F *>(sig->Clone((TString("rs_res_") + tag).Data()));
  res->SetDirectory(nullptr);
  res->Add(bkg, -1.0);

  // Background start from the fit-range edge farther from the injected peak.
  Int_t left_lo = res->FindBin(lo + 0.2), left_hi = res->FindBin(lo + 1.2);
  Int_t right_lo = res->FindBin(hi - 1.2), right_hi = res->FindBin(hi - 0.2);
  Double_t left_sum = 0, right_sum = 0;
  Int_t left_n = 0, right_n = 0;
  for (Int_t b = left_lo; b <= left_hi; b++) {
    left_sum += res->GetBinContent(b);
    left_n++;
  }
  for (Int_t b = right_lo; b <= right_hi; b++) {
    right_sum += res->GetBinContent(b);
    right_n++;
  }
  Double_t left = left_n ? left_sum / left_n : 0.0;
  Double_t right = right_n ? right_sum / right_n : 0.0;
  Double_t bkg0 =
      (std::fabs(e_inj - (lo + 0.7)) > std::fabs(e_inj - (hi - 0.7))) ? left
                                                                      : right;

  TF1 *model = nullptr;
  RateSubResult out =
      FitRateSubHist(res, lo, hi, tag, &model, 0.2, bkg0, rs_tail);
  // Production plots the rate-sub fit with a pull panel; the shared fit core
  // deliberately does not plot, so build the same figure here. SaveFigure
  // resolves the "fits" subdirectory against the cell dir from InitUtils.
  if (model && INJ_KEEP_TRIAL_PLOTS) {
    res->GetXaxis()->SetRangeUser(lo, hi);
    const Int_t kNPts = 400;
    Double_t xstep = (hi - lo) / (kNPts - 1);
    TGraph *total_graph = new TGraph(kNPts);
    TGraph *bkg_graph = new TGraph(kNPts);
    for (Int_t i = 0; i < kNPts; i++) {
      Double_t xv = lo + i * xstep;
      total_graph->SetPoint(i, xv, model->Eval(xv));
      bkg_graph->SetPoint(i, xv, out.bkg0 + out.slope * (xv - 68.75));
    }
    bkg_graph->SetLineStyle(2);
    std::vector<TGraph *> components = {bkg_graph};
    // No inlay label: the tag is already in the file name, and the old
    // "Rate-sub <tag>" box in the upper-right was just clutter.
    PlottingUtils::PlotFitWithResiduals(res, total_graph, components, lo, hi,
                                        "rate_fit_" + tag, "fits", "", kFALSE);
    delete total_graph;
    delete bkg_graph;
  }
  delete model;
  delete sig;
  delete bkg;
  delete res;
  return out;
}

UInt_t TrialSeed(Int_t pair_index, Int_t energy_index, Int_t shape_index,
                 Int_t trial, const TString &run_tag) {
  UInt_t h = 2166136261u;
  UInt_t parts[4] = {(UInt_t)(pair_index + 1), (UInt_t)(energy_index + 1),
                     (UInt_t)(shape_index + 1), (UInt_t)(trial + 1)};
  for (Int_t i = 0; i < 4; i++) {
    h ^= parts[i];
    h *= 16777619u;
  }
  for (Int_t i = 0; i < run_tag.Length(); i++) {
    h ^= (UInt_t)run_tag[i];
    h *= 16777619u;
  }
  return h ? h : 1u;
}

struct RateSubOutcome {
  TString method;
  Double_t mu = -1, err = -1;
  Bool_t valid = kFALSE, rail = kFALSE;
};

struct TrialOutcome {
  Int_t n_inj = 0;
  Double_t jitter_eV = 0;
  Double_t sim_mu = -1, sim_err = -1;
  Int_t sim_status = -999, sim_cov_qual = -999;
  Double_t sim_edm = -1;
  Bool_t sim_valid = kFALSE, sim_rail = kFALSE;
  std::vector<RateSubOutcome> ratesub; // one entry per INJ_RS_VARIANTS variant
};

void CleanTrialPlots(const TString &dir, const TString &trial_tag) {
  TSystemDirectory directory(dir, dir);
  TList *files = directory.GetListOfFiles();
  if (!files)
    return;
  TIter next(files);
  while (TSystemFile *file = static_cast<TSystemFile *>(next())) {
    TString name = file->GetName();
    if (name.EndsWith(".png") && name.Contains(trial_tag))
      gSystem->Unlink(dir + "/" + name);
  }
}

// Fully synthetic closure. Both channels are generated from the same family
// of PDFs and the same parameter values the fit uses: the two Pb peaks, the
// flat background, the shared low-side tails, and a Ge peak at E_inj. This
// differs from the "model" injection, which adds a saved-state Ge peak to a
// REAL background half; here the Pb/background data are model-consistent too.

struct SyntheticData {
  std::vector<Double_t> bkg, sig;
};

RooRealVar *MakeOwnedVar(std::vector<RooAbsArg *> &owned, const char *name,
                         Double_t value) {
  RooRealVar *v = new RooRealVar(name, name, value);
  owned.push_back(v);
  return v;
}

RooConstVar *MakeOwnedConst(std::vector<RooAbsArg *> &owned, const char *name,
                            Double_t value) {
  RooConstVar *v = new RooConstVar(name, name, value);
  owned.push_back(v);
  return v;
}

std::vector<Double_t> GenerateSyntheticChannel(const SavedState &s,
                                               const GeShape &shape,
                                               Bool_t with_ge, Double_t e_inj,
                                               Double_t n_inj,
                                               Bool_t use_low_lin) {
  const Double_t r_lin = use_low_lin ? shape.r_lin : 0.0;
  const Double_t tau = shape.tau_ratio;
  const Double_t lin_slope = shape.lin_slope;
  RooRealVar x("x", "x", s.range_lo, s.range_hi);
  std::vector<RooAbsArg *> owned;
  RooArgList pdfs, coefs;
  for (Int_t i = 1; i <= 2; i++) {
    TString suffix = TString::Format("%d", i);
    Double_t mu = s.Get(("bkg:Mu" + suffix).Data());
    Double_t sigma = s.Get(("bkg:Sigma" + suffix).Data(), s.Get("bkg:Sigma1"));
    Double_t yield = 0.5 * s.Get(("bkg:GausAmplitude" + suffix).Data());
    if (!(mu > 0) || !(sigma > 0) || !(yield > 0))
      continue;
    RooRealVar *m = MakeOwnedVar(owned, ("mu_pb" + suffix).Data(), mu);
    RooRealVar *sg = MakeOwnedVar(owned, ("sigma_pb" + suffix).Data(), sigma);
    RooRealVar *t = MakeOwnedVar(owned, ("tau_pb" + suffix).Data(), tau);
    RooRealVar *sl =
        MakeOwnedVar(owned, ("slope_pb" + suffix).Data(), lin_slope);
    RooGaussian *g = new RooGaussian(("g_pb" + suffix).Data(), "", x, *m, *sg);
    RooLowExpTail *e =
        new RooLowExpTail(("e_pb" + suffix).Data(), "", x, *m, *sg, *t);
    owned.push_back(g);
    owned.push_back(e);
    pdfs.add(*g);
    coefs.add(*MakeOwnedConst(owned, ("yg_pb" + suffix).Data(), yield));
    pdfs.add(*e);
    coefs.add(*MakeOwnedConst(owned, ("ye_pb" + suffix).Data(),
                              yield * InWindowExpRatio(shape, mu, sigma)));
    if (use_low_lin) {
      RooLowLinTail *l =
          new RooLowLinTail(("l_pb" + suffix).Data(), "", x, *m, *sg, *sl);
      owned.push_back(l);
      pdfs.add(*l);
      coefs.add(
          *MakeOwnedConst(owned, ("yl_pb" + suffix).Data(), yield * r_lin));
    }
  }
  if (with_ge) {
    Double_t r_exp = InWindowExpRatio(shape, e_inj, shape.sigma);
    Double_t yield = n_inj / (1.0 + r_exp + r_lin);
    RooRealVar *m = MakeOwnedVar(owned, "mu_ge", e_inj);
    RooRealVar *sg = MakeOwnedVar(owned, "sigma_ge", shape.sigma);
    RooRealVar *t = MakeOwnedVar(owned, "tau_ge", tau);
    RooRealVar *sl = MakeOwnedVar(owned, "slope_ge", lin_slope);
    RooGaussian *g = new RooGaussian("g_ge", "", x, *m, *sg);
    RooLowExpTail *e = new RooLowExpTail("e_ge", "", x, *m, *sg, *t);
    owned.push_back(g);
    owned.push_back(e);
    pdfs.add(*g);
    coefs.add(*MakeOwnedConst(owned, "yg_ge", yield));
    pdfs.add(*e);
    coefs.add(*MakeOwnedConst(owned, "ye_ge", yield * r_exp));
    if (use_low_lin) {
      RooLowLinTail *l = new RooLowLinTail("l_ge", "", x, *m, *sg, *sl);
      owned.push_back(l);
      pdfs.add(*l);
      coefs.add(*MakeOwnedConst(owned, "yl_ge", yield * r_lin));
    }
  }
  Double_t bkg_yield = 0.5 * s.Get("bkg:BkgConstant", 0.0);
  if (bkg_yield > 0) {
    RooRealVar *bkg_slope = MakeOwnedVar(owned, "slope_bkg", 0.0);
    RooPolynomial *poly =
        new RooPolynomial("bkg_poly", "", x, RooArgList(*bkg_slope));
    owned.push_back(poly);
    pdfs.add(*poly);
    coefs.add(*MakeOwnedConst(owned, "bkg_yield", bkg_yield));
  }

  std::vector<Double_t> events;
  {
    RooAddPdf pdf("synthetic_pdf", "", pdfs, coefs);
    RooDataSet *data = pdf.generate(RooArgSet(x), RooFit::Extended(kTRUE));
    events.reserve(data->numEntries());
    for (Int_t i = 0; i < data->numEntries(); i++)
      events.push_back(data->get(i)->getRealValue("x"));
    delete data;
  }
  for (size_t i = 0; i < owned.size(); i++)
    delete owned[i];
  return events;
}

SyntheticData GenerateSynthetic(const SavedState &s, const GeShape &shape,
                                Double_t e_inj, Double_t n_inj,
                                Bool_t use_low_lin) {
  SyntheticData out;
  out.bkg =
      GenerateSyntheticChannel(s, shape, kFALSE, e_inj, n_inj, use_low_lin);
  out.sig =
      GenerateSyntheticChannel(s, shape, kTRUE, e_inj, n_inj, use_low_lin);
  return out;
}

// Copy the real state into the synthetic truth: half-run Pb/background
// yields, the Ge yield for N_inj, and the jittered Ge start. Fills the bkg
// seed with the same scaled values.
Bool_t WriteSyntheticState(const TString &source, const TString &target,
                           const SavedState &s, const GeShape &shape,
                           Double_t ge_start, Double_t n_inj,
                           Bool_t use_low_lin, FitResult &seed) {
  std::ifstream in(source.Data());
  std::ofstream out(target.Data());
  if (!in.is_open() || !out.is_open())
    return kFALSE;
  const Double_t r_lin = use_low_lin ? shape.r_lin : 0.0;
  const Double_t ge_yield =
      n_inj / (1.0 + InWindowExpRatio(shape, ge_start, shape.sigma) + r_lin);
  std::string line;
  Bool_t replaced = kFALSE;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string name;
    Double_t value, error;
    Int_t fixed;
    if (!(ss >> name >> value >> error >> fixed)) {
      out << line << std::endl;
      continue;
    }
    if (name == "sig:Mu3") {
      value = ge_start;
      replaced = kTRUE;
    } else if (name == "bkg:GausAmplitude1" || name == "bkg:GausAmplitude2") {
      value *= 0.5;
    } else if (name == "bkg:BkgConstant") {
      value *= 0.5;
    } else if (name == "sig:GausAmplitude1") {
      value = 0.5 * s.Get("bkg:GausAmplitude1");
    } else if (name == "sig:GausAmplitude2") {
      value = 0.5 * s.Get("bkg:GausAmplitude2");
    } else if (name == "sig:BkgConstant") {
      value = 0.5 * s.Get("bkg:BkgConstant");
    } else if (name == "sig:GausAmplitude3") {
      value = ge_yield;
    }
    out << name << " " << std::setprecision(17) << value << " " << error << " "
        << fixed << std::endl;
  }
  seed = MakeBkgSeed(s);
  for (size_t i = 0; i < seed.peaks.size(); i++) {
    if (seed.peaks[i].gaus_amplitude > 0)
      seed.peaks[i].gaus_amplitude *= 0.5;
  }
  if (seed.bkg_constant > 0)
    seed.bkg_constant *= 0.5;
  return replaced;
}

TrialOutcome
RunOneTrial(UInt_t seed, Int_t trial, Double_t e_inj, const TString &shape,
            const std::vector<Double_t> &bkg_pool,
            const std::vector<Double_t> &sig_pool,
            const std::vector<Double_t> &am_events,
            const std::vector<Double_t> &pb_mus, const GeShape &ge_shape,
            const SavedState &saved_state, const FitResult &bkg_seed,
            const TString &real_state, const TString &fits_dir,
            const TString &trial_tag, Double_t n_inj_expected,
            Double_t live_sig, Double_t live_bkg, Double_t range_lo,
            Double_t range_hi, Bool_t use_low_lin, Bool_t synthetic) {
  TrialOutcome out;
  TRandom3 rng(seed);
  Double_t jitter_kev = rng.Uniform(-0.050, 0.050);
  out.jitter_eV = jitter_kev * 1000.0;

  std::vector<Double_t> bkg_events, sig_events;
  FitResult trial_seed = bkg_seed;
  FitResult trial_sig_seed = MakeSigSeed(saved_state, e_inj + jitter_kev);
  if (synthetic) {
    out.n_inj = n_inj_expected > 0 ? (Int_t)std::lround(n_inj_expected) : 0;
    SyntheticData syn = GenerateSynthetic(saved_state, ge_shape, e_inj,
                                          n_inj_expected, use_low_lin);
    bkg_events = syn.bkg;
    sig_events = syn.sig;
  } else {
    out.n_inj = n_inj_expected > 0 ? rng.Poisson(n_inj_expected) : 0;
    bkg_events = PoissonResample(bkg_pool, rng);
    sig_events = PoissonResample(sig_pool, rng);
    std::vector<Double_t> injected;
    if (shape == "model")
      injected = SampleModel(out.n_inj, e_inj, ge_shape);
    else
      injected = SampleEmpirical(out.n_inj, e_inj, ge_shape.sigma, am_events);
    sig_events.insert(sig_events.end(), injected.begin(), injected.end());
  }

  TString fit_label = "Inj_" + trial_tag;
  TString sig_label = "Sig_" + trial_tag;
  TString bkg_label = "Bkg_" + trial_tag;
  TString trial_state =
      fits_dir + "/" + fit_label + "_" + sig_label + ".simroofits";
  Bool_t state_ok =
      synthetic ? WriteSyntheticState(real_state, trial_state, saved_state,
                                      ge_shape, e_inj + jitter_kev,
                                      n_inj_expected, use_low_lin, trial_seed)
                : WriteTrialState(real_state, trial_state, e_inj + jitter_kev);
  if (!state_ok) {
    std::cerr << "INJ-FATAL: could not write trial state " << trial_state
              << std::endl;
    out.sim_valid = kFALSE;
    return out;
  }

  BkgGeSimResult sim =
      RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_label, sig_label,
                  fit_label, &trial_seed, &trial_sig_seed);
  if (!sim.sig_channel.peaks.empty()) {
    const PeakFitResult &ge = sim.sig_channel.peaks.back();
    out.sim_mu = ge.mu;
    out.sim_err = ge.mu_error;
    out.sim_status = sim.sig_channel.fit_status;
    out.sim_cov_qual = sim.sig_channel.cov_qual;
    out.sim_edm = sim.sig_channel.edm;
    out.sim_valid = sim.valid && std::isfinite(ge.mu) && ge.mu_error > 0 &&
                    sim.sig_channel.has_fit_diagnostics;
    Double_t margin = 0.01 * (range_hi - range_lo);
    out.sim_rail = (ge.mu <= range_lo + margin || ge.mu >= range_hi - margin);
  }

  // Run the cheap rate-sub fit against every (range, tail) variant on the SAME
  // injected data. The expensive simultaneous fit above is shared across them,
  // so the variants are directly comparable trial-for-trial.
  for (size_t vi = 0; vi < INJ_RS_VARIANTS.size(); vi++) {
    const RsVariant &var = INJ_RS_VARIANTS[vi];
    RateSubOutcome ro;
    ro.method = var.method;
    RateSubResult rs = FitRateSubtractedEvents(
        sig_events, bkg_events, live_sig, live_bkg, e_inj,
        trial_tag + "_" + var.method, var.lo, var.hi, var.tail);
    ro.mu = rs.mu;
    ro.err = rs.mu_err;
    ro.valid = rs.valid && std::isfinite(rs.mu) && rs.mu_err > 0;
    ro.rail = (rs.mu <= var.lo + 0.35 || rs.mu >= var.hi - 0.35);
    out.ratesub.push_back(ro);
  }

  // Per-trial fit-on-histogram and residual figures are removed after each
  // trial unless INJ_KEEP_TRIAL_PLOTS is set, to avoid thousands of PNGs in a
  // scan.
  if (!INJ_KEEP_TRIAL_PLOTS) {
    CleanTrialPlots(fits_dir, fit_label);
    CleanTrialPlots(fits_dir + "/residual_hists", fit_label);
  }
  return out;
}

TString ShapeTag(const TString &shape) {
  return shape == "model" ? "model" : "emp";
}

TString CellTag(const InjPair &pair, const TString &shape, Double_t e_inj) {
  TString e_tag = Form("E%.3fkeV", e_inj);
  e_tag.ReplaceAll(".", "_");
  return pair.key + "_" + ShapeTag(shape) + "_" + e_tag;
}

void WriteRow(std::ofstream &csv, const InjPair &pair, const TString &shape,
              Double_t e_inj, Int_t trial, const TrialOutcome &out,
              const TString &method, Double_t mu, Double_t mu_err, Int_t status,
              Int_t cov_qual, Double_t edm, Bool_t valid, Bool_t rail,
              Double_t jitter_eV) {
  Double_t residual_eV =
      (mu > 0 && std::isfinite(mu)) ? (mu - e_inj) * 1000.0 : 0.0;
  csv << pair.key.Data() << "," << shape.Data() << "," << std::fixed
      << std::setprecision(4) << e_inj << "," << trial << "," << out.n_inj
      << "," << std::setprecision(3) << jitter_eV << "," << method.Data() << ","
      << std::setprecision(8) << mu << "," << mu_err << ","
      << std::setprecision(4) << residual_eV << "," << status << "," << cov_qual
      << "," << std::setprecision(6) << edm << "," << (valid ? 1 : 0) << ","
      << (rail ? 1 : 0) << std::endl;
}

} // namespace

void InjectionTest() {
  const TString &shape_arg = INJ_SHAPE;
  const Int_t pair_index = INJ_PAIR_INDEX;
  const Int_t energy_index = INJ_ENERGY_INDEX;
  const Int_t n_trials = INJ_N_TRIALS;
  const Int_t trial_offset = INJ_TRIAL_OFFSET;
  const TString &run_tag = INJ_RUN_TAG;
  const TString &source_arg = INJ_SOURCE;
  const Bool_t synthetic = INJ_SYNTHETIC;
  const Bool_t use_low_lin = INJ_USE_LOW_LIN;
  gROOT->SetBatch(kTRUE);
  // One implementation of the hard-fit settings: the shared RunBkgGeSim reads
  // FIT_HARD. Unset it when the control is off so a stale shell variable cannot
  // silently change which fit this test runs.
  if (INJ_FIT_HARD)
    gSystem->Setenv("FIT_HARD", "1");
  else
    gSystem->Unsetenv("FIT_HARD");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // "halves": both channels are parity halves of the pair's BACKGROUND run, so
  //           the Pb content is identical in the two channels and the
  //           LinkPeakShape assumption is satisfied by construction.
  // "real"  : the real signal run is the signal channel and the real background
  //           run is the background channel, so the sig/bkg Pb mismatch (and
  //           the non-cancelling Pb lines in the subtraction) are preserved.
  const Bool_t real_source = (source_arg == "real");
  if (!real_source && source_arg != "halves") {
    std::cerr << "ERROR: source must be 'halves' or 'real', got '" << source_arg
              << "'" << std::endl;
    return;
  }

  if (TString(RESULT_TAG) != "gaintransfer") {
    std::cerr << "ERROR: injection test expects RESULT_TAG 'gaintransfer', "
                 "but the shared header says "
              << RESULT_TAG << std::endl;
    return;
  }
  if (shape_arg != "model" && shape_arg != "empirical" && shape_arg != "both") {
    std::cerr << "ERROR: shape must be model, empirical, or both" << std::endl;
    return;
  }

  std::cout << "INJ-CONFIG low_lin " << (use_low_lin ? 1 : 0) << " synthetic "
            << (synthetic ? 1 : 0) << " shape " << shape_arg << " pair "
            << pair_index << " energy " << energy_index << " trials "
            << n_trials << " tag " << run_tag << std::endl;
  const TString project_root = ProjectRoot();
  const TString base = project_root + "/plots/injection_" + run_tag;
  gSystem->mkdir(base + "/results", kTRUE);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, base,
                                project_root + "/root_files");

  std::cout << "INJ-SOURCE " << source_arg << " ("
            << (real_source ? "real sig/bkg runs"
                            : "parity halves of the background run")
            << ")" << std::endl;

  std::vector<InjPair> pairs = MakePairs();
  std::vector<TString> shapes;
  if (shape_arg == "both") {
    shapes.push_back("model");
    shapes.push_back("empirical");
  } else {
    shapes.push_back(shape_arg);
  }

  std::vector<Double_t> am_events;
  if (!synthetic) {
    am_events = LoadCalibrated(Constants::POSTREACTOR_AM241_20260113);
    if (am_events.empty()) {
      std::cerr << "ERROR: no calibrated Am-241 events" << std::endl;
      return;
    }
  }

  for (Int_t pi = 0; pi < (Int_t)pairs.size(); pi++) {
    if (pair_index >= 0 && pi != pair_index)
      continue;
    const InjPair &pair = pairs[pi];

    TString state_path = project_root + "/plots/calibration_low_" + RESULT_TAG +
                         "/fits/PostCal_" + pair.tag + "_" + pair.sig_run +
                         "_postcal_" + pair.tag + ".simroofits";
    SavedState saved = ReadSavedState(state_path);
    if (!saved.ok) {
      std::cerr << "INJ-FATAL: missing or incomplete real post-cal state "
                << state_path << std::endl;
      return;
    }
    GeShape ge_shape = ExtractGeShape(saved);
    if (!ge_shape.ok) {
      std::cerr << "INJ-FATAL: cannot extract Ge shape from " << state_path
                << std::endl;
      return;
    }
    FitResult bkg_seed = MakeBkgSeed(saved);

    // Pools each trial resamples from.
    //   halves - both channels are parity halves of the background run.
    //   real   - signal channel is the real signal run, background channel the
    //            real background run, so the sig/bkg Pb content differs as it
    //            does in the actual measurement.
    std::vector<Double_t> bkg_pool, sig_pool;
    if (!synthetic) {
      std::vector<Double_t> bkg_all = LoadCalibrated(pair.bkg_run);
      if (bkg_all.empty()) {
        std::cerr << "INJ-FATAL: no calibrated background events for "
                  << pair.bkg_run << std::endl;
        return;
      }
      if (real_source) {
        sig_pool = LoadCalibrated(pair.sig_run);
        if (sig_pool.empty()) {
          std::cerr << "INJ-FATAL: no calibrated signal events for "
                    << pair.sig_run << std::endl;
          return;
        }
        bkg_pool = bkg_all;
      } else {
        bkg_pool = ParityHalf(bkg_all, 0);
        sig_pool = ParityHalf(bkg_all, 1);
      }
    }
    Double_t t_bkg = RunLiveTime(pair.bkg_run);
    Double_t t_sig = RunLiveTime(pair.sig_run);
    if (!(t_bkg > 0 && t_sig > 0)) {
      std::cerr << "INJ-FATAL: missing live time for " << pair.bkg_run << " ("
                << t_bkg << ") or " << pair.sig_run << " (" << t_sig << ")"
                << std::endl;
      return;
    }
    // In "halves" the injected yield is rescaled to the half-run live time so
    // the test sees a line of realistic strength on a background-only channel.
    // In "real" the sig channel already contains the real Ge line, so injecting
    // the same yield again roughly doubles it: the recovered mu is a blend of
    // the real line and the injection, and only the METHOD DIFFERENCE
    // (sim - ratesub) at fixed E_inj is interpretable, not the absolute
    // residual.
    Double_t live_sig = 0, live_bkg = 0, n_inj_expected = 0;
    if (real_source) {
      live_sig = t_sig;
      live_bkg = t_bkg;
      n_inj_expected = ge_shape.area;
    } else {
      live_sig = 0.5 * t_bkg;
      live_bkg = 0.5 * t_bkg;
      n_inj_expected = ge_shape.area * live_sig / t_sig;
    }
    std::cout << "INJ-PAIR " << pair.label << " real_ge_area "
              << std::setprecision(1) << ge_shape.area << " sigma "
              << std::setprecision(5) << ge_shape.sigma << " t_bkg " << t_bkg
              << " t_sig " << t_sig << " n_inj_expected " << n_inj_expected
              << std::endl;

    for (Int_t ei = 0; ei < 9; ei++) {
      if (energy_index >= 0 && ei != energy_index)
        continue;
      Double_t e_inj = kInjEnergies[ei];
      for (Int_t si = 0; si < (Int_t)shapes.size(); si++) {
        const TString &shape = shapes[si];
        // The cell no longer encodes a rate-sub variant: one cell now holds the
        // shared simultaneous fit plus every rate-sub variant (see
        // INJ_RS_VARIANTS), and the variant is distinguished in the CSV method
        // column instead, so the variants cannot pool into one average.
        TString cell_tag =
            CellTag(pair, shape, e_inj) + (real_source ? "_real" : "");
        TString cell_dir = base + "/cells/" + cell_tag;
        TString fits_dir = cell_dir + "/fits";
        gSystem->mkdir(fits_dir, kTRUE);
        InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, cell_dir,
                                      project_root + "/root_files");
        // Include the trial range so a pilot (offset 0, N trials) and a later
        // top-up (offset N, M trials) can share a run tag without overwriting
        // each other; the summariser reads every CSV in results/.
        TString csv_path = base + "/results/injection_" + cell_tag +
                           Form("_t%d_n%d", trial_offset, n_trials) + ".csv";
        std::ofstream csv(csv_path.Data());
        if (!csv.is_open()) {
          std::cerr << "INJ-FATAL: cannot write " << csv_path << std::endl;
          return;
        }
        csv << "pair,shape,energy_keV,trial,n_inj,jitter_eV,method,mu_keV,"
               "mu_err_keV,residual_eV,status,covQual,edm,valid,rail"
            << std::endl;
        std::cout << "INJ-CELL-START " << cell_tag << std::endl;
        for (Int_t trial = trial_offset; trial < trial_offset + n_trials;
             trial++) {
          UInt_t seed = TrialSeed(pi, ei, si, trial, run_tag);
          gRandom->SetSeed(seed);
          RooRandom::randomGenerator()->SetSeed(seed);
          TString trial_tag = cell_tag + "_t" + Form("%d", trial);
          TrialOutcome out = RunOneTrial(
              seed, trial, e_inj, shape, bkg_pool, sig_pool, am_events,
              {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1}, ge_shape, saved,
              bkg_seed, state_path, fits_dir, trial_tag, n_inj_expected,
              live_sig, live_bkg, saved.range_lo, saved.range_hi, use_low_lin,
              synthetic);
          WriteRow(csv, pair, shape, e_inj, trial, out, "sim", out.sim_mu,
                   out.sim_err, out.sim_status, out.sim_cov_qual, out.sim_edm,
                   out.sim_valid, out.sim_rail, out.jitter_eV);
          // One row per rate-sub variant; the method column carries the variant
          // so the summarizer can group each (range, tail) setting on its own.
          for (size_t vi = 0; vi < out.ratesub.size(); vi++) {
            const RateSubOutcome &ro = out.ratesub[vi];
            WriteRow(csv, pair, shape, e_inj, trial, out, ro.method, ro.mu,
                     ro.err, 0, -1, -1, ro.valid, ro.rail, 0.0);
          }
          std::cout << "INJ-TRIAL " << trial_tag << " n_inj " << out.n_inj
                    << " sim_mu " << std::setprecision(8) << out.sim_mu
                    << " sim_err " << out.sim_err << " valid "
                    << (out.sim_valid ? 1 : 0) << " n_ratesub "
                    << out.ratesub.size() << std::endl;
        }
        csv.close();
        std::cout << "INJ-CELL-DONE " << cell_tag << std::endl;
      }
    }
  }
}
