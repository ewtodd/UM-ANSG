// The testbed of the CZT Ge fits: a Ge line of known energy in data whose
// other content is real or generated, fitted by production's own in-situ
// (Method 1) and rate-subtracted (Method 2) code in GainTransferFit.hpp.
//
// Each pair first gets a reference fit: production's post-cal in-situ fit,
// rerun on the real pair from its saved state in the testbed's own
// directory. Its converged parameters are the testbed's truth (the shared
// lineshape, the Pb lines, the continuum) and, scaled to what a trial's
// channels hold, the start of every trial's fit.
//
// Sources (INJ_SOURCE):
//   halves    both channels are halves of the background run, drawn afresh
//             each trial, the Ge line injected into one; both methods run.
//             Every trial shares the run's own fluctuations, which the
//             in-situ fit sees in both channels; only the rate subtraction
//             cancels them.
//   real      the real signal and background runs, a second Ge line
//             injected into the signal run on top of its own; only the
//             difference between the methods at one energy means anything.
//   synthetic both channels generated from the reference fit at the pair's
//             statistics, Ge at the injected energy; with
//             INJ_SYNTHETIC_ESCAPES each run's 95.902 keV line and its K
//             X-ray escapes under the window's lines too. In-situ fit only.
//   synthetic_halves  the halves procedure on INJ_N_BACKGROUNDS synthetic
//             background runs, each generated once at the background run's
//             statistics and held through its trials. The spread of their
//             mean residuals is what the shared fluctuation alone does to
//             the halves; the real run's halves residual beyond it would
//             be the real background's own effect. In-situ fit only.
// With INJ_FIT_WITHOUT_ESCAPES every trial is fitted a second time without
// the escape lines ("sim_noesc"), on the same events.
//
// Usage (no arguments; edit the controls below):
//   root -b -q 'InjectionTest.cpp+'
// then SummarizeInjection.cpp with the same run tag. Run pairs as separate
// processes (INJ_PAIR_INDEX), each with its own reference fit.
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
const TString INJ_SHAPE = "model"; // "model" | "empirical" | "both"
const Int_t INJ_PAIR_INDEX = -1;   // -1 = every pair
const Int_t INJ_ENERGY_INDEX = -1; // -1 = every injection energy
const Int_t INJ_N_TRIALS = 10;
const Int_t INJ_TRIAL_OFFSET = 0;
const TString INJ_RUN_TAG =
    "closure"; // results go to plots/testbed/injection_<tag>/
const TString INJ_SOURCE =
    "halves"; // "halves" | "real" | "synthetic" | "synthetic_halves"
// synthetic_halves: how many synthetic background runs, each the pool of
// INJ_N_TRIALS halves trials.
const Int_t INJ_N_BACKGROUNDS = 8;
// Synthetic data also carry each run's 95.902 keV line and, under everything
// else, its K X-ray escapes at the pair's escape fractions.
const Bool_t INJ_SYNTHETIC_ESCAPES = kTRUE;
// Fit every trial a second time without the escape lines ("sim_noesc"), on
// the same data, so the summary gives the escapes' effect trial by trial.
const Bool_t INJ_FIT_WITHOUT_ESCAPES = kTRUE;
// Keep the per-trial fit/pull PNGs and fit states. Off by default: a full scan
// makes thousands of files nobody reads, and painting them is a real cost.
// Turn on to inspect a single cell (set INJ_PAIR_INDEX/INJ_ENERGY_INDEX/
// INJ_N_TRIALS down too).
const Bool_t INJ_KEEP_TRIAL_PLOTS = kFALSE;

namespace {

const Double_t kInjEnergies[9] = {68.650, 68.700, 68.725, 68.750, 68.775,
                                  68.800, 68.850, 66.5,   70.5};

struct InjPair {
  TString tag;   // state-file and result key, e.g. CuShield10
  TString label; // plot label
  TString bkg_run;
  TString sig_run;
  TString am_run;
  // Where the pair's 95.902 keV line comes from, as in production.
  const EscapeFractions *escapes;
};

std::vector<InjPair> MakePairs() {
  std::vector<InjPair> pairs;
  pairs.push_back({"CuShield10", "Cu 10% (01/13)",
                   Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                   Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113, &ESCAPES_INSIDE});
  pairs.push_back({"CdShield10", "Cd 10% (01/13)",
                   Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                   Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113, &ESCAPES_OUTSIDE});
  pairs.push_back({"CdShield25", "Cd 25% (01/13)",
                   Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                   Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                   Constants::POSTREACTOR_AM241_20260113, &ESCAPES_OUTSIDE});
  return pairs;
}

TString ProjectRoot() { return Paths::ProjectRootOf(__FILE__); }

std::vector<Double_t> PbKaEnergies() {
  std::vector<Double_t> mus;
  mus.push_back((Double_t)E_PB_KA2);
  mus.push_back((Double_t)E_PB_KA1);
  return mus;
}

// The pair calibration's keV per raw unit, which turns the comb period into
// calibrated energy.
Double_t PairCalibrationSlope(const TString &pair_tag) {
  TFile *file = IO::OpenForReading(
      OUT_SUBDIR + "/calibration_function_low_20260113_" + pair_tag + ".root");
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

// --------------------------------------------------------------- fit states

struct ParamValue {
  Double_t value = 0, error = 0;
};

// An in-situ fit's window and parameters by name ("bkg:Mu1", "sig:Sigma3",
// ...), as a .simroofits file holds them.
struct FitState {
  Bool_t ok = kFALSE;
  Double_t range_lo = 0, range_hi = 0;
  std::map<std::string, ParamValue> params;
  Bool_t Has(const std::string &name) const { return params.count(name) > 0; }
  Double_t Get(const std::string &name, Double_t fallback = -1e30) const {
    std::map<std::string, ParamValue>::const_iterator it = params.find(name);
    return it == params.end() ? fallback : it->second.value;
  }
  Double_t Error(const std::string &name) const {
    std::map<std::string, ParamValue>::const_iterator it = params.find(name);
    return it == params.end() ? 0.0 : it->second.error;
  }
  void Set(const std::string &name, Double_t value, Double_t error) {
    params[name].value = value;
    params[name].error = error;
  }
};

TString Indexed(const TString &channel, const char *param, Int_t index) {
  return channel + ":" + param + TString::Format("%d", index);
}

Bool_t HasGeLine(const FitState &s) {
  return s.Has("sig:Mu3") && s.Has("sig:Sigma3") && s.Has("sig:GausAmplitude3");
}

FitState ReadFitState(const TString &path) {
  FitState s;
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
    s.Set(name, value, error);
  s.ok = HasGeLine(s);
  return s;
}

// The in-situ model's own free parameters (GainTransferFit.hpp's
// FitBkgGeSimFromSeed): every Pb width and tail is bkg:Ka1's, so the signal
// channel adds only its centroids, yields, continuum and the Ge width.
std::vector<std::string> FreeStateParams() {
  const char *names[] = {"bkg:Mu1",
                         "bkg:Sigma1",
                         "bkg:GausAmplitude1",
                         "bkg:LowExpTailAmplitude1",
                         "bkg:LowExpTailRatio1",
                         "bkg:Mu2",
                         "bkg:GausAmplitude2",
                         "bkg:BkgConstant",
                         "sig:Mu1",
                         "sig:GausAmplitude1",
                         "sig:Mu2",
                         "sig:GausAmplitude2",
                         "sig:Mu3",
                         "sig:Sigma3",
                         "sig:GausAmplitude3",
                         "sig:BkgConstant"};
  std::vector<std::string> out(names, names + sizeof(names) / sizeof(names[0]));
  return out;
}

// Written in Analysis-Utilities' current convention, whose marker line says
// the exponential-tail amplitudes are full-range ratios.
Bool_t WriteFitState(const FitState &s, const TString &path) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    std::cerr << "ERROR: cannot write fit state " << path << std::endl;
    return kFALSE;
  }
  out << std::setprecision(17);
  out << "RANGE " << s.range_lo << " " << s.range_hi << std::endl;
  out << "ExpTailAmplitudeConvention 1 0 1" << std::endl;
  std::vector<std::string> names = FreeStateParams();
  for (size_t i = 0; i < names.size(); i++) {
    if (!s.Has(names[i]))
      continue;
    out << names[i] << " " << s.Get(names[i]) << " " << s.Error(names[i])
        << " 0" << std::endl;
  }
  return kTRUE;
}

FitState FitStateOf(const BkgGeSimResult &fit, Double_t range_lo,
                    Double_t range_hi) {
  FitState s;
  s.range_lo = range_lo;
  s.range_hi = range_hi;
  const FitResult *channels[2] = {&fit.bkg_channel, &fit.sig_channel};
  for (Int_t c = 0; c < 2; c++) {
    const std::vector<FitParameterDiagnostic> &diags =
        channels[c]->parameter_diagnostics;
    for (size_t i = 0; i < diags.size(); i++)
      s.Set(diags[i].name, diags[i].value, diags[i].error);
  }
  s.ok = fit.valid && HasGeLine(s);
  return s;
}

// Multiplies a yield by f, its error by sqrt(f), as for a run f times as long.
void ScaleYield(FitState &s, const std::string &name, Double_t f) {
  if (s.Has(name))
    s.Set(name, f * s.Get(name), TMath::Sqrt(f) * s.Error(name));
}

Double_t PbSigma(const FitState &s, Int_t index) {
  return s.Get(Indexed("bkg", "Sigma", index).Data(), s.Get("bkg:Sigma1"));
}

// ------------------------------------------------------------- lineshape

// The Ge line: Gaussian plus the low exponential tail every line shares, its
// amplitude R the tail's full-range yield over the Gaussian's.
struct GeShape {
  Bool_t ok = kFALSE;
  Double_t sigma = 0;
  Double_t tail_r = 0;
  Double_t tau_over_sigma = 0;
  Double_t range_lo = 0, range_hi = 0;
  Double_t gaus_yield = 0;
  Double_t area = 0; // Gaussian plus tail, in the window
};

// The tail's yield inside the window over the Gaussian's, for a line at mu
// with width sigma.
Double_t InWindowExpRatio(const GeShape &g, Double_t mu, Double_t sigma) {
  return g.tail_r * LowExpTailWindowFactor(g.range_lo, g.range_hi, mu, sigma,
                                           g.tau_over_sigma);
}

// The in-window share of a Gaussian at mu with width sigma.
Double_t GaussianWindowShare(Double_t lo, Double_t hi, Double_t mu,
                             Double_t sigma) {
  return 0.5 * (TMath::Erf((hi - mu) / (TMath::Sqrt2() * sigma)) -
                TMath::Erf((lo - mu) / (TMath::Sqrt2() * sigma)));
}

GeShape GeShapeOf(const FitState &s) {
  GeShape g;
  if (!s.ok)
    return g;
  g.sigma = s.Get("sig:Sigma3");
  g.gaus_yield = s.Get("sig:GausAmplitude3");
  g.tail_r = s.Get("bkg:LowExpTailAmplitude1");
  g.tau_over_sigma = s.Get("bkg:LowExpTailRatio1");
  g.range_lo = s.range_lo;
  g.range_hi = s.range_hi;
  g.area =
      g.gaus_yield * (1.0 + InWindowExpRatio(g, s.Get("sig:Mu3"), g.sigma));
  g.ok = g.sigma > 0 && g.gaus_yield > 0 && g.tail_r >= 0 &&
         g.tau_over_sigma > 0 && g.range_hi > g.range_lo;
  return g;
}

// A channel's seed for RunBkgGeSim, from the state its fit starts from. The
// loaded state overrides most of it; what it sets is the tail the 95.902 keV
// line is measured with for the escapes (the background channel's last peak),
// and anything the state lacks (an old state has no signal-run Pb centroids,
// which then start from the background run's).
FitResult ChannelSeed(const FitState &s, const TString &channel) {
  const Int_t n_peaks = (channel == "sig") ? 3 : 2;
  const Double_t tail_r = s.Get("bkg:LowExpTailAmplitude1", 0.0);
  const Double_t tau = s.Get("bkg:LowExpTailRatio1", 1.0);
  FitResult seed;
  seed.peaks.resize(n_peaks);
  for (Int_t i = 1; i <= n_peaks; i++) {
    PeakFitResult &peak = seed.peaks[i - 1];
    const Bool_t is_ge = (i == 3);
    peak.mu = s.Get(Indexed(channel, "Mu", i).Data(),
                    s.Get(Indexed("bkg", "Mu", i).Data()));
    peak.sigma = is_ge ? s.Get("sig:Sigma3") : PbSigma(s, i);
    peak.gaus_amplitude = s.Get(Indexed(channel, "GausAmplitude", i).Data());
    peak.low_exp_tail_amplitude = tail_r * peak.gaus_amplitude;
    peak.low_exp_tail_ratio = tau;
  }
  seed.bkg_constant = s.Get((channel + ":BkgConstant").Data(), 0.0);
  seed.lin_bkg_slope = s.Get((channel + ":BkgSlope").Data(), 0.0);
  seed.valid = kTRUE;
  return seed;
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

// Each event a Poisson(mean) number of times: mean 1 resamples a run, mean
// 1/2 draws half of one. Two halves drawn this way are independent, and
// their difference averages to zero over trials; a fixed split would leave
// its own difference in every trial's rate-subtracted residual.
std::vector<Double_t> PoissonResample(const std::vector<Double_t> &events,
                                      Double_t mean, TRandom3 &rng) {
  std::vector<Double_t> out;
  out.reserve((size_t)(events.size() * (mean + 0.1)));
  for (size_t i = 0; i < events.size(); i++) {
    Int_t copies = rng.Poisson(mean);
    for (Int_t k = 0; k < copies; k++)
      out.push_back(events[i]);
  }
  return out;
}

// Where an injected line's events are generated: beyond both methods'
// windows (in-situ 65-80 keV, rate-sub 62-71.5 keV), so each sees the line's
// whole shape, its long low tail included.
const Double_t INJECTION_LO_KEV = 55.0, INJECTION_HI_KEV = 82.0;

// The Ge line as the reference fit has it, Gaussian plus the shared low
// exponential tail, with a Poisson number of events whose mean puts
// n_expected in the in-situ window. n_in_window returns how many landed
// there.
std::vector<Double_t> SampleModel(Double_t n_expected, Double_t e_inj,
                                  const GeShape &shape, TRandom3 &rng,
                                  Int_t &n_in_window) {
  std::vector<Double_t> out;
  n_in_window = 0;
  if (!(n_expected > 0))
    return out;
  const Double_t lo = INJECTION_LO_KEV, hi = INJECTION_HI_KEV;
  const Double_t sigma_ge = shape.sigma, rho = shape.tau_over_sigma;
  // Per unit whole-line Gaussian yield: the Gaussian and the tail in each
  // window.
  Double_t gaus_wide = GaussianWindowShare(lo, hi, e_inj, sigma_ge);
  Double_t tail_wide = gaus_wide * shape.tail_r *
                       LowExpTailWindowFactor(lo, hi, e_inj, sigma_ge, rho);
  Double_t gaus_window =
      GaussianWindowShare(shape.range_lo, shape.range_hi, e_inj, sigma_ge);
  Double_t tail_window = gaus_window * InWindowExpRatio(shape, e_inj, sigma_ge);
  Int_t n = rng.Poisson(n_expected * (gaus_wide + tail_wide) /
                        (gaus_window + tail_window));
  if (n <= 0)
    return out;
  RooRealVar x("x", "x", lo, hi);
  RooRealVar mu("mu", "mu", e_inj);
  RooRealVar sigma("sigma", "sigma", sigma_ge);
  RooRealVar tau("tau", "tau", rho);
  mu.setConstant(kTRUE);
  sigma.setConstant(kTRUE);
  tau.setConstant(kTRUE);
  RooGaussian gauss("gauss", "gauss", x, mu, sigma);
  RooLowExpTail lowexp("lowexp", "lowexp", x, mu, sigma, tau);
  RooConstVar f_gaus("f_gaus", "f_gaus", gaus_wide / (gaus_wide + tail_wide));
  RooAddPdf total("total", "total", RooArgList(gauss, lowexp),
                  RooArgList(f_gaus));
  RooDataSet *data = total.generate(RooArgSet(x), n);
  out.reserve(n);
  for (Int_t i = 0; i < data->numEntries(); i++) {
    Double_t e = data->get(i)->getRealValue("x");
    out.push_back(e);
    if (e >= shape.range_lo && e < shape.range_hi)
      n_in_window++;
  }
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

// Production's rate subtraction on event lists: rate histograms one comb period
// wide at N_RATESUB_PHASES starting points, the residual's 95.902 keV line for
// its escapes, the Ge tail fixed, through the same FitRateSubPhases the
// measurement uses. The first starting point's fit is drawn with the trial
// plots.
RateSubResult FitRateSubtractedEvents(const std::vector<Double_t> &sig_events,
                                      const std::vector<Double_t> &bkg_events,
                                      Double_t live_sig, Double_t live_bkg,
                                      const RateSubShape &tail_shape,
                                      Double_t comb_bin, const TString &tag) {
  RateSubShape shape = tail_shape;
  if (EscapeLinesRequested()) {
    std::vector<TH1F *> capture_sigs =
        RateHistsFromEvents(sig_events, live_sig, "inj_capture_sig_" + tag,
                            88.0, 104.0, comb_bin, N_RATESUB_PHASES);
    std::vector<TH1F *> capture_bkgs =
        RateHistsFromEvents(bkg_events, live_bkg, "inj_capture_bkg_" + tag,
                            88.0, 104.0, comb_bin, N_RATESUB_PHASES);
    SetResidualCaptureLine(capture_sigs, capture_bkgs, tag, shape);
    for (size_t k = 0; k < capture_sigs.size(); k++)
      delete capture_sigs[k];
    for (size_t k = 0; k < capture_bkgs.size(); k++)
      delete capture_bkgs[k];
  }
  std::vector<TH1F *> sigs = RateHistsFromEvents(
      sig_events, live_sig, "inj_rate_sig_" + tag, RATESUB_LO - comb_bin,
      RATESUB_HI + comb_bin, comb_bin, N_RATESUB_PHASES);
  std::vector<TH1F *> bkgs = RateHistsFromEvents(
      bkg_events, live_bkg, "inj_rate_bkg_" + tag, RATESUB_LO - comb_bin,
      RATESUB_HI + comb_bin, comb_bin, N_RATESUB_PHASES);
  TF1 *model = nullptr;
  TH1F *res = nullptr;
  RateSubResult out = FitRateSubPhases(sigs, bkgs, shape, tag, &model, &res);
  for (size_t k = 0; k < sigs.size(); k++)
    delete sigs[k];
  for (size_t k = 0; k < bkgs.size(); k++)
    delete bkgs[k];
  if (model && res && INJ_KEEP_TRIAL_PLOTS) {
    const Double_t lo = RATESUB_LO, hi = RATESUB_HI;
    res->GetXaxis()->SetRangeUser(lo, hi);
    res->GetYaxis()->SetTitle(
        Form("Counts/s / %.0f eV", 1000.0 * res->GetBinWidth(1)));
    const Int_t kNPts = 400;
    Double_t xstep = (hi - lo) / (kNPts - 1);
    TGraph *total_graph = new TGraph(kNPts);
    TGraph *bkg_graph = new TGraph(kNPts);
    for (Int_t i = 0; i < kNPts; i++) {
      Double_t xv = lo + i * xstep;
      total_graph->SetPoint(i, xv, model->Eval(xv));
      bkg_graph->SetPoint(i, xv,
                          model->GetParameter(5) +
                              model->GetParameter(6) * (xv - 68.75));
    }
    bkg_graph->SetLineStyle(2);
    std::vector<TGraph *> components = {bkg_graph};
    PlottingUtils::PlotFitWithResiduals(res, total_graph, components, lo, hi,
                                        "rate_fit_" + tag, "fits", "", kFALSE);
    delete total_graph;
    delete bkg_graph;
  }
  delete model;
  delete res;
  return out;
}

// The rate-sub fit's fixed Ge tail: the reference fit's shared tail, as a
// height ratio for the TF1 model.
RateSubShape TestbedRateSubShape(const GeShape &g) {
  PeakFitResult tail_source;
  tail_source.gaus_amplitude = 1.0;
  tail_source.low_exp_tail_amplitude = g.tail_r;
  tail_source.low_exp_tail_ratio = g.tau_over_sigma;
  RateSubShape shape;
  shape.ge_tail_height_ratio = TailHeightRatio(tail_source);
  shape.ge_tau_over_sigma = g.tau_over_sigma;
  return shape;
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

// One fit's Ge centroid and its quality.
struct FitOutcome {
  Double_t mu = -1, err = -1;
  Int_t status = -999, cov_qual = -999;
  Double_t edm = -1;
  Bool_t valid = kFALSE, rail = kFALSE;
};

struct TrialOutcome {
  Int_t n_inj = 0;
  Double_t jitter_eV = 0;
  FitOutcome sim;       // the in-situ fit as production runs it
  FitOutcome sim_noesc; // the same, without the escape lines
  FitOutcome ratesub;   // production's rate subtraction
  Bool_t have_sim_noesc = kFALSE, have_ratesub = kFALSE;
};

// A trial's PNGs, and with them its fit states, are removed after each trial
// unless INJ_KEEP_TRIAL_PLOTS is set.
void CleanTrialFiles(const TString &dir, const TString &trial_tag) {
  TSystemDirectory directory(dir, dir);
  TList *files = directory.GetListOfFiles();
  if (!files)
    return;
  TIter next(files);
  while (TSystemFile *file = static_cast<TSystemFile *>(next())) {
    TString name = file->GetName();
    Bool_t trial_file = name.EndsWith(".png") || name.EndsWith(".simroofits");
    Bool_t this_trial =
        name.Contains(trial_tag + "_") || name.Contains(trial_tag + ".");
    if (trial_file && this_trial)
      gSystem->Unlink(dir + "/" + name);
  }
}

// ------------------------------------------------------- the 95.902 keV line

// Where CaptureLineModel's linear background is centred.
const Double_t CAPTURE_MODEL_X0 = 96.0;

// A run's 95.902 keV line: its Gaussian yield, centroid and width with the
// shared tail, and the continuum under it, c + s (x - CAPTURE_MODEL_X0) per
// keV, from production's capture-line model over its fixed window.
struct CaptureTruth {
  Bool_t ok = kFALSE;
  Double_t gaus_yield = 0, mu = 0, sigma = 0;
  Double_t continuum_per_kev = 0, continuum_slope = 0;
};

CaptureTruth MeasureCaptureTruth(const std::vector<Double_t> &events,
                                 const GeShape &g, Double_t comb_bin) {
  CaptureTruth truth;
  const Int_t n_bins = (Int_t)std::lround(16.0 / comb_bin);
  TH1D h("capture_truth", "", n_bins, 88.0, 88.0 + n_bins * comb_bin);
  h.SetDirectory(nullptr);
  for (size_t i = 0; i < events.size(); i++)
    h.Fill(events[i]);
  PeakFitResult tail_source;
  tail_source.gaus_amplitude = 1.0;
  tail_source.low_exp_tail_amplitude = g.tail_r;
  tail_source.low_exp_tail_ratio = g.tau_over_sigma;
  Double_t tail_height = TailHeightRatio(tail_source);
  Int_t first = h.FindBin(CAPTURE_FIT_LO), last = h.FindBin(CAPTURE_FIT_HI);
  Int_t peak_bin = first;
  for (Int_t b = first; b <= last; b++)
    if (h.GetBinContent(b) > h.GetBinContent(peak_bin))
      peak_bin = b;
  Double_t level =
      0.5 * (h.GetBinContent(first + 1) + h.GetBinContent(last - 1));
  Double_t yield0 = (h.GetBinContent(peak_bin) - level) * 0.9 *
                    TMath::Sqrt(2.0 * TMath::Pi()) / comb_bin;
  TF1 f("capture_truth_fit", CaptureLineModel, CAPTURE_FIT_LO, CAPTURE_FIT_HI,
        8);
  f.SetParameters(yield0, h.GetBinCenter(peak_bin), 0.9, tail_height,
                  g.tau_over_sigma, level, 0.0, comb_bin);
  f.FixParameter(3, tail_height);
  f.FixParameter(4, g.tau_over_sigma);
  f.FixParameter(7, comb_bin);
  f.SetParLimits(2, 0.4, 1.5);
  TFitResultPtr fr = h.Fit(&f, "S L Q R N");
  if (fr.Get() && fr->IsValid() && f.GetParameter(0) > 0) {
    truth.gaus_yield = f.GetParameter(0);
    truth.mu = f.GetParameter(1);
    truth.sigma = f.GetParameter(2);
    truth.continuum_per_kev = f.GetParameter(5) / comb_bin;
    truth.continuum_slope = f.GetParameter(6) / comb_bin;
    truth.ok = kTRUE;
  }
  return truth;
}

// --------------------------------------------------------- the reference fit

// What a pair's trials take from its real runs: the reference fit's state,
// the Ge lineshape, and each run's 95.902 keV line.
struct PairReference {
  Bool_t ok = kFALSE;
  FitState state;
  GeShape ge;
  CaptureTruth bkg_capture, sig_capture;
};

FitOutcome OutcomeOf(const BkgGeSimResult &sim, Double_t range_lo,
                     Double_t range_hi) {
  FitOutcome out;
  if (sim.sig_channel.peaks.empty())
    return out;
  const PeakFitResult &ge = sim.sig_channel.peaks.back();
  out.mu = ge.mu;
  out.err = ge.mu_error;
  out.status = sim.sig_channel.fit_status;
  out.cov_qual = sim.sig_channel.cov_qual;
  out.edm = sim.sig_channel.edm;
  out.valid = sim.valid && std::isfinite(ge.mu) && ge.mu_error > 0 &&
              sim.sig_channel.has_fit_diagnostics;
  Double_t margin = 0.01 * (range_hi - range_lo);
  out.rail = (ge.mu <= range_lo + margin || ge.mu >= range_hi - margin);
  return out;
}

// Production's in-situ fit from the state `start`, which the caller has put
// where the fitter looks for it (<plots dir>/fits/<fit_label>_<sig_label>
// .simroofits).
BkgGeSimResult RunInSituFromState(const FitState &start,
                                  const std::vector<Double_t> &bkg_events,
                                  const std::vector<Double_t> &sig_events,
                                  const TString &bkg_label,
                                  const TString &sig_label,
                                  const TString &fit_label,
                                  const InSituEscapes &escapes) {
  FitResult bkg_seed = ChannelSeed(start, "bkg");
  FitResult sig_seed = ChannelSeed(start, "sig");
  return RunBkgGeSim(bkg_events, sig_events, PbKaEnergies(), bkg_label,
                     sig_label, fit_label, &bkg_seed, &sig_seed, escapes);
}

Bool_t CopyFile(const TString &source, const TString &target) {
  std::ifstream in(source.Data(), std::ios::binary);
  std::ofstream out(target.Data(), std::ios::binary);
  if (!in.is_open() || !out.is_open())
    return kFALSE;
  out << in.rdbuf();
  return kTRUE;
}

// Production's post-cal fit of the pair, on the real runs, from production's
// saved state, which is copied and never written. As production does, the
// 95.902 keV line is measured with a converged shared tail: here the same fit
// without the escapes, which production's pre-cal fit stands in for.
PairReference RunReferenceFit(const InjPair &pair,
                              const TString &production_state,
                              const TString &reference_dir,
                              const std::vector<Double_t> &bkg_all,
                              const std::vector<Double_t> &sig_all,
                              Double_t comb_bin) {
  PairReference ref;
  FitState production = ReadFitState(production_state);
  if (!production.ok) {
    std::cerr << "INJ-FATAL: missing or incomplete post-cal state "
              << production_state << std::endl;
    return ref;
  }
  const TString fits_dir = reference_dir + "/fits";
  gSystem->mkdir(fits_dir, kTRUE);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, reference_dir,
                                ProjectRoot() + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  const TString bkg_label = pair.bkg_run + "_postcal_" + pair.tag;
  const TString sig_label = pair.sig_run + "_postcal_" + pair.tag;

  const TString first_label = "PostCalNoEsc_" + pair.tag;
  if (!CopyFile(production_state,
                fits_dir + "/" + first_label + "_" + sig_label + ".simroofits"))
    return ref;
  BkgGeSimResult first =
      RunInSituFromState(production, bkg_all, sig_all, bkg_label, sig_label,
                         first_label, InSituEscapes());
  FitState converged =
      FitStateOf(first, production.range_lo, production.range_hi);
  BkgGeSimResult fit = first;
  if (EscapeLinesRequested() && converged.ok) {
    InSituEscapes escapes;
    escapes.fractions = pair.escapes;
    escapes.units_per_kev = 1.0;
    escapes.comb_period = comb_bin;
    const TString label = "PostCal_" + pair.tag;
    if (!WriteFitState(converged, fits_dir + "/" + label + "_" + sig_label +
                                      ".simroofits"))
      return ref;
    fit = RunInSituFromState(converged, bkg_all, sig_all, bkg_label, sig_label,
                             label, escapes);
  }
  ref.state = FitStateOf(fit, production.range_lo, production.range_hi);
  ref.ge = GeShapeOf(ref.state);
  if (!ref.state.ok || !ref.ge.ok) {
    std::cerr << "INJ-FATAL: the reference fit of " << pair.label
              << " did not converge" << std::endl;
    return ref;
  }
  ref.bkg_capture = MeasureCaptureTruth(bkg_all, ref.ge, comb_bin);
  ref.sig_capture = MeasureCaptureTruth(sig_all, ref.ge, comb_bin);

  FitOutcome ge = OutcomeOf(fit, production.range_lo, production.range_hi);
  std::cout << "INJ-REFERENCE " << pair.label << ": Ge " << std::fixed
            << std::setprecision(4) << ge.mu << " +/- " << ge.err
            << " keV (status " << ge.status << ", chi2/ndf bkg "
            << std::setprecision(3) << fit.bkg_channel.reduced_chi2 << ", sig "
            << fit.sig_channel.reduced_chi2 << "), sigma "
            << std::setprecision(4) << ref.ge.sigma << ", tail R "
            << ref.ge.tail_r << " tau/sigma " << ref.ge.tau_over_sigma
            << ", Ge area " << std::setprecision(0) << ref.ge.area << std::endl;
  const CaptureTruth *captures[2] = {&ref.bkg_capture, &ref.sig_capture};
  const char *runs[2] = {"background", "signal"};
  for (Int_t c = 0; c < 2; c++)
    std::cout << "INJ-CAPTURE " << pair.label << " " << runs[c]
              << " run: " << (captures[c]->ok ? "" : "NOT FOUND, ")
              << std::fixed << std::setprecision(4) << captures[c]->mu
              << " keV, sigma " << captures[c]->sigma << ", Gaussian yield "
              << std::setprecision(0) << captures[c]->gaus_yield
              << ", continuum " << captures[c]->continuum_per_kev
              << " per keV, slope " << std::setprecision(1)
              << captures[c]->continuum_slope << " per keV^2" << std::endl;
  ref.ok = kTRUE;
  return ref;
}

// --------------------------------------------------------------- synthetic

// Both channels generated from the reference fit, each with its own run's
// content at that run's statistics: the two Pb lines, the flat continuum and,
// in the signal channel, the Ge line at E_inj, all with the shared tail.
// With escapes, each channel also carries its run's 95.902 keV line on its
// continuum between 85 and 104 keV and, in the window, the line's eight K
// X-ray escapes at the escape fractions; production's fit then has a line to
// measure and escapes to remove, and the same events fitted without them
// show what ignoring them costs.

struct SyntheticData {
  std::vector<Double_t> bkg, sig;
};

RooRealVar *MakeOwnedVar(std::vector<RooAbsArg *> &owned, const TString &name,
                         Double_t value) {
  RooRealVar *v = new RooRealVar(name, name, value);
  owned.push_back(v);
  return v;
}

RooConstVar *MakeOwnedConst(std::vector<RooAbsArg *> &owned,
                            const TString &name, Double_t value) {
  RooConstVar *v = new RooConstVar(name, name, value);
  owned.push_back(v);
  return v;
}

// One line of the fit's family for a generator over x's range: of a whole-line
// Gaussian yield, the share inside that range, and the shared tail's matching
// yield there.
void AddTailedLine(std::vector<RooAbsArg *> &owned, RooArgList &pdfs,
                   RooArgList &coefs, RooRealVar &x, const TString &name,
                   Double_t mu, Double_t sigma, Double_t whole_line_yield,
                   const GeShape &shape) {
  const Double_t lo = x.getMin(), hi = x.getMax();
  const Double_t gaus_yield =
      whole_line_yield * GaussianWindowShare(lo, hi, mu, sigma);
  const Double_t tail_yield =
      gaus_yield * shape.tail_r *
      LowExpTailWindowFactor(lo, hi, mu, sigma, shape.tau_over_sigma);
  RooRealVar *m = MakeOwnedVar(owned, "mu_" + name, mu);
  RooRealVar *s = MakeOwnedVar(owned, "sigma_" + name, sigma);
  RooRealVar *t = MakeOwnedVar(owned, "tau_" + name, shape.tau_over_sigma);
  RooGaussian *g = new RooGaussian("g_" + name, "", x, *m, *s);
  RooLowExpTail *e = new RooLowExpTail("e_" + name, "", x, *m, *s, *t);
  owned.push_back(g);
  owned.push_back(e);
  pdfs.add(*g);
  coefs.add(*MakeOwnedConst(owned, "yg_" + name, gaus_yield));
  pdfs.add(*e);
  coefs.add(*MakeOwnedConst(owned, "ye_" + name, tail_yield));
}

// A continuum c + s (x - x0) per keV over x's range, c > 0.
void AddLinearContinuum(std::vector<RooAbsArg *> &owned, RooArgList &pdfs,
                        RooArgList &coefs, RooRealVar &x, const TString &name,
                        Double_t c, Double_t s, Double_t x0) {
  const Double_t lo = x.getMin(), hi = x.getMax();
  // RooPolynomial's 1 + a x, normalised over the range.
  RooRealVar *a = MakeOwnedVar(owned, "a_" + name, s / (c - s * x0));
  RooPolynomial *poly = new RooPolynomial(name, "", x, RooArgList(*a));
  owned.push_back(poly);
  pdfs.add(*poly);
  Double_t yield =
      c * (hi - lo) + 0.5 * s * ((hi - x0) * (hi - x0) - (lo - x0) * (lo - x0));
  coefs.add(*MakeOwnedConst(owned, "y_" + name, yield));
}

std::vector<Double_t> GenerateFromComponents(RooRealVar &x,
                                             const RooArgList &pdfs,
                                             const RooArgList &coefs) {
  RooAddPdf pdf("synthetic_pdf", "", pdfs, coefs);
  RooDataSet *data = pdf.generate(RooArgSet(x), RooFit::Extended(kTRUE));
  std::vector<Double_t> events;
  events.reserve(data->numEntries());
  for (Int_t i = 0; i < data->numEntries(); i++)
    events.push_back(data->get(i)->getRealValue(x.GetName()));
  delete data;
  return events;
}

void DeleteOwned(std::vector<RooAbsArg *> &owned) {
  for (size_t i = owned.size(); i > 0; i--)
    delete owned[i - 1];
  owned.clear();
}

// The synthetic window reaches this far beyond the fit window on each side,
// as real data do: cut at the fit window's edges, the edge bins of the fit's
// binned chi2 would be half empty.
const Double_t SYNTHETIC_MARGIN_KEV = 1.0;

// One channel's in-situ window; n_inj 0 leaves the Ge line out, escapes null
// the escapes. The state's yields, and n_inj, are the fit window's.
std::vector<Double_t> GenerateSyntheticWindow(const PairReference &ref,
                                              const TString &channel,
                                              Double_t e_inj, Double_t n_inj,
                                              const CaptureTruth &capture,
                                              const EscapeFractions *escapes) {
  const FitState &s = ref.state;
  const GeShape &shape = ref.ge;
  RooRealVar x("x", "x", s.range_lo - SYNTHETIC_MARGIN_KEV,
               s.range_hi + SYNTHETIC_MARGIN_KEV);
  std::vector<RooAbsArg *> owned;
  RooArgList pdfs, coefs;
  for (Int_t i = 1; i <= 2; i++) {
    Double_t mu = s.Get(Indexed(channel, "Mu", i).Data());
    Double_t sigma = PbSigma(s, i);
    Double_t yield = s.Get(Indexed(channel, "GausAmplitude", i).Data());
    if (mu > 0 && yield > 0)
      AddTailedLine(
          owned, pdfs, coefs, x, TString::Format("pb%d", i), mu, sigma,
          yield / GaussianWindowShare(s.range_lo, s.range_hi, mu, sigma),
          shape);
  }
  // The escapes with the Pb Ka1 shape the fit draws them with, and their
  // whole-line Gaussian yields.
  if (escapes && capture.ok) {
    Double_t sigma = PbSigma(s, 2);
    for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
      Double_t whole_line = capture.gaus_yield *
                            escapes->series[ESCAPE_SERIES[i]] *
                            ESCAPE_SERIES_SHARE[i];
      if (whole_line > 0)
        AddTailedLine(owned, pdfs, coefs, x, TString::Format("esc%d", i + 1),
                      capture.mu - ESCAPE_XRAY_KEV[i], sigma, whole_line,
                      shape);
    }
  }
  if (n_inj > 0) {
    Double_t gaus_in_window =
        n_inj / (1.0 + InWindowExpRatio(shape, e_inj, shape.sigma));
    AddTailedLine(owned, pdfs, coefs, x, "ge", e_inj, shape.sigma,
                  gaus_in_window / GaussianWindowShare(s.range_lo, s.range_hi,
                                                       e_inj, shape.sigma),
                  shape);
  }
  Double_t continuum = s.Get((channel + ":BkgConstant").Data(), 0.0);
  if (continuum > 0)
    AddLinearContinuum(owned, pdfs, coefs, x, "bkg_flat",
                       continuum / (s.range_hi - s.range_lo), 0.0, 0.0);
  std::vector<Double_t> events = GenerateFromComponents(x, pdfs, coefs);
  DeleteOwned(owned);
  return events;
}

// A run's 95.902 keV line with the shared tail, on its continuum, between 85
// and 104 keV: what production measures the escapes' yield and position
// from.
std::vector<Double_t> GenerateCaptureRegion(const GeShape &shape,
                                            const CaptureTruth &capture) {
  RooRealVar x("x_capture", "x", 85.0, 104.0);
  std::vector<RooAbsArg *> owned;
  RooArgList pdfs, coefs;
  AddTailedLine(owned, pdfs, coefs, x, "capture", capture.mu, capture.sigma,
                capture.gaus_yield, shape);
  AddLinearContinuum(owned, pdfs, coefs, x, "capture_continuum",
                     capture.continuum_per_kev, capture.continuum_slope,
                     CAPTURE_MODEL_X0);
  std::vector<Double_t> events = GenerateFromComponents(x, pdfs, coefs);
  DeleteOwned(owned);
  return events;
}

SyntheticData GenerateSynthetic(const PairReference &ref, Double_t e_inj,
                                Double_t n_inj,
                                const EscapeFractions *escapes) {
  SyntheticData out;
  out.bkg =
      GenerateSyntheticWindow(ref, "bkg", e_inj, 0.0, ref.bkg_capture, escapes);
  out.sig = GenerateSyntheticWindow(ref, "sig", e_inj, n_inj, ref.sig_capture,
                                    escapes);
  if (escapes) {
    if (ref.bkg_capture.ok) {
      std::vector<Double_t> region =
          GenerateCaptureRegion(ref.ge, ref.bkg_capture);
      out.bkg.insert(out.bkg.end(), region.begin(), region.end());
    }
    if (ref.sig_capture.ok) {
      std::vector<Double_t> region =
          GenerateCaptureRegion(ref.ge, ref.sig_capture);
      out.sig.insert(out.sig.end(), region.begin(), region.end());
    }
  }
  return out;
}

// A synthetic background run at the background run's statistics: the
// reference fit's background channel and its 95.902 keV region, as the
// synthetic source makes them, without a Ge line.
std::vector<Double_t> GenerateBackgroundRun(const PairReference &ref,
                                            const EscapeFractions *escapes) {
  const EscapeFractions *with = INJ_SYNTHETIC_ESCAPES ? escapes : nullptr;
  std::vector<Double_t> events =
      GenerateSyntheticWindow(ref, "bkg", 0.0, 0.0, ref.bkg_capture, with);
  if (with && ref.bkg_capture.ok) {
    std::vector<Double_t> region =
        GenerateCaptureRegion(ref.ge, ref.bkg_capture);
    events.insert(events.end(), region.begin(), region.end());
  }
  return events;
}

// Halves trials draw both channels from one pool; synthetic sources lack
// the 62-65 keV the rate subtraction needs and use the model lineshape.
Bool_t HalvesLike(const TString &source) {
  return source == "halves" || source == "synthetic_halves";
}
Bool_t SyntheticSource(const TString &source) {
  return source == "synthetic" || source == "synthetic_halves";
}

// ------------------------------------------------------------------ trials

// What a pair's trials share.
struct PairSetup {
  InjPair pair;
  TString source;
  PairReference ref;
  // What halves (both the background run) and real trials resample.
  std::vector<Double_t> bkg_pool, sig_pool;
  std::vector<Double_t> am_events;     // the empirical lineshape
  Double_t live_sig = 0, live_bkg = 0; // the channels' live times, for Method 2
  Double_t n_inj_expected = 0;
  Double_t comb_bin = 0;
  RateSubShape rate_sub_shape;
};

// The reference state scaled to what the trial's channels hold, with the Ge
// line started at ge_start.
FitState StartState(const PairSetup &setup, Double_t ge_start) {
  FitState start = setup.ref.state;
  if (HalvesLike(setup.source)) {
    ScaleYield(start, "bkg:GausAmplitude1", 0.5);
    ScaleYield(start, "bkg:GausAmplitude2", 0.5);
    ScaleYield(start, "bkg:BkgConstant", 0.5);
    const char *shared[5] = {"Mu1", "Mu2", "GausAmplitude1", "GausAmplitude2",
                             "BkgConstant"};
    for (Int_t i = 0; i < 5; i++) {
      std::string from = std::string("bkg:") + shared[i];
      start.Set(std::string("sig:") + shared[i], start.Get(from),
                start.Error(from));
    }
    if (start.Has("bkg:BkgSlope"))
      start.Set("sig:BkgSlope", start.Get("bkg:BkgSlope"),
                start.Error("bkg:BkgSlope"));
  }
  const GeShape &ge = setup.ref.ge;
  Double_t ge_gaus =
      setup.n_inj_expected / (1.0 + InWindowExpRatio(ge, ge_start, ge.sigma));
  if (setup.source == "real")
    ge_gaus += setup.ref.state.Get("sig:GausAmplitude3");
  start.Set("sig:Mu3", ge_start, setup.ref.state.Error("sig:Mu3"));
  start.Set("sig:GausAmplitude3", ge_gaus, TMath::Sqrt(ge_gaus));
  return start;
}

TrialOutcome RunOneTrial(const PairSetup &setup, const TString &shape,
                         Double_t e_inj, UInt_t seed, const TString &fits_dir,
                         const TString &trial_tag) {
  TrialOutcome out;
  TRandom3 rng(seed);
  Double_t jitter_kev = rng.Uniform(-0.050, 0.050);
  out.jitter_eV = jitter_kev * 1000.0;

  std::vector<Double_t> bkg_events, sig_events;
  if (setup.source == "synthetic") {
    out.n_inj = (Int_t)std::lround(setup.n_inj_expected);
    SyntheticData syn =
        GenerateSynthetic(setup.ref, e_inj, setup.n_inj_expected,
                          INJ_SYNTHETIC_ESCAPES ? setup.pair.escapes : nullptr);
    bkg_events.swap(syn.bkg);
    sig_events.swap(syn.sig);
  } else {
    const Double_t share = HalvesLike(setup.source) ? 0.5 : 1.0;
    bkg_events = PoissonResample(setup.bkg_pool, share, rng);
    sig_events = PoissonResample(setup.sig_pool, share, rng);
    std::vector<Double_t> injected;
    if (shape == "model") {
      injected = SampleModel(setup.n_inj_expected, e_inj, setup.ref.ge, rng,
                             out.n_inj);
    } else {
      out.n_inj =
          setup.n_inj_expected > 0 ? rng.Poisson(setup.n_inj_expected) : 0;
      injected = SampleEmpirical(out.n_inj, e_inj, setup.ref.ge.sigma,
                                 setup.am_events);
    }
    sig_events.insert(sig_events.end(), injected.begin(), injected.end());
  }

  const FitState start = StartState(setup, e_inj + jitter_kev);
  InSituEscapes with_escapes;
  with_escapes.fractions =
      EscapeLinesRequested() ? setup.pair.escapes : nullptr;
  with_escapes.units_per_kev = 1.0;
  with_escapes.comb_period = setup.comb_bin;
  const Bool_t fit_without = INJ_FIT_WITHOUT_ESCAPES && with_escapes.fractions;
  const TString sig_label = "Sig_" + trial_tag;
  std::vector<TString> fit_labels;
  fit_labels.push_back("Inj_" + trial_tag);
  if (fit_without)
    fit_labels.push_back("InjNoEsc_" + trial_tag);
  for (size_t f = 0; f < fit_labels.size(); f++) {
    if (!WriteFitState(start, fits_dir + "/" + fit_labels[f] + "_" + sig_label +
                                  ".simroofits")) {
      std::cerr << "INJ-FATAL: could not write the trial state for "
                << fit_labels[f] << std::endl;
      return out;
    }
  }

  const Double_t lo = start.range_lo, hi = start.range_hi;
  out.sim = OutcomeOf(RunInSituFromState(start, bkg_events, sig_events,
                                         "Bkg_" + trial_tag, sig_label,
                                         fit_labels[0], with_escapes),
                      lo, hi);
  if (fit_without) {
    out.sim_noesc = OutcomeOf(
        RunInSituFromState(start, bkg_events, sig_events, "Bkg_" + trial_tag,
                           sig_label, fit_labels[1], InSituEscapes()),
        lo, hi);
    out.have_sim_noesc = kTRUE;
  }

  // The synthetic channels hold only the in-situ window and the 95.902 keV
  // region, not the 62-65 keV the rate subtraction also fits.
  if (!SyntheticSource(setup.source)) {
    RateSubResult rs = FitRateSubtractedEvents(
        sig_events, bkg_events, setup.live_sig, setup.live_bkg,
        setup.rate_sub_shape, setup.comb_bin, trial_tag + "_ratesub");
    out.ratesub.mu = rs.mu;
    out.ratesub.err = rs.mu_err;
    out.ratesub.valid = rs.valid && std::isfinite(rs.mu) && rs.mu_err > 0;
    out.ratesub.rail =
        (rs.mu <= RATESUB_LO + 1.01 || rs.mu >= RATESUB_HI - 1.01);
    out.have_ratesub = kTRUE;
  }

  if (!INJ_KEEP_TRIAL_PLOTS) {
    CleanTrialFiles(fits_dir, trial_tag);
    CleanTrialFiles(fits_dir + "/residual_hists", trial_tag);
  }
  return out;
}

TString ShapeTag(const TString &shape) {
  return shape == "model" ? "model" : "emp";
}

TString CellTag(const InjPair &pair, const TString &shape, Double_t e_inj) {
  TString e_tag = Form("E%.3fkeV", e_inj);
  e_tag.ReplaceAll(".", "_");
  return pair.tag + "_" + ShapeTag(shape) + "_" + e_tag;
}

void WriteRow(std::ofstream &csv, const InjPair &pair, const TString &shape,
              Double_t e_inj, Int_t trial, Int_t n_inj, Double_t jitter_eV,
              const TString &method, const FitOutcome &fit) {
  Double_t residual_eV =
      (fit.mu > 0 && std::isfinite(fit.mu)) ? (fit.mu - e_inj) * 1000.0 : 0.0;
  csv << pair.tag.Data() << "," << shape.Data() << "," << std::fixed
      << std::setprecision(4) << e_inj << "," << trial << "," << n_inj << ","
      << std::setprecision(3) << jitter_eV << "," << method.Data() << ","
      << std::setprecision(8) << fit.mu << "," << fit.err << ","
      << std::setprecision(4) << residual_eV << "," << fit.status << ","
      << fit.cov_qual << "," << std::setprecision(6) << fit.edm << ","
      << (fit.valid ? 1 : 0) << "," << (fit.rail ? 1 : 0) << std::endl;
}

} // namespace

void InjectionTest() {
  const TString &source = INJ_SOURCE;
  const Int_t pair_index = INJ_PAIR_INDEX;
  const Int_t energy_index = INJ_ENERGY_INDEX;
  const Int_t n_trials = INJ_N_TRIALS;
  const Int_t trial_offset = INJ_TRIAL_OFFSET;
  const TString &run_tag = INJ_RUN_TAG;
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if (source != "halves" && source != "real" && source != "synthetic" &&
      source != "synthetic_halves") {
    std::cerr << "ERROR: source must be halves, real, synthetic or "
                 "synthetic_halves, got '"
              << source << "'" << std::endl;
    return;
  }
  if (INJ_SHAPE != "model" && INJ_SHAPE != "empirical" && INJ_SHAPE != "both") {
    std::cerr << "ERROR: shape must be model, empirical, or both" << std::endl;
    return;
  }
  // Synthetic data are generated from the fit's own model.
  std::vector<TString> shapes;
  if (INJ_SHAPE == "both" && !SyntheticSource(source)) {
    shapes.push_back("model");
    shapes.push_back("empirical");
  } else if (SyntheticSource(source)) {
    shapes.push_back("model");
  } else {
    shapes.push_back(INJ_SHAPE);
  }

  std::cout << "INJ-CONFIG tag " << RESULT_TAG << " source " << source
            << " synthetic escapes " << (INJ_SYNTHETIC_ESCAPES ? 1 : 0)
            << " fit without escapes " << (INJ_FIT_WITHOUT_ESCAPES ? 1 : 0)
            << " shape " << INJ_SHAPE << " pair " << pair_index << " energy "
            << energy_index << " trials " << n_trials << " from "
            << trial_offset << " run " << run_tag << std::endl;
  const TString project_root = ProjectRoot();
  const TString base = project_root + "/plots/testbed/injection_" + run_tag;
  gSystem->mkdir(base + "/results", kTRUE);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, base,
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);

  std::vector<InjPair> pairs = MakePairs();
  for (Int_t pi = 0; pi < (Int_t)pairs.size(); pi++) {
    if (pair_index >= 0 && pi != pair_index)
      continue;
    PairSetup setup;
    setup.pair = pairs[pi];
    setup.source = source;
    const InjPair &pair = setup.pair;

    Double_t slope = PairCalibrationSlope(pair.tag);
    Double_t t_bkg = RunLiveTime(pair.bkg_run);
    Double_t t_sig = RunLiveTime(pair.sig_run);
    if (!(slope > 0) || !(t_bkg > 0 && t_sig > 0)) {
      std::cerr << "INJ-FATAL: no calibration (" << slope << ") or live time ("
                << t_bkg << ", " << t_sig << ") for " << pair.label
                << std::endl;
      return;
    }
    setup.comb_bin = COMB_PERIOD_RAW_KEV * slope;
    std::vector<Double_t> bkg_all = LoadCalibrated(pair.bkg_run);
    std::vector<Double_t> sig_all = LoadCalibrated(pair.sig_run);
    if (bkg_all.empty() || sig_all.empty()) {
      std::cerr << "INJ-FATAL: no calibrated events for " << pair.label
                << std::endl;
      return;
    }

    const TString production_state =
        project_root + "/plots/czt/calibration_low_" + RESULT_TAG +
        "/fits/PostCal_" + pair.tag + "_" + pair.sig_run + "_postcal_" +
        pair.tag + ".simroofits";
    const TString reference_dir =
        base + "/reference/" + pair.tag +
        TString::Format("_e%d_t%d", energy_index, trial_offset);
    setup.ref = RunReferenceFit(pair, production_state, reference_dir, bkg_all,
                                sig_all, setup.comb_bin);
    if (!setup.ref.ok)
      return;
    setup.rate_sub_shape = TestbedRateSubShape(setup.ref.ge);

    // halves: the injected line has the signal run's Ge rate, over half the
    // background run's live time. real: the signal run's own Ge again, so
    // the recovered line is a blend of the two and only the method difference
    // at one energy is meaningful. synthetic: the signal run's Ge.
    if (source == "synthetic_halves") {
      setup.live_bkg = 0.5 * t_bkg;
      setup.live_sig = 0.5 * t_bkg;
      setup.n_inj_expected = setup.ref.ge.area * setup.live_sig / t_sig;
    } else if (source == "halves") {
      setup.bkg_pool = bkg_all;
      setup.sig_pool.swap(bkg_all);
      setup.live_bkg = 0.5 * t_bkg;
      setup.live_sig = 0.5 * t_bkg;
      setup.n_inj_expected = setup.ref.ge.area * setup.live_sig / t_sig;
    } else if (source == "real") {
      setup.bkg_pool.swap(bkg_all);
      setup.sig_pool.swap(sig_all);
      setup.live_bkg = t_bkg;
      setup.live_sig = t_sig;
      setup.n_inj_expected = setup.ref.ge.area;
    } else {
      setup.n_inj_expected = setup.ref.ge.area;
    }
    std::vector<Double_t>().swap(bkg_all);
    std::vector<Double_t>().swap(sig_all);
    if (!SyntheticSource(source) && INJ_SHAPE != "model") {
      setup.am_events = LoadCalibrated(pair.am_run);
      if (setup.am_events.empty()) {
        std::cerr << "INJ-FATAL: no calibrated Am-241 events" << std::endl;
        return;
      }
    }
    std::cout << "INJ-PAIR " << pair.label << " t_bkg " << std::fixed
              << std::setprecision(1) << t_bkg << " t_sig " << t_sig
              << " n_inj_expected " << setup.n_inj_expected << " comb bin "
              << std::setprecision(5) << setup.comb_bin << std::endl;

    const Int_t n_backgrounds =
        (source == "synthetic_halves") ? INJ_N_BACKGROUNDS : 1;
    std::vector<Double_t> background_means;
    for (Int_t bi = 0; bi < n_backgrounds; bi++) {
      std::vector<Double_t> residuals_eV;
      if (source == "synthetic_halves") {
        RooRandom::randomGenerator()->SetSeed(
            TrialSeed(pi, 1000 + bi, 0, 0, run_tag));
        setup.bkg_pool = GenerateBackgroundRun(setup.ref, pair.escapes);
        setup.sig_pool = setup.bkg_pool;
        std::cout << "INJ-BACKGROUND " << pair.label << " b" << bi << ": "
                  << setup.bkg_pool.size() << " synthetic events" << std::endl;
      }
      for (Int_t ei = 0; ei < 9; ei++) {
        if (energy_index >= 0 && ei != energy_index)
          continue;
        Double_t e_inj = kInjEnergies[ei];
        for (Int_t si = 0; si < (Int_t)shapes.size(); si++) {
          const TString &shape = shapes[si];
          TString cell_tag =
              CellTag(pair, shape, e_inj) +
              (source == "halves" ? TString("") : TString("_") + source) +
              (source == "synthetic_halves" ? TString(Form("_b%d", bi))
                                            : TString(""));
          TString cell_dir = base + "/cells/" + cell_tag;
          TString fits_dir = cell_dir + "/fits";
          gSystem->mkdir(fits_dir, kTRUE);
          InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, cell_dir,
                                        project_root + "/root_files");
          PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
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
            TrialOutcome out =
                RunOneTrial(setup, shape, e_inj, seed, fits_dir, trial_tag);
            WriteRow(csv, pair, shape, e_inj, trial, out.n_inj, out.jitter_eV,
                     "sim", out.sim);
            if (out.sim.valid)
              residuals_eV.push_back(1000.0 * (out.sim.mu - e_inj));
            if (out.have_sim_noesc)
              WriteRow(csv, pair, shape, e_inj, trial, out.n_inj, out.jitter_eV,
                       "sim_noesc", out.sim_noesc);
            if (out.have_ratesub)
              WriteRow(csv, pair, shape, e_inj, trial, out.n_inj, 0.0,
                       "ratesub", out.ratesub);
            std::cout << "INJ-TRIAL " << trial_tag << " n_inj " << out.n_inj
                      << " sim " << std::setprecision(5) << out.sim.mu
                      << " +/- " << out.sim.err << " (status " << out.sim.status
                      << ", valid " << (out.sim.valid ? 1 : 0) << ")";
            if (out.have_sim_noesc)
              std::cout << "  without escapes " << out.sim_noesc.mu
                        << " (status " << out.sim_noesc.status << ")";
            if (out.have_ratesub)
              std::cout << "  ratesub " << out.ratesub.mu;
            std::cout << std::endl;
          }
          csv.close();
          std::cout << "INJ-CELL-DONE " << cell_tag << std::endl;
        }
      }
      if (source == "synthetic_halves" && residuals_eV.size() > 1) {
        Double_t mean = 0, var = 0;
        for (size_t k = 0; k < residuals_eV.size(); k++)
          mean += residuals_eV[k];
        mean /= residuals_eV.size();
        for (size_t k = 0; k < residuals_eV.size(); k++)
          var += (residuals_eV[k] - mean) * (residuals_eV[k] - mean);
        var /= (residuals_eV.size() - 1);
        background_means.push_back(mean);
        std::cout << "INJ-BACKGROUND " << pair.label << " b" << bi
                  << ": in-situ residual " << std::fixed << std::setprecision(1)
                  << mean << " +- " << TMath::Sqrt(var / residuals_eV.size())
                  << " eV over " << residuals_eV.size() << " valid trials"
                  << std::endl;
      }
    }
    if (background_means.size() > 1) {
      Double_t mean = 0, var = 0;
      for (size_t k = 0; k < background_means.size(); k++)
        mean += background_means[k];
      mean /= background_means.size();
      for (size_t k = 0; k < background_means.size(); k++)
        var += (background_means[k] - mean) * (background_means[k] - mean);
      var /= (background_means.size() - 1);
      std::cout << "INJ-BACKGROUNDS " << pair.label << ": mean residual "
                << std::fixed << std::setprecision(1) << mean
                << " eV, spread between backgrounds " << TMath::Sqrt(var)
                << " eV over " << background_means.size() << " backgrounds"
                << std::endl;
    }
  }
}
