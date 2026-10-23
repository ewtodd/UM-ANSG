#ifndef GAINTRANSFER_FIT_HPP
#define GAINTRANSFER_FIT_HPP

// Shared fitting core for the 73mGe gain-transfer analysis.
//
// Both CalibrateAndFit.cpp (the production pipeline) and InjectionTest.cpp (the
// closure test) include this, so the injection test exercises the SAME
// simultaneous-fit and rate-subtracted-fit code that produces the measurement.
// Before this existed, InjectionTest.cpp carried its own copy in
// CalibrationLow.cpp and the two could disagree for reasons unrelated to the
// physics.

#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "RooFitPhotopeakPdfs.hpp"
#include "RooFitUtils.hpp"

#include <RooRealVar.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------- constants

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_79 = 79.6142;
const Float_t E_BA133_81 = 80.9979;
const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;

const Float_t E_NEAR_GE = 69;

// Anchor reference-energy uncertainties (keV), folded into the cal-fit x-errors
// so they propagate into the pol1 parameter covariance and thus the Ge energy.
// Pb Ka are known to <1 eV (Deslattes).
//
// dE_LINE remains the conservative ~10 eV lump, now used ONLY for the Ba-133
// lines of the 15th/16th line-cal days.
// TODO: replace with per-line evaluated Ba-133 uncertainties.
const Float_t dE_PB_KA = 0.001;
const Float_t dE_AM241 = 0.0001;
const Float_t dE_LINE = 0.010;

const Float_t CAL_RANGE_LOW = 65;
const Float_t CAL_RANGE_HIGH = 80;
const Float_t RATESUB_LO = 62.0;
const Float_t RATESUB_HI = 71.5;

const Bool_t LINK_PB_KA_SIGMA = kTRUE;

const Bool_t USE_FLAT_BKG = kTRUE;

// LINESHAPE=legacy reproduces the lineshape used before 2026-10-22 (tail
// capped at 8 sigma, Gaussian rate-sub fit) under its own tag, so the effect
// of the lineshape change can be measured on its own.
inline Bool_t LegacyLineshape() {
  const char *env = getenv("LINESHAPE");
  return env && TString(env) == "legacy";
}

// FIT_HARD=1: Minuit strategy 2, a final HESSE and 200k calls in RunBkgGeSim.
inline Bool_t HardFitRequested() {
  const char *env = getenv("FIT_HARD");
  return env && *env && atoi(env) != 0;
}

// PB_ANCHOR=signal: each in-situ signal run is calibrated on its own Pb Ka
// lines instead of inheriting the background run's. The two runs of a pair sit
// at different rates, and the detector reading moves with rate by tens of eV.
inline Bool_t SignalPbAnchorRequested() {
  const char *env = getenv("PB_ANCHOR");
  return env && TString(env) == "signal";
}

// At 8 the low-exp tail railed on every line. Up to 20 it settles inside
// (tau/sigma 9.5 Am, 10.5 Ba, 10.7 Pb in-situ) and the Pb Ka chi2 drops from
// 2.6 to 1.1 per degree of freedom.
const Double_t TAIL_RATIO_MAX = LegacyLineshape() ? 8.0 : 20.0;

// Fixed when the macro loads, from the three switches above, so each
// configuration writes its own calibrated trees, plots and results.
const TString RESULT_TAG = TString(LegacyLineshape() ? "tau8" : "tau20") +
                           (HardFitRequested() ? "_hard" : "") +
                           (SignalPbAnchorRequested() ? "_sigpb" : "");
const TString OUT_SUBDIR = "calibrated_" + RESULT_TAG;

// Which result set this run writes. RESULT_TAG names the default one;
// RESULT_TAG_OVERRIDE lets the same macro produce a second set (e.g. a
// harder-convergence re-run) instead of keeping a duplicate source file.
// It changes OUTPUT paths only -- the calibrated inputs always come from
// OUT_SUBDIR, because overriding the tag must not change which calibrated
// trees are read, only where the results land.
inline TString RuntimeResultTag() {
  const char *env = getenv("RESULT_TAG_OVERRIDE");
  if (env && *env)
    return TString(env);
  return RESULT_TAG;
}

// ------------------------------------------------------------------- events

inline std::vector<Double_t> LoadCalibrated(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading(OUT_SUBDIR + "/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << OUT_SUBDIR << "/" << input_name
              << ".root" << std::endl;
    return events;
  }
  TTree *tree = static_cast<TTree *>(file->Get("calibrated_tree"));
  if (tree)
    events = RooFitUtils::LoadEventsFromTree(tree, "depositedEnergykeV");
  file->Close();
  delete file;
  return events;
}

// Counts-per-second histogram of an already-calibrated event list. Sumw2 so the
// propagated errors survive the subtraction.
inline TH1F *RateHistFromEvents(const std::vector<Double_t> &events,
                                Double_t live_time, const TString &name,
                                Double_t lo, Double_t hi) {
  Int_t nbins = (Int_t)std::lround((hi - lo) / Constants::BIN_WIDTH_KEV);
  TH1F *hist =
      new TH1F(name, ";Deposited Energy [keV];Rate [counts/s]", nbins, lo, hi);
  hist->SetDirectory(nullptr);
  hist->Sumw2();
  for (size_t i = 0; i < events.size(); i++)
    hist->Fill(events[i]);
  if (live_time > 0)
    hist->Scale(1.0 / live_time);
  return hist;
}

// ------------------------------------------------------------ simultaneous

struct BkgGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid = kFALSE;
};

// bkg channel: the Pb Ka peaks (1 or 2). sig channel: the same peaks plus the
// Ge line at E_NEAR_GE. One low-side tail and one Pb shape are shared across
// the two channels; see the comments in the body for why.
inline BkgGeSimResult
RunBkgGeSim(const std::vector<Double_t> &bkg_events,
            const std::vector<Double_t> &sig_events,
            const std::vector<Double_t> &bkg_peak_mus, const TString &bkg_label,
            const TString &sig_label, const TString &fit_label,
            const FitResult *precomputed_bkg_seed = nullptr,
            const FitResult *precomputed_sig_seed = nullptr) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();
  if (n_bkg < 1 || n_bkg > 2 || bkg_events.empty() || sig_events.empty())
    return out;

  // No high-side exponential tail: it never takes an interior value across
  // the hand-fitted runs, either switching off or railing at both bounds at
  // once.
  const Bool_t kFlatBkg = USE_FLAT_BKG;
  const Bool_t kStep = kFALSE;
  const Bool_t kLowExp = kTRUE;
  const Bool_t kLowLin = kFALSE;
  const Bool_t kHighExp = kFALSE;

  FitResult seed;
  if (precomputed_bkg_seed) {
    seed = *precomputed_bkg_seed;
  } else {
    RooFitUtils bkg_fitter(bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                           Constants::BIN_WIDTH_KEV, kFlatBkg, kStep, kLowExp,
                           kLowLin, kHighExp);
    bkg_fitter.SetTailRatioMax(TAIL_RATIO_MAX);
    bkg_fitter.SetInteractive();
    TString seed_label = fit_label + "_BkgSeed";
    if (n_bkg == 1)
      seed = bkg_fitter.FitSinglePeak(bkg_label, seed_label);
    else
      seed = bkg_fitter.FitDoublePeak(bkg_label, seed_label, bkg_peak_mus[0],
                                      bkg_peak_mus[1]);
    if (!seed.valid) {
      std::cerr << "ERROR: bkg-only seed fit failed for " << bkg_label
                << std::endl;
      return out;
    }
  }

  // The saved .simroofits carries the hand-tuned starting point a cold start
  // cannot reach; SetRefitAfterLoad makes it SEED a real fit rather than be
  // returned as the result, which once produced post-cal rows holding day-one
  // constants at chi2 7-9.
  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  sim.SetInteractive();
  sim.SetRefitAfterLoad();
  if (HardFitRequested()) {
    sim.SetFitStrategy(2);
    sim.SetFitFinalHesse(kTRUE);
    sim.SetFitMaxCalls(200000);
  }

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back((Double_t)E_NEAR_GE);
  std::vector<Bool_t> bkg_fixed(n_bkg, kFALSE);
  std::vector<Bool_t> sig_fixed(n_bkg + 1, kFALSE);

  std::vector<Bool_t> sig_step(n_bkg + 1, kTRUE);
  sig_step[n_bkg] = kFALSE;

  sim.AddChannel("bkg", bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                 Constants::BIN_WIDTH_KEV, n_bkg, bkg_peak_mus, kFlatBkg, kStep,
                 kLowExp, kLowLin, kHighExp, bkg_fixed, kFALSE, kTRUE, kFALSE);
  sim.AddChannel("sig", sig_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                 Constants::BIN_WIDTH_KEV, n_bkg + 1, sig_mus, kFlatBkg, kStep,
                 kLowExp, kLowLin, kHighExp, sig_fixed, kFALSE, kTRUE, kFALSE,
                 sig_step);
  // Share ONE low-side tail (LowExp + LowLin) across both Pb-Ka lines and the
  // Ge peak, tied to bkg:Ka1. Fit independently, Ka1 inflates its own tail to
  // absorb the inter-peak fill so neither Ka tail is trustworthy, and the Ge
  // tail rails to zero, biasing its centroid low. Ge mu/sigma/yield stay free.
  // Pushed BEFORE LinkPeakShape so they win first-match on the tail params.
  // Route every target straight to bkg:Ka1 -- chaining via Ka2 hits a
  // ResolveOrCreate quirk that orphans sig:Ka2's tail into a free parameter.
  // Only link components that exist: with kLowLin off the LowLinTail params
  // are never built and linking to them aborts the channel build.
  const char *tp[4] = {"LowExpTailAmplitude", "LowExpTailRatio",
                       "LowLinTailAmplitude", "LowLinTailSlope"};
  const Int_t n_tp = kLowLin ? 4 : 2;
  TString ge_idx = TString::Format("%d", n_bkg + 1);
  for (Int_t t = 0; t < n_tp; t++) {
    if (n_bkg == 2)
      sim.LinkParameter(TString("bkg:") + tp[t] + "2",
                        TString("bkg:") + tp[t] + "1");
    for (Int_t i = 1; i < n_bkg; i++)
      sim.LinkParameter(TString("sig:") + tp[t] + TString::Format("%d", i + 1),
                        TString("bkg:") + tp[t] + "1");
    sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                      TString("bkg:") + tp[t] + "1");
  }

  if (LINK_PB_KA_SIGMA && n_bkg == 2)
    sim.LinkParameter("bkg:Sigma2", "bkg:Sigma1");

  // LinkPeakShape ties Mu along with the shape. With PB_ANCHOR=signal the
  // signal run's Pb lines need their own centroids, so the same links are made
  // one by one without Mu.
  const char *shape_without_mu[8] = {"Sigma",
                                     "StepAmplitude",
                                     "LowExpTailAmplitude",
                                     "LowExpTailRatio",
                                     "LowLinTailAmplitude",
                                     "LowLinTailSlope",
                                     "HighExpTailAmplitude",
                                     "HighExpTailRatio"};
  for (Int_t i = 0; i < n_bkg; i++) {
    if (!SignalPbAnchorRequested()) {
      sim.LinkPeakShape("sig", i, "bkg", i);
      continue;
    }
    TString idx = TString::Format("%d", i + 1);
    for (Int_t p = 0; p < 8; p++)
      sim.LinkParameter(TString("sig:") + shape_without_mu[p] + idx,
                        TString("bkg:") + shape_without_mu[p] + idx);
  }

  sim.SeedChannel("bkg", seed);
  if (precomputed_sig_seed)
    sim.SeedChannel("sig", *precomputed_sig_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_label, fit_label);
  if (results.size() < 2)
    return out;
  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;
  return out;
}

// ------------------------------------------------------------- rate-sub fit

struct RateSubResult {
  Float_t mu = -1, mu_err = -1, chi2 = -1;
  Float_t sigma = -1, bkg0 = 0, slope = 0;
  Float_t tail_amp = -1, tail_ratio = -1; // -1 when no tail was fitted
  Bool_t valid = kFALSE;
};

// Gaussian core plus the library's low-energy exponential tail on a linear
// background. Parameters: [0] gaussian height, [1] mu, [2] sigma, [3] tail
// height, [4] tau/sigma, [5] linear-tail height, [6] linear-tail slope,
// [7] bkg0, [8] bkg slope. Uses FittingFunctions, i.e. the same tail shape the
// production peak fits use, so this is not a reimplementation.
inline Double_t RateSubTailModel(Double_t *x, Double_t *p) {
  Double_t gaus[3] = {p[1], p[2], p[0]};
  Double_t tail[6] = {p[1], p[2], p[3], p[4], p[5], p[6]};
  return FittingFunctions::Gaussian(x, gaus) +
         FittingFunctions::LowTail(x, tail) + p[7] + p[8] * (x[0] - 68.75);
}

// Sentinel for "no background hint": compute it from the residual's left edge.
const Double_t kNoBkgHint = -1e30;

// Fit the ALREADY-SUBTRACTED residual histogram. Callers own the histogram and
// the model: pass out_model to keep the TF1 for drawing (the production macro
// does), otherwise it is deleted here.
//
// chi2 fit, not FittingUtils: its Poisson option assumes integer counts, but
// this histogram holds counts/second and goes negative after subtraction. The
// histogram must have Sumw2 set so propagated errors exist.
// Linear background: the 01/16 residual runs -2 to +2 counts/s across the
// window and a flat model drags the centroid high.
inline RateSubResult FitRateSubHist(TH1F *res, Double_t lo, Double_t hi,
                                    const TString &tag,
                                    TF1 **out_model = nullptr,
                                    Double_t mu_margin = 1.0,
                                    Double_t bkg0_hint = kNoBkgHint,
                                    Bool_t use_low_exp_tail = kFALSE) {
  RateSubResult out;
  if (!res)
    return out;

  Double_t bkg0 = 0.0;
  if (bkg0_hint > kNoBkgHint / 2) {
    bkg0 = bkg0_hint;
  } else {
    Double_t edge_lo =
        res->Integral(res->FindBin(lo + 0.2), res->FindBin(66.0));
    Double_t nedge = res->FindBin(66.0) - res->FindBin(lo + 0.2) + 1;
    bkg0 = (nedge > 0) ? edge_lo / nedge : 0.0;
  }
  Int_t pk_bin = res->GetMaximumBin();
  Double_t pk_val = res->GetBinContent(pk_bin) - bkg0;
  if (pk_val < 0)
    pk_val = 0;
  Double_t sigma0 = 0.85;
  static Int_t model_counter = 0;
  TString mname = Form("ratesub_model_%s_%d", tag.Data(), ++model_counter);
  TF1 *model = nullptr;
  if (use_low_exp_tail) {
    // Gaussian + library low-energy tail on a linear background.
    //   [0] gaus height  [1] mu  [2] sigma  [3] tail height  [4] tau/sigma
    //   [5] lin-tail height  [6] lin-tail slope  [7] bkg0  [8] bkg slope
    model = new TF1(mname, RateSubTailModel, lo, hi, 9);
    model->SetParNames("GausHeight", "Mu", "Sigma", "TailHeight",
                       "TauOverSigma", "LinTailHeight", "LinTailSlope", "Bkg0",
                       "BkgSlope");
    model->SetParameters(pk_val, E_NEAR_GE, sigma0, 0.05 * pk_val, 5.0, 0.0,
                         0.0, bkg0, 0.0);
    model->FixParameter(5, 0.0);
    model->FixParameter(6, 0.0);
    model->SetParLimits(3, 0.0, 4.0 * (pk_val > 0 ? pk_val : 1.0));
    model->SetParLimits(4, 1.0, TAIL_RATIO_MAX);
  } else {
    // Gaussian core on a linear background. The subtraction leaves a peak on a
    // near-flat residual, so no tail is fitted here.
    //   [0] area   [1] mu   [2] sigma   [3] bkg0   [4] bkg slope
    const char *kGausLinBkg =
        "[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*(x-68.75)";
    model = new TF1(mname, kGausLinBkg, lo, hi);
    model->SetParNames("Area", "Mu", "Sigma", "Bkg0", "BkgSlope");
    model->SetParameters(pk_val * sigma0 * 2.5066, E_NEAR_GE, sigma0, bkg0,
                         0.0);
  }
  model->SetParLimits(1, lo + mu_margin, hi - mu_margin);
  model->SetParLimits(2, 0.2, 3.0);
  TFitResultPtr fr = res->Fit(model, "S Q R N");
  Bool_t ok = fr.Get() && fr->IsValid() && model->GetNDF() > 0;
  if (ok) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.sigma = model->GetParameter(2);
    if (use_low_exp_tail) {
      out.tail_amp = model->GetParameter(3);
      out.tail_ratio = model->GetParameter(4);
      out.bkg0 = model->GetParameter(7);
      out.slope = model->GetParameter(8);
    } else {
      out.bkg0 = model->GetParameter(3);
      out.slope = model->GetParameter(4);
    }
    out.valid = kTRUE;
  }
  if (out_model)
    *out_model = model;
  else
    delete model;
  return out;
}

// ------------------------------------------- rate-sub fit, shared lineshape

// The residual alone cannot fix the Ge tail: it is degenerate with the linear
// background. So the rate-sub fit borrows the tail from the same day's
// reference lines, which puts its centroid on the same definition as the
// calibration anchors and the in-situ fit.

// Tail/Gaussian HEIGHT ratio, which the TF1 model below needs, from the tail's
// full-range yield ratio R = k * 2 rho exp(1/(2 rho^2)) / sqrt(2 pi). R is what
// RooFitUtils reports (low_exp_tail_amplitude / gaus_amplitude) and does not
// depend on the fit window.
inline Double_t TailHeightRatio(const PeakFitResult &pk) {
  if (!(pk.gaus_amplitude > 0) || !(pk.low_exp_tail_amplitude > 0) ||
      !(pk.low_exp_tail_ratio > 0))
    return 0.0;
  Double_t rho = pk.low_exp_tail_ratio;
  Double_t full_range_ratio = pk.low_exp_tail_amplitude / pk.gaus_amplitude;
  return full_range_ratio * TMath::Sqrt(2.0 * TMath::Pi()) /
         (2.0 * rho * TMath::Exp(1.0 / (2.0 * rho * rho)));
}

// f_T / f_G: the fraction of a low exponential tail inside [lo, hi] over the
// fraction of its Gaussian, evaluated by RooFitUtils' own factor. A full-range
// ratio R times this is the tail's in-window yield ratio for that peak.
inline Double_t LowExpTailWindowFactor(Double_t lo, Double_t hi, Double_t mu,
                                       Double_t sigma, Double_t rho) {
  RooRealVar m("window_factor_mu", "", mu);
  RooRealVar s("window_factor_sigma", "", sigma);
  RooRealVar r("window_factor_rho", "", rho);
  RooRealVar a("window_factor_lo", "", lo);
  RooRealVar b("window_factor_hi", "", hi);
  RooExpTailWindowFactor factor("window_factor", "", m, s, r, a, b, kFALSE);
  return factor.getVal();
}

// What the rate-sub fit takes as given: the Ge tail.
struct RateSubShape {
  Double_t ge_tail_height_ratio = 0;
  Double_t ge_tau_over_sigma = 1;
};

// Ge Gaussian + tail on a linear background. Parameters:
//   [0] Ge height  [1] Ge mu  [2] Ge sigma  [3] Ge tail/gaus height ratio
//   [4] Ge tau/sigma  [5] bkg0  [6] bkg slope
inline Double_t RateSubShapedModel(Double_t *x, Double_t *p) {
  Double_t ge_gaus[3] = {p[1], p[2], p[0]};
  Double_t ge_tail[6] = {p[1], p[2], p[0] * p[3], p[4], 0.0, 0.0};
  return FittingFunctions::Gaussian(x, ge_gaus) +
         FittingFunctions::LowTail(x, ge_tail) + p[5] + p[6] * (x[0] - 68.75);
}

// Fit the already-subtracted residual over the production window with the
// borrowed shape, seeded from the plain Gaussian fit.
inline RateSubResult FitRateSubHistShaped(TH1F *res, const RateSubShape &shape,
                                          const TString &tag,
                                          TF1 **out_model = nullptr) {
  RateSubResult out;
  if (!res)
    return out;
  TF1 *seed = nullptr;
  RateSubResult gaus =
      FitRateSubHist(res, RATESUB_LO, RATESUB_HI, tag + "_seed", &seed);
  if (!gaus.valid || !seed) {
    delete seed;
    return out;
  }
  const char *names[7] = {"GeHeight",   "GeMu", "GeSigma", "GeTailRatio",
                          "GeTauRatio", "Bkg0", "BkgSlope"};
  static Int_t model_counter = 0;
  TString mname = Form("ratesub_shaped_%s_%d", tag.Data(), ++model_counter);
  TF1 *model = new TF1(mname, RateSubShapedModel, RATESUB_LO, RATESUB_HI, 7);
  for (Int_t i = 0; i < 7; i++)
    model->SetParName(i, names[i]);
  Double_t sigma0 = seed->GetParameter(2);
  model->SetParameter(0, seed->GetParameter(0) /
                             (sigma0 * TMath::Sqrt(2.0 * TMath::Pi())));
  model->SetParameter(1, seed->GetParameter(1));
  model->SetParameter(2, sigma0);
  model->FixParameter(3, shape.ge_tail_height_ratio);
  model->FixParameter(4, shape.ge_tau_over_sigma);
  model->SetParameter(5, seed->GetParameter(3));
  model->SetParameter(6, seed->GetParameter(4));
  model->SetParLimits(1, RATESUB_LO + 1.0, RATESUB_HI - 1.0);
  model->SetParLimits(2, 0.2, 3.0);
  delete seed;

  TFitResultPtr fr = res->Fit(model, "S Q R N");
  if (fr.Get() && fr->IsValid() && model->GetNDF() > 0) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.sigma = model->GetParameter(2);
    out.tail_amp = shape.ge_tail_height_ratio;
    out.tail_ratio = shape.ge_tau_over_sigma;
    out.bkg0 = model->GetParameter(5);
    out.slope = model->GetParameter(6);
    out.valid = kTRUE;
  }
  if (out_model)
    *out_model = model;
  else
    delete model;
  return out;
}

#endif // GAINTRANSFER_FIT_HPP
