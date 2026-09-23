#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <RtypesCore.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

// calibration
const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_79 = 79.6142;
const Float_t E_BA133_81 = 80.9979;
const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;
const Float_t E_PB_KB1 = 84.936;

// just used for seeding fits
const Float_t E_PB_KB2 = 87.32;
const Float_t E_GE_73M = 68.752;

// Anchor reference-energy uncertainties (keV), folded into the cal-fit x-errors
// so they propagate into the pol1 parameter covariance and thus the Ge energy.
// Pb Ka are known to <1 eV (Deslattes); Kb1 carries the ~2:1-vs-1.9:1 intensity
// -blend ambiguity (~6 eV).
//
// dE_LINE remains the conservative ~10 eV lump, now used ONLY for the Ba-133
// lines of the 15th/16th line-cal days.
// TODO: replace with per-line evaluated Ba-133 uncertainties.
const Float_t dE_PB_KA = 0.001;
const Float_t dE_PB_KB1 = 0.006;
const Float_t dE_AM241 = 0.0001;
const Float_t dE_LINE = 0.010;

const Float_t CAL_RANGE_LOW = 61.5;
const Float_t CAL_RANGE_HIGH = 80.5;

const Bool_t TAIL_SHARED = kTRUE;

const Int_t REF_CAL_DEGREE = 2;
const Bool_t LINK_PB_KA_SIGMA = kTRUE;

// Constrain the pre-cal Pb-Ka spacing to tabulated * local gain (the doublet
// sum fixes the gain, the difference is constrained, so it does not feed on
// itself). OFF: the 14-23 eV excess is 1.2-1.9 sigma against this
// measurement's own ~12 eV separation uncertainty -- ordinary scatter.
const Bool_t CONSTRAIN_PRECAL_SEPARATION = kFALSE;
const Bool_t CONSTRAIN_PB_KA_SEPARATION = kFALSE;

const Double_t dE_PB_KA_SEP = 0.0015;
const Bool_t LOCK_POSTCAL_PB = kFALSE;

// High-side exponential tail. OFF: it never takes an interior value across
// the hand-fitted runs, either switching off or railing at both bounds at
// once. A 100 keV decay across a 19 keV window is a flat pedestal, degenerate
// with BkgConstant.
const Bool_t USE_HIGH_EXP_TAIL = kFALSE;

const Bool_t USE_FLAT_BKG = kTRUE;

// Step (charge-collection shelf) on the calibration-peak fits. OFF: tested
// on Am-241, the step floats to 73 while the tail decay is unchanged to four
// digits, so the tails are not standing in for a missing shelf.
const Bool_t USE_STEP_CAL_PEAKS = kFALSE;
const Double_t TAIL_RATIO_MAX = 8.0;

const Bool_t USE_AFFINE_TRANSFER = kFALSE;

// Rate-subtraction tail, held fixed. Nominal is no tail; set FRAC to the
// in-situ lineshape's tail fraction to measure the systematic.
const Double_t RATESUB_TAIL_FRAC = 0.0;
const Double_t RATESUB_TAIL_TAU = 1.0;

// Rate-subtract the 13th's pairs too, for one method across all six
// datasets. OFF: those rows land 110-150 eV below their own in-situ values
// for reasons not understood, blowing the method systematic to 80 eV.
const Bool_t RATESUB_ALL_DAYS = kFALSE;
const Bool_t RATESUB_TAIL_SCAN = kFALSE;

const TString RESULT_TAG = "amxfer";

const TString RESULT_VARIANT = USE_HIGH_EXP_TAIL ? "_hitail" : "_nohitail";
const TString OUT_SUBDIR = "calibrated_low_" + RESULT_TAG;

TString G_CSV_PREFIX = "";

struct CalibrationData {
  std::vector<Float_t> mu, mu_errors, calibration_values_keV, reduced_chi2;
  std::vector<Float_t> energy_errors;
  std::vector<TString> run_names;
};

struct BkgGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid = kFALSE;
};

struct DayResult {
  std::vector<TString> cal_labels;
  std::vector<TF1 *> cal_funcs;
  std::vector<TString> ge_labels;
  std::vector<Float_t> ge_mus, ge_errs, ge_chi2s;
  std::vector<Float_t> ge_cal_errs;
  // The OFFSET part of the reference calibration error: its value at the
  // best-constrained energy. The remainder, sqrt(ref^2 - off^2), is the
  // multiplicative term, and the gain transfer is a pure scale factor.
  std::vector<Float_t> ge_cal_offs;
  // The TRANSFER part of ge_cal_errs, parallel to ge_*. ge_cal_err is
  // hypot(reference_cal_err, gain_err), so the reference part -- COMMON to
  // every run on the day -- is sqrt(cal^2 - gain^2) and the transfer part is
  // independent per run. The combiner needs that split to build the covariance
  // correctly instead of assuming the whole calibration is correlated or none
  // of it is.
  std::vector<Float_t> ge_gain_errs;
  std::vector<Float_t> ge_bkg_errs;
  std::vector<TString> rs_labels;
  std::vector<Float_t> rs_mus, rs_errs, rs_chi2s, rs_cal_errs, rs_gain_errs;
  std::vector<Float_t> rs_cal_offs;
};

struct CalFit {
  TF1 *func = nullptr;
  Double_t p0 = 0, p1 = 1;
  Int_t npar = 2;
  Double_t var_p0 = 0, var_p1 = 0, cov_p0p1 = 0;
  Double_t var_p2 = 0, cov_p0p2 = 0, cov_p1p2 = 0;
};

struct PairCalResult {
  TF1 *cal_func = nullptr;
  Double_t ge_cal_err = 0; // calibration sigma propagated to the Ge energy
  Double_t ge_cal_off = 0; // its value at the best-constrained energy

  Double_t ge_bkg_err = 0; // |cal(precal Ge) - postcal Ge|; see file-top note
  BkgGeSimResult postcal;
  TString pair_tag;
  std::vector<Float_t> pb_mus_precal;
  std::vector<Float_t> pb_mu_errs_precal;
  CalFit ref_calfit;
  Bool_t has_ref_cal = kFALSE;
  Double_t ge_gain_err = 0;
};

struct LineCalConfig {
  TString date_label;
  TString am_run;
  Float_t am_lo, am_hi;
  TString ba_run;
  std::vector<TString> day_datasets;
  TString postcal_bkg;
  TString postcal_sig;
  TString ge_label;
  Float_t rs_lo = 62.0;
  Float_t rs_hi = 71.5;
  Bool_t rs_linear_bkg = kTRUE;
};

std::vector<Double_t> LoadFiltered(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    return events;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    std::vector<Double_t> ch =
        RooFitUtils::LoadEventsFromTree(tree, "energykeV");
    events.insert(events.end(), ch.begin(), ch.end());
  }
  file->Close();
  delete file;
  return events;
}

std::vector<Double_t> LoadCalibrated(const TString &input_name) {
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

FitResult FitCalPeak(const std::vector<Double_t> &events,
                     const TString &input_name, const TString &peak_name,
                     Float_t fit_low, Float_t fit_high, Bool_t double_peak,
                     Double_t mu1, Double_t mu2, Bool_t interactive,
                     Bool_t link_sigma = kFALSE, Bool_t use_step = kFALSE) {
  if (events.empty())
    return {};
  RooFitUtils fitter(events, fit_low, fit_high, Constants::BIN_WIDTH_KEV, kTRUE,
                     use_step, kTRUE, kTRUE, kTRUE);
  fitter.SetTailRatioMax(TAIL_RATIO_MAX);
  if (interactive)
    fitter.SetInteractive();
  if (double_peak)
    return fitter.FitDoublePeak(input_name, peak_name, mu1, mu2, link_sigma);
  return fitter.FitSinglePeak(input_name, peak_name);
}

BkgGeSimResult
RunBkgGeSim(const std::vector<Double_t> &bkg_events,
            const std::vector<Double_t> &sig_events,
            const std::vector<Double_t> &bkg_peak_mus, const TString &bkg_label,
            const TString &sig_label, const TString &fit_label,
            Bool_t interactive, Bool_t lock_bkg_peaks = kFALSE,
            const FitResult *precomputed_bkg_seed = nullptr,
            const FitResult *precomputed_sig_seed = nullptr,
            Bool_t share_tail = kTRUE, Float_t fit_lo = CAL_RANGE_LOW,
            Float_t fit_hi = CAL_RANGE_HIGH, Bool_t force_sim_fit = kFALSE,
            Double_t constrain_sep_delta = -1.0,
            Double_t constrain_sep_sigma = -1.0) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();
  if (n_bkg < 1 || n_bkg > 2 || bkg_events.empty() || sig_events.empty())
    return out;

  const Bool_t kFlatBkg = USE_FLAT_BKG;
  const Bool_t kStep = kFALSE;
  const Bool_t kLowExp = kTRUE;
  const Bool_t kLowLin = kTRUE;
  const Bool_t kHighExp = USE_HIGH_EXP_TAIL;

  FitResult seed;
  if (precomputed_bkg_seed) {
    seed = *precomputed_bkg_seed;
  } else {
    RooFitUtils bkg_fitter(bkg_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                           kFlatBkg, kStep, kLowExp, kLowLin, kHighExp);
    bkg_fitter.SetTailRatioMax(TAIL_RATIO_MAX);
    if (interactive)
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

  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  if (interactive) {
    sim.SetInteractive();
    if (force_sim_fit)
      sim.SetRefitAfterLoad();
  }

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back((Double_t)E_GE_73M);
  std::vector<Bool_t> bkg_fixed(n_bkg, lock_bkg_peaks);
  std::vector<Bool_t> sig_fixed(n_bkg + 1, kFALSE);
  for (Int_t i = 0; i < n_bkg; i++)
    sig_fixed[i] = lock_bkg_peaks;

  std::vector<Bool_t> sig_step(n_bkg + 1, kTRUE);
  sig_step[n_bkg] = kFALSE;

  sim.AddChannel("bkg", bkg_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg, bkg_peak_mus, kFlatBkg, kStep, kLowExp, kLowLin,
                 kHighExp, bkg_fixed, kFALSE, kTRUE, kFALSE);
  sim.AddChannel("sig", sig_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg + 1, sig_mus, kFlatBkg, kStep, kLowExp, kLowLin,
                 kHighExp, sig_fixed, kFALSE, kTRUE, kFALSE, sig_step);
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
  if (share_tail) {
    for (Int_t t = 0; t < n_tp; t++) {
      if (n_bkg == 2)
        sim.LinkParameter(TString("bkg:") + tp[t] + "2",
                          TString("bkg:") + tp[t] + "1");
      for (Int_t i = 1; i < n_bkg; i++)
        sim.LinkParameter(TString("sig:") + tp[t] +
                              TString::Format("%d", i + 1),
                          TString("bkg:") + tp[t] + "1");
      sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                        TString("bkg:") + tp[t] + "1");
    }
  } else {
    for (Int_t t = 0; t < n_tp; t++)
      sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                        TString("bkg:") + tp[t] + "1");
  }

  if (constrain_sep_delta > 0 && n_bkg == 2)
    sim.ConstrainPeakSeparation("bkg", 1, 0, constrain_sep_delta,
                                constrain_sep_sigma > 0 ? constrain_sep_sigma
                                                        : dE_PB_KA_SEP);

  if (LINK_PB_KA_SIGMA && n_bkg == 2)
    sim.LinkParameter("bkg:Sigma2", "bkg:Sigma1");

  for (Int_t i = 0; i < n_bkg; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  sim.SeedChannel("bkg", seed);
  if (precomputed_sig_seed)
    sim.SeedChannel("sig", *precomputed_sig_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_label, fit_label);
  if (results.size() < 2)
    return out;
  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;

  if (!G_CSV_PREFIX.IsNull()) {
    sim.DumpChannelCSV("bkg", G_CSV_PREFIX + "_bkg");
    sim.DumpChannelCSV("sig", G_CSV_PREFIX + "_sig");
  }
  return out;
}

void AddCalPoint(CalibrationData &cal_data, const TString &name, Float_t mu,
                 Float_t mu_err, Float_t true_e, Float_t chi2,
                 Float_t energy_err = 0) {
  cal_data.run_names.push_back(name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_err);
  cal_data.calibration_values_keV.push_back(true_e);
  cal_data.reduced_chi2.push_back(chi2);
  cal_data.energy_errors.push_back(energy_err);
}

// Cramer-Rao floor on a fitted centroid, sigma/sqrt(N). A fit reporting less
// has a broken covariance, not a better measurement: Pb-Kb1 and Am-241 were
// reaching the weighted cal fit with ~0 error and outvoting every other
// anchor. A floor only -- a peak that trips it still wants re-tuning.
Float_t GuardedMuError(const TString &name, const PeakFitResult &pk,
                       const FitResult &fr) {
  if (!(pk.sigma > 0) || !(pk.gaus_amplitude > 0))
    return pk.mu_error;
  const Float_t floor_err = pk.sigma / TMath::Sqrt(pk.gaus_amplitude);
  if (pk.mu_error > floor_err)
    return pk.mu_error;
  std::cerr << "WARNING: " << name << " centroid error " << pk.mu_error * 1000
            << " eV is below the Poisson floor " << floor_err * 1000
            << " eV (sigma " << pk.sigma << ", N " << pk.gaus_amplitude
            << ", covQual " << fr.cov_qual << ") -- raising to the floor. "
            << "Re-tune this fit's saved state under TAIL_RATIO_MAX."
            << std::endl;
  return floor_err;
}

void PrintCalSummary(const CalibrationData &cal_data,
                     const TString &date_label) {
  std::cout << "Calibration points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); i++) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0)
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
    std::cout << std::endl;
  }
}

Double_t CalEnergyError(const CalFit &c, Double_t mu) {
  Double_t var = c.var_p0 + mu * mu * c.var_p1 + 2.0 * mu * c.cov_p0p1;
  if (c.npar >= 3) {
    Double_t mu2 = mu * mu;
    var += mu2 * mu2 * c.var_p2 + 2.0 * mu2 * c.cov_p0p2 +
           2.0 * mu2 * mu * c.cov_p1p2;
  }
  return (var > 0) ? std::sqrt(var) : 0.0;
}

// The calibration uncertainty at its best-constrained energy, the minimum of
// CalEnergyError over the anchor region. For a straight line this is the
// intercept error at the pivot; sqrt(cal^2 - off^2) is then the
// multiplicative term. The same definition is applied to the pol2.
Double_t CalOffsetError(const CalFit &c) {
  Double_t best = 0;
  for (Double_t x = 50.0; x <= 100.0; x += 0.01) {
    Double_t e = CalEnergyError(c, x);
    if (e > 0 && (best <= 0 || e < best))
      best = e;
  }
  return best;
}

// Parameter covariance as (J^T W J)^-1 with the effective variances
// ey^2 + (f'(x) ex)^2 at the fitted parameters, exact for the linearised
// problem. On raw keV the pol2 parameters are 99% anticorrelated and Minuit's
// numerical Hesse moved the propagated Ge error between 3.2 and 4.6 eV with
// minimizer and strategy; this gives 5.5 eV independent of both.
void FillCalCovariance(CalFit &c, const std::vector<Float_t> &x,
                       const std::vector<Float_t> &ex,
                       const std::vector<Float_t> &ey, Bool_t fix_p0) {
  std::vector<Int_t> free_par;
  for (Int_t k = (fix_p0 ? 1 : 0); k < c.npar; k++)
    free_par.push_back(k);
  Int_t nf = (Int_t)free_par.size();
  TMatrixD H(nf, nf);
  for (size_t i = 0; i < x.size(); i++) {
    Double_t xi = x[i];
    Double_t fp = c.func->Derivative(xi);
    Double_t s2 = (Double_t)ey[i] * ey[i] + fp * fp * (Double_t)ex[i] * ex[i];
    for (Int_t a = 0; a < nf; a++)
      for (Int_t b = 0; b < nf; b++)
        H(a, b) += std::pow(xi, free_par[a]) * std::pow(xi, free_par[b]) / s2;
  }
  H.Invert();
  Double_t cov[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (Int_t a = 0; a < nf; a++)
    for (Int_t b = 0; b < nf; b++)
      cov[free_par[a]][free_par[b]] = H(a, b);
  c.var_p0 = cov[0][0];
  c.var_p1 = cov[1][1];
  c.cov_p0p1 = cov[0][1];
  if (c.npar >= 3) {
    c.var_p2 = cov[2][2];
    c.cov_p0p2 = cov[0][2];
    c.cov_p1p2 = cov[1][2];
  }
  for (Int_t k = 0; k < c.npar; k++)
    c.func->SetParError(k, std::sqrt(std::max(0.0, cov[k][k])));
}

CalFit CreateAndSavePol1Cal(const CalibrationData &cal_data,
                            const TString &date_label, Bool_t fix_p0 = kFALSE,
                            Float_t p0_value = 0, Int_t degree = 1) {
  Int_t n = (Int_t)cal_data.mu.size();
  std::vector<Float_t> ex(cal_data.mu_errors);
  std::vector<Float_t> ey(n);
  for (Int_t i = 0; i < n; i++) {
    Float_t e = (i < (Int_t)cal_data.energy_errors.size())
                    ? cal_data.energy_errors[i]
                    : 0.0f;
    ey[i] = (e > 0)
                ? e
                : 1e-4f; // tiny floor so the point isn't infinitely weighted
    if (!(ex[i] > 0))
      ex[i] = 1e-4f;
  }
  TGraphErrors *graph = new TGraphErrors(n, cal_data.mu.data(),
                                         cal_data.calibration_values_keV.data(),
                                         ex.data(), ey.data());
  Double_t draw_lo = CAL_RANGE_LOW, draw_hi = CAL_RANGE_HIGH;
  for (Int_t i = 0; i < n; i++) {
    draw_lo = std::min({draw_lo, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
    draw_hi = std::max({draw_hi, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
  }
  draw_lo -= 3.0;
  draw_hi += 3.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      graph, kBlue, "; Precalibrated Energy [keV]; Deposited Energy [keV]");
  graph->GetXaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetYaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetMarkerStyle(5);
  graph->SetMarkerSize(2);
  graph->Draw("AP");

  TF1 *cal = new TF1("cal_low_" + date_label, (degree >= 2) ? "pol2" : "pol1",
                     draw_lo, draw_hi);
  cal->SetParameter(0, 0);
  cal->SetParameter(1, 1);
  if (degree >= 2)
    cal->SetParameter(2, 0);
  if (fix_p0)
    cal->FixParameter(0, p0_value);
  cal->SetNpx(1000);
  TFitResultPtr fr = graph->Fit(cal, "E S Q");
  cal->Draw("SAME");

  CalFit out;
  out.func = cal;
  out.npar = cal->GetNpar();
  out.p0 = cal->GetParameter(0);
  out.p1 = cal->GetParameter(1);
  FillCalCovariance(out, cal_data.mu, ex, ey, fix_p0);
  if (fr.Get() && fr->Ndf() > 0)
    std::cout << "CAL-FIT " << date_label << " (pol" << (out.npar - 1) << ", "
              << n << " anchors): chi2/ndf = " << std::fixed
              << std::setprecision(2) << fr->Chi2() / fr->Ndf()
              << "   sigma_E at best-constrained point = "
              << std::setprecision(1) << CalOffsetError(out) * 1000.0 << " eV"
              << std::endl;

  PlottingUtils::SaveFigure(canvas, "calibration_low_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return out;
}

struct AffineResult {
  Double_t a = 0.0; // offset [keV]
  Double_t b = 1.0; // gain
  Double_t var_a = 0, var_b = 0, cov_ab = 0;
  Double_t chi2_ndf = 0;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

AffineResult PbAffineVsReference(const std::vector<Float_t> &mu_run,
                                 const std::vector<Float_t> &err_run,
                                 const std::vector<Float_t> &mu_ref,
                                 const std::vector<Float_t> &err_ref) {
  AffineResult out;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  std::vector<Double_t> x, y, ex, ey;
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    x.push_back(mu_ref[i]);
    y.push_back(mu_run[i]);
    ex.push_back(i < err_ref.size() && err_ref[i] > 0 ? err_ref[i] : 1e-4);
    ey.push_back(i < err_run.size() && err_run[i] > 0 ? err_run[i] : 1e-4);
  }
  out.n = (Int_t)x.size();
  if (out.n < 3)
    return out; // need a spare point for a meaningful residual

  TGraphErrors g((Int_t)x.size(), x.data(), y.data(), ex.data(), ey.data());
  TF1 lin("pb_affine", "pol1", x.front() - 5.0, x.back() + 5.0);
  lin.SetParameters(0.0, 1.0);
  TFitResultPtr fr = g.Fit(&lin, "E S Q N");
  if (!fr.Get())
    return out;

  out.a = lin.GetParameter(0);
  out.b = lin.GetParameter(1);
  out.var_a = fr->CovMatrix(0, 0);
  out.var_b = fr->CovMatrix(1, 1);
  out.cov_ab = fr->CovMatrix(0, 1);
  out.chi2_ndf = (fr->Ndf() > 0) ? fr->Chi2() / fr->Ndf() : 0.0;
  if (out.chi2_ndf > 1.0) {
    out.var_a *= out.chi2_ndf;
    out.var_b *= out.chi2_ndf;
    out.cov_ab *= out.chi2_ndf;
  }
  out.valid = kTRUE;
  return out;
}

CalFit AffineScaleCal(const CalFit &ref, const AffineResult &t,
                      const TString &label) {
  TF1 *src = ref.func;
  Int_t npar = src->GetNpar();
  Double_t lo = 0, hi = 0;
  src->GetRange(lo, hi);
  TF1 *scaled = new TF1("cal_low_" + label, src->GetExpFormula(), lo, hi);

  Double_t p0 = src->GetParameter(0);
  Double_t p1 = (npar > 1) ? src->GetParameter(1) : 0.0;
  Double_t p2 = (npar > 2) ? src->GetParameter(2) : 0.0;
  Double_t a = t.a, b = (t.b != 0 ? t.b : 1.0);

  scaled->SetParameter(0, p0 - p1 * a / b + p2 * a * a / (b * b));
  if (npar > 1)
    scaled->SetParameter(1, p1 / b - 2.0 * p2 * a / (b * b));
  if (npar > 2)
    scaled->SetParameter(2, p2 / (b * b));
  scaled->SetNpx(1000);

  CalFit out;
  out.func = scaled;
  out.npar = npar;
  out.p0 = scaled->GetParameter(0);
  out.p1 = scaled->GetParameter(1);
  out.var_p0 = ref.var_p0;
  out.var_p1 = ref.var_p1 / (b * b);
  out.cov_p0p1 = ref.cov_p0p1 / b;
  if (npar >= 3) {
    Double_t b2 = b * b;
    out.var_p2 = ref.var_p2 / (b2 * b2);
    out.cov_p0p2 = ref.cov_p0p2 / b2;
    out.cov_p1p2 = ref.cov_p1p2 / (b2 * b);
  }
  return out;
}

Double_t AffineTransferEnergyError(const CalFit &ref_cal, Double_t mu_raw,
                                   const AffineResult &t) {
  if (!ref_cal.func || !t.valid || t.b == 0)
    return 0.0;
  Double_t z = (mu_raw - t.a) / t.b;
  Double_t fp = ref_cal.func->Derivative(z);
  Double_t dda = -fp / t.b;
  Double_t ddb = -fp * z / t.b;
  Double_t var =
      dda * dda * t.var_a + ddb * ddb * t.var_b + 2.0 * dda * ddb * t.cov_ab;
  return (var > 0) ? std::sqrt(var) : 0.0;
}

struct GainResult {
  Double_t g = 1.0; // error-weighted mean relative gain
  Double_t sigma_g = 0.0;
  Double_t chi2_ndf = 0.0;
  Int_t n = 0;
};

GainResult PbGainVsReferenceWithError(const std::vector<Float_t> &mu_run,
                                      const std::vector<Float_t> &err_run,
                                      const std::vector<Float_t> &mu_ref,
                                      const std::vector<Float_t> &err_ref) {
  GainResult out;
  Double_t sw = 0, swg = 0;
  std::vector<Double_t> gs, ws;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    Double_t g = (Double_t)mu_run[i] / (Double_t)mu_ref[i];
    Double_t fr = (err_run[i] > 0) ? err_run[i] / mu_run[i] : 0.0;
    Double_t ff = (err_ref[i] > 0) ? err_ref[i] / mu_ref[i] : 0.0;
    Double_t rel = std::sqrt(fr * fr + ff * ff);
    Double_t sg = (rel > 0) ? g * rel : 1e-9;
    Double_t w = 1.0 / (sg * sg);
    sw += w;
    swg += w * g;
    gs.push_back(g);
    ws.push_back(w);
  }
  if (sw <= 0)
    return out;
  out.g = swg / sw;
  out.sigma_g = std::sqrt(1.0 / sw);
  out.n = (Int_t)gs.size();
  if (out.n > 1) {
    Double_t chi2 = 0;
    for (size_t i = 0; i < gs.size(); i++)
      chi2 += ws[i] * (gs[i] - out.g) * (gs[i] - out.g);
    out.chi2_ndf = chi2 / (out.n - 1);
    if (out.chi2_ndf > 1.0)
      out.sigma_g *= std::sqrt(out.chi2_ndf);
  }
  return out;
}

CalFit GainScaleCal(const CalFit &ref, Double_t g, const TString &label) {
  TF1 *src = ref.func;
  Int_t npar = src->GetNpar();
  Double_t lo = 0, hi = 0;
  src->GetRange(lo, hi);
  TF1 *scaled = new TF1("cal_low_" + label, src->GetExpFormula(), lo, hi);
  Double_t gk = 1.0;
  for (Int_t k = 0; k < npar; k++) {
    scaled->SetParameter(k, src->GetParameter(k) / gk);
    gk *= g;
  }
  scaled->SetNpx(1000);
  CalFit out;
  out.func = scaled;
  out.p0 = scaled->GetParameter(0);
  out.p1 = scaled->GetParameter(1);
  out.npar = npar;
  out.var_p0 = ref.var_p0;
  out.var_p1 = ref.var_p1 / (g * g);
  out.cov_p0p1 = ref.cov_p0p1 / g;
  if (npar >= 3) {
    Double_t g2 = g * g;
    out.var_p2 = ref.var_p2 / (g2 * g2);
    out.cov_p0p2 = ref.cov_p0p2 / g2;
    out.cov_p1p2 = ref.cov_p1p2 / (g2 * g);
  }
  if (npar > 0)
    scaled->SetParError(0, std::sqrt(std::max(0.0, out.var_p0)));
  if (npar > 1)
    scaled->SetParError(1, std::sqrt(std::max(0.0, out.var_p1)));
  return out;
}

Double_t GainTransferEnergyError(const CalFit &ref_cal, Double_t mu_raw,
                                 Double_t g, Double_t sigma_g) {
  if (!(sigma_g > 0) || !(g > 0) || !ref_cal.func)
    return 0.0;
  Double_t fprime = ref_cal.func->Derivative(mu_raw / g);
  return std::fabs(-mu_raw / (g * g) * fprime) * sigma_g;
}

Float_t LoadRunStartTime(const TString &input_name) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name
              << ".root for run start time" << std::endl;
    if (file) {
      file->Close();
      delete file;
    }
    return 0;
  }
  Float_t t = 0;
  TTree *hdr = static_cast<TTree *>(file->Get("bef_header_tree"));
  if (hdr && hdr->GetEntries() > 0) {
    UInt_t timeStart = 0;
    hdr->SetBranchAddress("timeStart", &timeStart);
    hdr->GetEntry(0);
    t = (Float_t)timeStart;
  }
  file->Close();
  delete file;
  return t;
}

void ApplyPol1Cal(const std::vector<TString> &input_names, TF1 *cal,
                  const TString &date_label) {
  TFile *cal_file = IO::OpenForWriting(
      OUT_SUBDIR + "/calibration_function_low_" + date_label + ".root");
  cal->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  for (size_t i = 0; i < input_names.size(); i++) {
    TString input_name = input_names[i];
    TFile *in_file = IO::OpenForReading("filtered/" + input_name + ".root");
    if (!in_file || in_file->IsZombie()) {
      std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
                << std::endl;
      delete in_file;
      continue;
    }
    TFile *out_file =
        IO::OpenForWriting(OUT_SUBDIR + "/" + input_name + ".root");

    Float_t deposited_energy = 0;
    TTree *cal_tree = new TTree("calibrated_tree", "Calibrated events");
    cal_tree->SetDirectory(out_file);
    cal_tree->Branch("depositedEnergykeV", &deposited_energy,
                     "depositedEnergykeV/F");

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TTree *tree = static_cast<TTree *>(
          in_file->Get(Form("crystal%d_filtered_tree", c)));
      if (!tree)
        continue;
      Float_t energy = 0;
      tree->SetBranchAddress("energykeV", &energy);
      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        deposited_energy = cal->Eval(energy);
        cal_tree->Fill();
      }
    }

    out_file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    std::cout << "Calibrated " << input_name << " -> " << OUT_SUBDIR << "/"
              << input_name << ".root" << std::endl;
    out_file->Close();
    delete out_file;
    in_file->Close();
    delete in_file;
  }
}

void AddGeResult(DayResult &result, const BkgGeSimResult &r,
                 const TString &label, Float_t cal_err = 0,
                 Float_t gain_err = 0, Float_t bkg_err = 0,
                 Float_t cal_off = 0) {
  if (!r.valid || r.sig_channel.peaks.empty())
    return;
  const PeakFitResult &ge = r.sig_channel.peaks.back();
  result.ge_labels.push_back(label);
  result.ge_mus.push_back(ge.mu);
  result.ge_errs.push_back(ge.mu_error);
  result.ge_chi2s.push_back(r.sig_channel.reduced_chi2);
  result.ge_cal_errs.push_back(cal_err);
  result.ge_cal_offs.push_back(cal_off);
  result.ge_gain_errs.push_back(gain_err);
  result.ge_bkg_errs.push_back(bkg_err);
}

Float_t LoadCrystalLiveTime(TFile *file, Int_t crystal) {
  TParameter<Float_t> *p = static_cast<TParameter<Float_t> *>(
      file->Get(Form("LiveTime_Filtered_Crystal%d_s", crystal)));
  if (!p) {
    std::cerr << "WARNING: missing LiveTime_Filtered_Crystal" << crystal << "_s"
              << std::endl;
    return 0;
  }
  return p->GetVal();
}

TH1F *BuildRateHist(const TString &input_name, TF1 *cal, const TString &hname,
                    Float_t lo, Float_t hi) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    delete file;
    return nullptr;
  }
  Int_t nbins = (Int_t)std::lround((hi - lo) / Constants::BIN_WIDTH_KEV);
  TH1F *sum = new TH1F(hname, ";Deposited Energy [keV]; Rate [counts/s]", nbins,
                       lo, hi);
  sum->SetDirectory(nullptr);
  sum->Sumw2();

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t lt = LoadCrystalLiveTime(file, c);
    if (lt <= 0)
      continue;
    TH1F *ch = new TH1F(hname + Form("_c%d", c), "", nbins, lo, hi);
    ch->SetDirectory(nullptr);
    ch->Sumw2();
    Float_t energy = 0;
    tree->SetBranchAddress("energykeV", &energy);
    Int_t n = (Int_t)tree->GetEntries();
    for (Int_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      ch->Fill(cal->Eval(energy));
    }
    ch->Scale(1.0 / lt);
    sum->Add(ch);
    delete ch;
  }
  file->Close();
  delete file;
  return sum;
}

struct RateSubResult {
  Float_t mu = -1, mu_err = -1, chi2 = -1;
  Bool_t valid = kFALSE;
};

RateSubResult
FitRateSubtractedGe(const TString &sig_run, const TString &bkg_run, TF1 *cal,
                    const TString &tag, Bool_t interactive, Float_t lo = 62.0,
                    Float_t hi = 71.5, Bool_t linear_bkg = kFALSE,
                    Bool_t free_tail = kFALSE, Double_t tail_frac = -1.0,
                    Double_t tail_tau = -1.0, Double_t bkg_scale = 1.0) {
  RateSubResult out;

  TH1F *sig = BuildRateHist(sig_run, cal, "rate_sig_" + tag, lo, hi);
  TH1F *bkg = BuildRateHist(bkg_run, cal, "rate_bkg_" + tag, lo, hi);
  if (!sig || !bkg) {
    delete sig;
    delete bkg;
    return out;
  }
  TH1F *res = static_cast<TH1F *>(sig->Clone("rate_residual_" + tag));
  res->SetDirectory(nullptr);
  res->Add(bkg, -bkg_scale);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(res, kBlack);
  PlottingUtils::SaveFigure(canvas, "rate_residual_" + tag, "fits",
                            PlotSaveOptions::kLINEAR);

  // chi2 fit, not FittingUtils: its Poisson option assumes integer counts,
  // but this histogram holds counts/second and goes negative after
  // subtraction. BuildRateHist calls Sumw2(), so propagated errors exist.
  // Linear background: the 01/16 residual runs -2 to +2 counts/s across the
  // window and a flat model drags the centroid high.
  Double_t edge_lo = res->Integral(res->FindBin(lo + 0.2), res->FindBin(66.0));
  Double_t nedge = res->FindBin(66.0) - res->FindBin(lo + 0.2) + 1;
  Double_t bkg0 = (nedge > 0) ? edge_lo / nedge : 0.0;
  Int_t pk_bin = res->GetMaximumBin();
  Double_t pk_val = res->GetBinContent(pk_bin) - bkg0;
  Double_t sigma0 = 0.85;
  // Gaussian core PLUS a low-side exponential tail, as an exponentially
  // modified Gaussian (the exponential convolved with the resolution, which is
  // what the tail physically is).
  //
  // The tail is NOT optional. A bare Gaussian biases the centroid LOW on a
  // tailed peak, and the in-situ model this is compared against carries tails.
  // Fitting a bare Gaussian here is very likely why the rate-sub sat ~46 eV
  // below the in-situ once the in-situ shape was freed. Both methods must
  // describe the same lineshape or the method systematic is a comparison of
  // two different biases.
  //
  //   [0] area   [1] mu   [2] sigma   [3] bkg0   [4] bkg slope
  //   [5] tail fraction   [6] tail decay length tau [keV]
  const char *kEmgLowTail =
      "(1-[5])*[0]*TMath::Gaus(x,[1],[2],1)"
      " + [5]*[0]/(2*[6])*TMath::Exp((x-[1])/[6] + [2]*[2]/(2*[6]*[6]))"
      "   *TMath::Erfc((x-[1])/(TMath::Sqrt2()*[2]) + [2]/(TMath::Sqrt2()*[6]))"
      " + [3] + [4]*(x-68.75)";
  TF1 *model = new TF1("ratesub_model_" + tag, kEmgLowTail, lo, hi);
  model->SetParNames("Area", "Mu", "Sigma", "Bkg0", "BkgSlope", "TailFrac",
                     "TailTau");
  model->SetParameters(pk_val * sigma0 * 2.5066, E_GE_73M, sigma0, bkg0, 0.0,
                       0.20, 1.0);
  model->SetParLimits(1, lo + 1.0, hi - 1.0);
  model->SetParLimits(2, 0.2, 3.0);
  if (free_tail) {
    model->SetParLimits(5, 0.0, 0.6);
    model->SetParLimits(6, 0.2, 5.0);
  } else {
    model->FixParameter(5, tail_frac >= 0 ? tail_frac : RATESUB_TAIL_FRAC);
    model->FixParameter(6, tail_tau > 0 ? tail_tau : RATESUB_TAIL_TAU);
  }
  if (!linear_bkg)
    model->FixParameter(4, 0.0);
  TFitResultPtr fr = res->Fit(model, "S Q R");
  Bool_t ok = fr.Get() && fr->IsValid() && model->GetNDF() > 0;
  if (ok) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.valid = kTRUE;
    std::cout << "RATESUB-FIT " << tag << ": mu = " << std::fixed
              << std::setprecision(4) << out.mu << " +/- " << out.mu_err
              << " keV   sigma = " << model->GetParameter(2)
              << "   tailFrac = " << model->GetParameter(5)
              << "   tailTau = " << model->GetParameter(6)
              << "   bkg0 = " << model->GetParameter(3)
              << "   slope = " << model->GetParameter(4)
              << "   k = " << bkg_scale
              << "   chi2/ndf = " << std::setprecision(2) << out.chi2
              << std::endl;
  } else {
    std::cerr << "WARNING: rate-sub chi2 fit failed for " << tag << std::endl;
  }
  const Int_t kNPts = 400;
  Double_t xstep = (hi - lo) / (kNPts - 1);
  TGraph *total_graph = new TGraph(kNPts);
  TGraph *bkg_graph = new TGraph(kNPts);
  for (Int_t i = 0; i < kNPts; i++) {
    Double_t xv = lo + i * xstep;
    total_graph->SetPoint(i, xv, model->Eval(xv));
    bkg_graph->SetPoint(
        i, xv, model->GetParameter(3) + model->GetParameter(4) * (xv - 68.75));
  }
  bkg_graph->SetLineStyle(2);
  std::vector<TGraph *> components = {bkg_graph};

  res->GetXaxis()->SetRangeUser(lo, hi);
  Double_t ymin = 1e30, ymax = -1e30;
  for (Int_t b = res->FindBin(lo); b <= res->FindBin(hi - 1e-6); b++) {
    Double_t v = res->GetBinContent(b);
    Double_t e = res->GetBinError(b);
    if (v - e < ymin)
      ymin = v - e;
    if (v + e > ymax)
      ymax = v + e;
  }
  Double_t pad = 0.12 * (ymax - ymin);
  res->SetMinimum(ymin - pad);
  res->SetMaximum(ymax + pad);
  PlottingUtils::PlotFitWithResiduals(res, total_graph, components, lo, hi,
                                      "ratesub_fit_" + tag, "fits",
                                      "Rate-subtracted " + tag, kFALSE);
  delete sig;
  delete bkg;
  delete res;
  return out;
}

PairCalResult ProcessPbAnchoredPair(
    const TString &bkg_run, const TString &sig_run, const TString &date_label,
    const TString &pair_tag, Bool_t interactive, const TString &am_run,
    Bool_t is_reference, Float_t ref_p0, const CalFit *ref_cal = nullptr,
    const std::vector<Float_t> *ref_pb_mus = nullptr,
    const std::vector<Float_t> *ref_pb_errs = nullptr) {
  PairCalResult out;
  out.pair_tag = pair_tag;

  std::vector<Double_t> bkg_events = LoadFiltered(bkg_run);
  std::vector<Double_t> sig_events = LoadFiltered(sig_run);
  std::vector<Double_t> pb_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  TString cal_label = date_label + "_" + pair_tag;

  FitResult am;
  Bool_t have_am = kFALSE;
  if (is_reference && !am_run.IsNull()) {
    std::vector<Double_t> am_events = LoadFiltered(am_run);
    am = FitCalPeak(am_events, am_run, "Am_59.5keV", 51, 71, kFALSE, 0, 0,
                    interactive, kFALSE, USE_STEP_CAL_PEAKS);
    if (am.valid && !am.peaks.empty()) {
      const PeakFitResult &p = am.peaks.at(0);
      std::cout << "AM-SHAPE step=" << (USE_STEP_CAL_PEAKS ? "on " : "off")
                << "  chi2/ndf " << am.reduced_chi2 << "  sigma " << p.sigma
                << "  step " << p.step_amplitude << "  lowExpAmp "
                << p.low_exp_tail_amplitude << "  lowExpTau "
                << p.low_exp_tail_ratio << "  lowLinAmp "
                << p.low_lin_tail_amplitude << "  lowLinSlope "
                << p.low_lin_tail_slope << std::endl;
    }
    have_am = am.valid;
    if (!have_am)
      std::cerr << "WARNING: Am-241 fit failed for " << am_run
                << "; reference line will use Pb peaks only" << std::endl;
  }

  BkgGeSimResult pre =
      RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_run, sig_run,
                  "PreCal_" + pair_tag, interactive, kFALSE, nullptr, nullptr,
                  TAIL_SHARED, CAL_RANGE_LOW, CAL_RANGE_HIGH, kTRUE);
  if (!pre.valid)
    return out;

  if (pre.bkg_channel.peaks.size() >= 2) {
    Double_t mu_lo = pre.bkg_channel.peaks.at(0).mu;
    Double_t mu_hi = pre.bkg_channel.peaks.at(1).mu;
    Double_t g = (mu_lo + mu_hi) / (Double_t)(E_PB_KA2 + E_PB_KA1);
    Double_t delta_raw = (Double_t)(E_PB_KA1 - E_PB_KA2) * g;
    Double_t sig_raw = TMath::Sqrt(dE_PB_KA_SEP * dE_PB_KA_SEP +
                                   (delta_raw * 5e-4) * (delta_raw * 5e-4));
    std::cout << "PRECAL-SEP " << pair_tag << ": fitted " << std::fixed
              << std::setprecision(5) << (mu_hi - mu_lo) << "  expected "
              << delta_raw << " +/- " << sig_raw << "  (excess "
              << std::setprecision(1) << ((mu_hi - mu_lo) - delta_raw) * 1000
              << " eV)" << std::setprecision(4) << std::endl;
    if (CONSTRAIN_PRECAL_SEPARATION) {
      BkgGeSimResult pre2 =
          RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_run, sig_run,
                      "PreCal_" + pair_tag, interactive, kFALSE, nullptr,
                      nullptr, TAIL_SHARED, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                      kTRUE, delta_raw, sig_raw);
      if (pre2.valid)
        pre = pre2;
      else
        std::cerr << "WARNING: constrained pre-cal fit failed for " << pair_tag
                  << "; keeping the unconstrained result." << std::endl;
    }
  }

  for (size_t i = 0; i < pre.bkg_channel.peaks.size(); i++) {
    out.pb_mus_precal.push_back(pre.bkg_channel.peaks[i].mu);
    out.pb_mu_errs_precal.push_back(GuardedMuError(
        pair_tag + " Pb-Ka", pre.bkg_channel.peaks[i], pre.bkg_channel));
  }

  // Pb K-beta doublet (Kb1,3 group ~84.77 and Kb2 group ~87.32): same
  // FitDoublePeak method that seeds the K-alpha, run on the bkg spectrum in the
  // K-beta window. Fit as a doublet to cleanly separate the peaks, but anchor
  // only on Kb1 -- three Pb points (Ka2, Ka1, Kb1) spanning 12 keV per pair.
  // On the reference pair they are joined by Am; on the others they set the
  // gain.
  // link_sigma = kTRUE: the Pb K-beta1/K-beta2 groups share one detector
  // resolution, so tie their widths. Stops the weaker Kb2 from floating its
  // sigma and trading against the continuum, which would skew the Kb1 anchor.
  FitResult kb = FitCalPeak(bkg_events, bkg_run, "Pb_Kbeta", 81, 91, kTRUE,
                            E_PB_KB1, E_PB_KB2, interactive, kTRUE);

  if (kb.valid && !kb.peaks.empty()) {
    out.pb_mus_precal.push_back(kb.peaks.at(0).mu);
    out.pb_mu_errs_precal.push_back(
        GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb));
  }

  CalibrationData cal_data;
  AddCalPoint(cal_data, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
              GuardedMuError(pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0),
                             pre.bkg_channel),
              E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
  AddCalPoint(cal_data, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
              GuardedMuError(pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1),
                             pre.bkg_channel),
              E_PB_KA1, -1, dE_PB_KA);
  if (kb.valid)
    AddCalPoint(cal_data, pair_tag + " Pb-Kb1", kb.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb),
                E_PB_KB1, kb.reduced_chi2, dE_PB_KB1);
  else
    std::cerr << "WARNING: Pb K-beta fit failed for " << bkg_run
              << "; calibration falls back to K-alpha only" << std::endl;

  Bool_t transfer_mode =
      (ref_cal != nullptr && ref_pb_mus != nullptr && ref_pb_mus->size() >= 2 &&
       pre.bkg_channel.peaks.size() >= 2);
  Bool_t am_anchor_mode = (!transfer_mode && have_am);

  CalFit calfit;
  GainResult gain;
  AffineResult affine;

  if (transfer_mode) {
    std::vector<Float_t> mu_run = {(Float_t)pre.bkg_channel.peaks.at(0).mu,
                                   (Float_t)pre.bkg_channel.peaks.at(1).mu};
    std::vector<Float_t> err_run = {
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka2",
                                pre.bkg_channel.peaks.at(0), pre.bkg_channel),
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka1",
                                pre.bkg_channel.peaks.at(1), pre.bkg_channel)};
    if (kb.valid && !kb.peaks.empty()) {
      mu_run.push_back((Float_t)kb.peaks.at(0).mu);
      err_run.push_back(
          (Float_t)GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb));
    }
    std::vector<Float_t> er =
        ref_pb_errs ? *ref_pb_errs
                    : std::vector<Float_t>(ref_pb_mus->size(), 0.f);

    if (USE_AFFINE_TRANSFER)
      affine = PbAffineVsReference(mu_run, err_run, *ref_pb_mus, er);
    if (affine.valid) {
      calfit = AffineScaleCal(*ref_cal, affine, cal_label);
      std::cout << "AFFINE-TRANSFER " << cal_label << ": a = " << std::fixed
                << std::setprecision(6) << affine.a << " +/- "
                << std::sqrt(std::max(0.0, affine.var_a))
                << "   b = " << affine.b << " +/- "
                << std::sqrt(std::max(0.0, affine.var_b))
                << "  (n = " << affine.n
                << ", chi2/ndf = " << std::setprecision(2) << affine.chi2_ndf
                << ")" << std::endl;
    } else {
      gain = PbGainVsReferenceWithError(mu_run, err_run, *ref_pb_mus, er);
      calfit = GainScaleCal(*ref_cal, gain.g, cal_label);
      std::cout << "GAIN-TRANSFER (fallback) " << cal_label
                << ": g = " << std::fixed << std::setprecision(6) << gain.g
                << " +/- " << gain.sigma_g << "  (n = " << gain.n
                << ", chi2/ndf = " << std::setprecision(2) << gain.chi2_ndf
                << ")" << std::endl;
    }
  } else if (am_anchor_mode) {
    CalibrationData cal_am;
    AddCalPoint(cal_am, pair_tag + " Am-241 59.5", am.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Am-241", am.peaks.at(0), am),
                E_AM241, am.reduced_chi2, dE_AM241);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Ka2",
                               pre.bkg_channel.peaks.at(0), pre.bkg_channel),
                E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
                GuardedMuError(pair_tag + " Pb-Ka1",
                               pre.bkg_channel.peaks.at(1), pre.bkg_channel),
                E_PB_KA1, -1, dE_PB_KA);
    if (kb.valid)
      AddCalPoint(cal_am, pair_tag + " Pb-Kb1", kb.peaks.at(0).mu,
                  GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb),
                  E_PB_KB1, kb.reduced_chi2, dE_PB_KB1);
    PrintCalSummary(cal_am, cal_label);
    calfit = CreateAndSavePol1Cal(cal_am, cal_label, kFALSE, 0, REF_CAL_DEGREE);
  } else {
    PrintCalSummary(cal_data, cal_label);
    calfit = CreateAndSavePol1Cal(cal_data, cal_label, kFALSE, 0);
  }
  TF1 *cal = calfit.func;

  CalFit pb_probe = calfit;
  if (have_am)
    pb_probe = CreateAndSavePol1Cal(cal_data, cal_label + "_PbOnly", kFALSE, 0);

  if (have_am) {
    out.ref_calfit = calfit;
    out.has_ref_cal = kTRUE;
  }
  out.cal_func = cal;
  if (!pre.sig_channel.peaks.empty()) {
    Double_t mu_ge = pre.sig_channel.peaks.back().mu;
    out.ge_cal_err = CalEnergyError(calfit, mu_ge);
    out.ge_cal_off = CalOffsetError(calfit);

    if (transfer_mode) {
      out.ge_gain_err =
          affine.valid
              ? AffineTransferEnergyError(*ref_cal, mu_ge, affine)
              : GainTransferEnergyError(*ref_cal, mu_ge, gain.g, gain.sigma_g);
      out.ge_cal_err = std::hypot(out.ge_cal_err, out.ge_gain_err);

      std::cout << "GAIN-CAL " << cal_label << ": transfer term " << std::fixed
                << std::setprecision(1) << out.ge_gain_err * 1000.0
                << " eV   total cal " << out.ge_cal_err * 1000.0 << " eV"
                << std::endl;
    }
  }

  if (pre.bkg_channel.peaks.size() >= 2 && !pre.sig_channel.peaks.empty()) {
    Double_t mu_pb1 = pre.bkg_channel.peaks[0].mu;
    Double_t mu_pb2 = pre.bkg_channel.peaks[1].mu;
    Double_t denom = mu_pb2 - mu_pb1;
    if (std::fabs(denom) > 0) {
      Double_t ge_raw = pre.sig_channel.peaks.back().mu;
      Double_t r_ge_pb = (ge_raw - mu_pb1) / denom;
      std::cout << "R_GEPB " << cal_label << ": raw Ka2 " << std::fixed
                << std::setprecision(4) << mu_pb1 << "  raw Ka1 " << mu_pb2
                << "  raw Ge " << ge_raw << "  R = " << std::setprecision(6)
                << r_ge_pb << std::endl;
    }
  }

  // Cross-checks on the pair that also has Am. All lines are evaluated at the
  // SAME raw Ge centroid, so the fit cancels and only the calibration differs
  // -- this is the one place the two effects are cleanly separable.
  //
  // The nominal line already contains Am, so the LINEARITY test is scored
  // against a separately rebuilt Pb-only line: Am must play no part in the
  // line whose Am prediction is being tested, or the check is circular.
  if (have_am && !pre.sig_channel.peaks.empty()) {
    CalibrationData cal_am;
    AddCalPoint(cal_am, pair_tag + " Am-241 59.5", am.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Am-241", am.peaks.at(0), am),
                E_AM241, am.reduced_chi2, dE_AM241);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Ka2",
                               pre.bkg_channel.peaks.at(0), pre.bkg_channel),
                E_PB_KA2, -1, dE_PB_KA);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
                GuardedMuError(pair_tag + " Pb-Ka1",
                               pre.bkg_channel.peaks.at(1), pre.bkg_channel),
                E_PB_KA1, -1, dE_PB_KA);
    CalFit cal_check =
        CreateAndSavePol1Cal(cal_am, cal_label + "_AmCheck", kFALSE, 0);
    Double_t ge_pre = pre.sig_channel.peaks.back().mu;
    Double_t e_nominal = cal->Eval(ge_pre);
    Double_t e_pbonly = pb_probe.func->Eval(ge_pre);
    Double_t e_amka = cal_check.func->Eval(ge_pre);
    std::cout << "VALIDATION " << pair_tag << " (Ge precal " << std::fixed
              << std::setprecision(4) << ge_pre << "):  nominal -> "
              << e_nominal << " keV   Pb-only -> " << e_pbonly
              << " keV   Am+Ka -> " << e_amka << " keV   (nominal - Pb-only) "
              << std::setprecision(1) << (e_nominal - e_pbonly) * 1000.0
              << " eV   (nominal - Am+Ka) " << (e_nominal - e_amka) * 1000.0
              << " eV" << std::endl;

    Double_t am_raw = am.peaks.at(0).mu;
    Double_t am_pred = pb_probe.func->Eval(am_raw);
    Double_t sig_cal = CalEnergyError(pb_probe, am_raw);
    Double_t sig_amfit = pb_probe.p1 * am.peaks.at(0).mu_error;
    Double_t sig_pred = std::sqrt(sig_cal * sig_cal + sig_amfit * sig_amfit);
    Double_t resid = am_pred - E_AM241;
    std::cout << "LINEARITY " << pair_tag
              << ": Pb-only line predicts Am = " << std::fixed
              << std::setprecision(4) << am_pred << " +/- " << sig_pred
              << " keV   (truth " << E_AM241 << ", residual "
              << std::setprecision(1) << resid * 1000.0
              << " eV = " << std::setprecision(2)
              << (sig_pred > 0 ? resid / sig_pred : 0) << " sigma)"
              << std::endl;
  }

  std::vector<TString> pair_runs = {bkg_run, sig_run};
  if (have_am)
    pair_runs.push_back(am_run);
  ApplyPol1Cal(pair_runs, cal, cal_label);

  std::vector<Double_t> bc = LoadCalibrated(bkg_run);
  std::vector<Double_t> sc = LoadCalibrated(sig_run);

  Double_t cal_p0 = cal->GetParameter(0);
  Double_t cal_p1 = cal->GetParameter(1);
  FitResult bkg_seed_post = pre.bkg_channel;
  FitResult sig_seed_post = pre.sig_channel;
  if (cal_p1 > 0) {
    Double_t s_b = pre.bkg_channel.lin_bkg_slope;
    Double_t s_s = pre.sig_channel.lin_bkg_slope;
    bkg_seed_post.lin_bkg_slope = s_b / (cal_p1 - s_b * cal_p0);
    sig_seed_post.lin_bkg_slope = s_s / (cal_p1 - s_s * cal_p0);
    bkg_seed_post.bkg_constant = pre.bkg_channel.bkg_constant / cal_p1;
    sig_seed_post.bkg_constant = pre.sig_channel.bkg_constant / cal_p1;
  }
  // Map the seed peak mus through the calibration so the Pb and Ge peaks start
  // at their calibrated positions, and convert the SHAPE into calibrated units
  // as well: sigma and the tail lengths were fitted in RAW units and carry the
  // local derivative f'(mu) ~ 0.989-0.995, so handing them over unconverted
  // makes every width ~1% too large (sigma wide by 7-10 eV).
  //
  // Scaling under E = f(x) with s = f'(mu), evaluated per peak because the
  // reference calibration is a pol2 and the derivative varies across the range:
  //   sigma, *TailRatio (lengths in x)  -> multiply by s
  //   LowLinTailSlope   (per unit x)    -> divide by s
  //   *TailAmplitude, step (ratios)     -> unchanged
  auto convert_peaks = [&cal](const FitResult &src, FitResult &dst) {
    for (size_t i = 0; i < dst.peaks.size() && i < src.peaks.size(); i++) {
      Double_t mu_raw = src.peaks[i].mu;
      Double_t s = cal->Derivative(mu_raw);
      if (!(s > 0))
        s = 1.0;
      dst.peaks[i].mu = cal->Eval(mu_raw);
      dst.peaks[i].sigma = src.peaks[i].sigma * s;
      dst.peaks[i].low_exp_tail_ratio = src.peaks[i].low_exp_tail_ratio * s;
      dst.peaks[i].high_exp_tail_ratio = src.peaks[i].high_exp_tail_ratio * s;
      dst.peaks[i].low_lin_tail_slope = src.peaks[i].low_lin_tail_slope / s;
    }
  };
  convert_peaks(pre.bkg_channel, bkg_seed_post);
  convert_peaks(pre.sig_channel, sig_seed_post);

  // Reports only. sigma and the tail parameters are fitted in RAW units and
  // seeded into a fit of CALIBRATED data, so anything carrying units of energy
  // is wrong by s = f'(mu), the local derivative: sigma and tail ratios scale
  // by s, low_lin_tail_slope by 1/s, amplitudes are dimensionless.
  if (!pre.sig_channel.peaks.empty()) {
    const PeakFitResult &ge_pre = pre.sig_channel.peaks.back();
    Double_t s = cal->Derivative(ge_pre.mu);
    std::cout << "UNIT-CHECK " << cal_label
              << " (Ge peak, f'(mu) = " << std::fixed << std::setprecision(6)
              << s << "):" << std::endl;
    std::cout << "  sigma            raw " << std::setprecision(5)
              << ge_pre.sigma << "   used " << ge_pre.sigma << "   should be "
              << ge_pre.sigma * s << "   (delta "
              << (ge_pre.sigma - ge_pre.sigma * s) * 1000.0 << " eV)"
              << std::endl;
    std::cout << "  lowExpTailRatio  raw " << ge_pre.low_exp_tail_ratio
              << "   used " << ge_pre.low_exp_tail_ratio << "   should be "
              << ge_pre.low_exp_tail_ratio * s << std::endl;
    std::cout << "  lowLinTailSlope  raw " << ge_pre.low_lin_tail_slope
              << "   used " << ge_pre.low_lin_tail_slope << "   should be "
              << (s != 0 ? ge_pre.low_lin_tail_slope / s : 0.0) << std::endl;
    std::cout << "  highExpTailRatio raw " << ge_pre.high_exp_tail_ratio
              << "   used " << ge_pre.high_exp_tail_ratio << "   should be "
              << ge_pre.high_exp_tail_ratio * s << std::endl;
  }

  // Pb peaks left FREE post-cal (lock_bkg_peaks = kFALSE): where the calibrated
  // Pb doublet actually lands is a blind check at the energy of interest. The
  // whole chain (precal -> cal -> postcal) uses ONE tail model (TAIL_SHARED),
  // so this run is internally consistent; the alternative tail model is a
  // separate macro and the lineshape systematic is taken in the combiner.
  // force_sim_fit = kTRUE, same as the pre-cal. The saved .simroofits carries
  // the hand-tuned starting point a cold start cannot reach; SetRefitAfterLoad
  // makes it SEED a real fit rather than be returned as the result, which once
  // produced post-cal rows holding day-one constants at chi2 7-9.
  out.postcal = RunBkgGeSim(
      bc, sc, pb_mus, bkg_run + "_postcal_" + pair_tag,
      sig_run + "_postcal_" + pair_tag, "PostCal_" + pair_tag, kTRUE,
      LOCK_POSTCAL_PB, &bkg_seed_post, &sig_seed_post, TAIL_SHARED,
      CAL_RANGE_LOW, CAL_RANGE_HIGH, kTRUE,
      CONSTRAIN_PB_KA_SEPARATION ? (Double_t)(E_PB_KA1 - E_PB_KA2) : -1.0,
      dE_PB_KA_SEP);

  if (out.postcal.valid) {
    std::cout << "Post-cal FREE Pb mu for " << cal_label << ":" << std::endl;
    for (size_t i = 0;
         i < out.postcal.bkg_channel.peaks.size() && i < pb_mus.size(); i++)
      std::cout << "  Pb peak " << i << ": " << std::fixed
                << std::setprecision(4) << out.postcal.bkg_channel.peaks[i].mu
                << " keV (truth " << pb_mus[i] << ", shift "
                << (out.postcal.bkg_channel.peaks[i].mu - pb_mus[i]) << ")"
                << std::endl;

    if (!pre.sig_channel.peaks.empty()) {
      Double_t ge_precal_mapped = cal->Eval(pre.sig_channel.peaks.back().mu);
      Double_t ge_postcal = out.postcal.sig_channel.peaks.back().mu;
      out.ge_bkg_err = std::fabs(ge_precal_mapped - ge_postcal);
      std::cout << "BACKGROUND " << cal_label
                << ": cal(precal Ge)=" << std::fixed << std::setprecision(4)
                << ge_precal_mapped << "  postcal Ge=" << ge_postcal
                << "  |diff|=" << std::setprecision(1)
                << out.ge_bkg_err * 1000.0 << " eV" << std::endl;
    }
  }
  return out;
}

void AppendPairToDay(DayResult &result, const PairCalResult &p,
                     const TString &date_label, const TString &ge_label) {
  if (p.cal_func) {
    result.cal_labels.push_back(date_label + " " + p.pair_tag);
    result.cal_funcs.push_back(p.cal_func);
  }
  AddGeResult(result, p.postcal, ge_label, p.ge_cal_err, p.ge_gain_err,
              p.ge_bkg_err, p.ge_cal_off);
}

void PrintPrecalPbComparison(const TString &date_label,
                             const std::vector<PairCalResult *> &pairs) {
  std::cout << "Pre-cal Pb mu comparison for " << date_label << std::endl;
  std::cout << std::left << std::setw(14) << "Pair" << std::right
            << std::setw(14) << "Pb-Ka2 mu" << std::setw(14) << "shift vs true"
            << std::setw(14) << "Pb-Ka1 mu" << std::setw(14) << "shift vs true"
            << std::endl;
  for (size_t i = 0; i < pairs.size(); i++) {
    const PairCalResult *p = pairs[i];
    if (!p || p->pb_mus_precal.size() < 2)
      continue;
    Float_t mu1 = p->pb_mus_precal[0];
    Float_t mu2 = p->pb_mus_precal[1];
    std::cout << std::left << std::setw(14) << p->pair_tag << std::right
              << std::fixed << std::setprecision(4) << std::setw(14) << mu1
              << std::setw(14) << (mu1 - E_PB_KA2) << std::setw(14) << mu2
              << std::setw(14) << (mu2 - E_PB_KA1) << std::endl;
  }
}

DayResult ProcessLineCalDay(const LineCalConfig &cfg, Bool_t interactive) {
  CalibrationData cal_data;

  std::vector<Double_t> am_events = LoadFiltered(cfg.am_run);
  FitResult am = FitCalPeak(am_events, cfg.am_run, "Am_59.5keV", cfg.am_lo,
                            cfg.am_hi, kFALSE, 0, 0, interactive);
  if (am.valid)
    AddCalPoint(cal_data, "Am-241 59.5", am.peaks.at(0).mu,
                am.peaks.at(0).mu_error, E_AM241, am.reduced_chi2, dE_AM241);

  const Bool_t USE_BA133_53 = kFALSE;

  std::vector<Double_t> ba_events = LoadFiltered(cfg.ba_run);
  if (USE_BA133_53) {
    FitResult ba53 = FitCalPeak(ba_events, cfg.ba_run, "Ba_53.16keV", 48, 58,
                                kFALSE, 0, 0, interactive);
    if (ba53.valid)
      AddCalPoint(cal_data, "Ba-133 53.16", ba53.peaks.at(0).mu,
                  ba53.peaks.at(0).mu_error, E_BA133_53, ba53.reduced_chi2,
                  dE_LINE);
  }

  FitResult ba81 =
      FitCalPeak(ba_events, cfg.ba_run, "Ba_80.98keV", 75, 90, kTRUE,
                 E_BA133_79, E_BA133_81, interactive, kTRUE);
  if (ba81.valid)
    AddCalPoint(cal_data, "Ba-133 80.98", ba81.peaks.at(1).mu,
                ba81.peaks.at(1).mu_error, E_BA133_81, ba81.reduced_chi2,
                dE_LINE);

  PrintCalSummary(cal_data, cfg.date_label);
  CalFit calfit = CreateAndSavePol1Cal(cal_data, cfg.date_label);
  TF1 *cal = calfit.func;
  Double_t mu_at_ge =
      (calfit.p1 != 0) ? (E_GE_73M - calfit.p0) / calfit.p1 : E_GE_73M;
  Double_t ge_cal_err = CalEnergyError(calfit, mu_at_ge);
  Double_t ge_cal_off = CalOffsetError(calfit);
  ApplyPol1Cal(cfg.day_datasets, cal, cfg.date_label);

  DayResult result;
  result.cal_labels.push_back(cfg.date_label);
  result.cal_funcs.push_back(cal);

  RateSubResult rs =
      FitRateSubtractedGe(cfg.postcal_sig, cfg.postcal_bkg, cal, cfg.date_label,
                          kFALSE, cfg.rs_lo, cfg.rs_hi, cfg.rs_linear_bkg);
  if (rs.valid) {
    result.rs_labels.push_back(cfg.ge_label + " [rate-sub]");
    result.rs_mus.push_back(rs.mu);
    result.rs_errs.push_back(rs.mu_err);
    result.rs_chi2s.push_back(rs.chi2);
    result.rs_cal_errs.push_back(ge_cal_err);
    result.rs_cal_offs.push_back(ge_cal_off);
    result.rs_gain_errs.push_back(0.0f);
    std::cout << "Rate-sub Ge mu for " << cfg.date_label << ": " << std::fixed
              << std::setprecision(4) << rs.mu << " +/- " << rs.mu_err
              << " keV (chi2/ndf = " << std::setprecision(3) << rs.chi2 << ")"
              << std::endl;
  }
  return result;
}

struct ShieldPair {
  TString bkg_run;
  TString sig_run;
  TString tag;
  TString ge_label;
};

DayResult ProcessDay_20260113(Bool_t interactive, Float_t &ref_p0_out,
                              PairCalResult &ref_pair_out) {
  DayResult result;

  std::vector<ShieldPair> pairs;
  ShieldPair cu;
  cu.bkg_run = Constants::CUSHIELDBACKGROUND_10PERCENT_20260113;
  cu.sig_run = Constants::CUSHIELDSIGNAL_10PERCENT_20260113;
  cu.tag = "CuShield10";
  cu.ge_label = "Cu Shield Signal 10% (01/13)";
  pairs.push_back(cu);
  ShieldPair cd10;
  cd10.bkg_run = Constants::CDSHIELDBACKGROUND_10PERCENT_20260113;
  cd10.sig_run = Constants::CDSHIELDSIGNAL_10PERCENT_20260113;
  cd10.tag = "CdShield10";
  cd10.ge_label = "Cd Shield Signal 10% (01/13)";
  pairs.push_back(cd10);
  ShieldPair cd25;
  cd25.bkg_run = Constants::CDSHIELDBACKGROUND_25PERCENT_20260113;
  cd25.sig_run = Constants::CDSHIELDSIGNAL_25PERCENT_20260113;
  cd25.tag = "CdShield25";
  cd25.ge_label = "Cd Shield Signal 25% (01/13)";
  pairs.push_back(cd25);

  TString am_run = Constants::POSTREACTOR_AM241_20260113;
  Float_t am_time = LoadRunStartTime(am_run);
  Int_t ref_idx = 0;
  Float_t best_dt = 1e30f;
  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    Float_t dt = TMath::Abs(LoadRunStartTime(pairs[i].bkg_run) - am_time);
    if (dt < best_dt) {
      best_dt = dt;
      ref_idx = i;
    }
  }
  std::cout << "Reference pair for 13/14 drift: " << pairs[ref_idx].tag
            << " (Am from " << am_run << ")" << std::endl;

  std::vector<PairCalResult> results(pairs.size());

  results[ref_idx] = ProcessPbAnchoredPair(
      pairs[ref_idx].bkg_run, pairs[ref_idx].sig_run, "20260113",
      pairs[ref_idx].tag, interactive, am_run, kTRUE, 0);
  ref_p0_out = results[ref_idx].cal_func
                   ? results[ref_idx].cal_func->GetParameter(0)
                   : 0;

  const CalFit *ref_cal =
      results[ref_idx].has_ref_cal ? &results[ref_idx].ref_calfit : nullptr;
  const std::vector<Float_t> *ref_mus = &results[ref_idx].pb_mus_precal;
  const std::vector<Float_t> *ref_errs = &results[ref_idx].pb_mu_errs_precal;

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    if (i == ref_idx)
      continue;
    results[i] = ProcessPbAnchoredPair(
        pairs[i].bkg_run, pairs[i].sig_run, "20260113", pairs[i].tag,
        interactive, "", kFALSE, ref_p0_out, ref_cal, ref_mus, ref_errs);
  }

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++)
    AppendPairToDay(result, results[i], "20260113", pairs[i].ge_label);

  if (RATESUB_ALL_DAYS)
    for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
      if (!results[i].cal_func)
        continue;
      RateSubResult rs = FitRateSubtractedGe(
          pairs[i].sig_run, pairs[i].bkg_run, results[i].cal_func,
          "20260113_" + pairs[i].tag, kFALSE, 62.0, 71.5, kTRUE);
      if (!rs.valid)
        continue;
      result.rs_labels.push_back(pairs[i].ge_label + " [rate-sub]");
      result.rs_mus.push_back(rs.mu);
      result.rs_errs.push_back(rs.mu_err);
      result.rs_chi2s.push_back(rs.chi2);
      result.rs_cal_errs.push_back(results[i].ge_cal_err);
      result.rs_cal_offs.push_back(results[i].ge_cal_off);
      result.rs_gain_errs.push_back(results[i].ge_gain_err);
      std::cout << "Rate-sub Ge mu for 20260113_" << pairs[i].tag << ": "
                << std::fixed << std::setprecision(4) << rs.mu << " +/- "
                << rs.mu_err << " keV (chi2/ndf = " << std::setprecision(3)
                << rs.chi2 << ")" << std::endl;
    }

  std::vector<PairCalResult *> ptrs;
  for (Int_t i = 0; i < (Int_t)pairs.size(); i++)
    ptrs.push_back(&results[i]);
  PrintPrecalPbComparison("20260113", ptrs);
  ref_pair_out = results[ref_idx];
  return result;
}

DayResult ProcessDay_20260114(Bool_t interactive, Float_t ref_p0,
                              const PairCalResult &ref_pair) {
  DayResult result;
  const CalFit *ref_cal = ref_pair.has_ref_cal ? &ref_pair.ref_calfit : nullptr;
  PairCalResult cu = ProcessPbAnchoredPair(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, "20260114", "CuShield10",
      interactive, "", kFALSE, ref_p0, ref_cal, &ref_pair.pb_mus_precal,
      &ref_pair.pb_mu_errs_precal);
  AppendPairToDay(result, cu, "20260114", "Cu Shield Signal 10% (01/14)");

  if (cu.cal_func) {
    if (RATESUB_TAIL_SCAN) {
      const Double_t scan_k[] = {0.94, 0.97, 1.00, 1.03, 1.06};
      for (size_t k = 0; k < sizeof(scan_k) / sizeof(scan_k[0]); k++)
        FitRateSubtractedGe(
            Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, cu.cal_func,
            TString::Format("20260114_k%03d", (Int_t)(scan_k[k] * 100)), kFALSE,
            62.0, 71.5, kTRUE, kFALSE, 0.0, 1.0, scan_k[k]);
      const Double_t scan_fr[] = {0.0, 0.05, 0.10, 0.20, 0.30};
      for (size_t k = 0; k < sizeof(scan_fr) / sizeof(scan_fr[0]); k++)
        FitRateSubtractedGe(
            Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, cu.cal_func,
            TString::Format("20260114_fix%03d", (Int_t)(scan_fr[k] * 100)),
            kFALSE, 62.0, 71.5, kTRUE, kFALSE, scan_fr[k], 1.0);
      FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                          Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                          cu.cal_func, "20260114_freetail_flat", kFALSE, 62.0,
                          71.5, kFALSE, kTRUE);
      FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                          Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                          cu.cal_func, "20260114_freetail_lin", kFALSE, 62.0,
                          71.5, kTRUE, kTRUE);
    }
    RateSubResult rs =
        FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                            cu.cal_func, "20260114", kFALSE, 62.0, 71.5, kTRUE);
    if (rs.valid) {
      result.rs_labels.push_back("Cu Shield Signal 10% (01/14) [rate-sub]");
      result.rs_mus.push_back(rs.mu);
      result.rs_errs.push_back(rs.mu_err);
      result.rs_chi2s.push_back(rs.chi2);
      result.rs_cal_errs.push_back(cu.ge_cal_err);
      result.rs_cal_offs.push_back(cu.ge_cal_off);
      result.rs_gain_errs.push_back(cu.ge_gain_err);
      std::cout << "Rate-sub Ge mu for 20260114: " << std::fixed
                << std::setprecision(4) << rs.mu << " +/- " << rs.mu_err
                << " keV (chi2/ndf = " << std::setprecision(3) << rs.chi2 << ")"
                << std::endl;
    }
  }
  return result;
}

DayResult ProcessDay_20260115(Bool_t interactive) {
  LineCalConfig cfg;
  cfg.date_label = "20260115";
  cfg.am_run = Constants::POSTREACTOR_AM241_20260115;
  cfg.am_lo = 51;
  cfg.am_hi = 71;
  cfg.ba_run = Constants::POSTREACTOR_BA133_20260115;
  cfg.day_datasets = {Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::POSTREACTOR_AM241_20260115,
                      Constants::POSTREACTOR_BA133_20260115,
                      Constants::SHUTTERCLOSED_20260115};
  cfg.postcal_bkg = Constants::NOSHIELDBACKGROUND_5PERCENT_20260115;
  cfg.postcal_sig = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;
  cfg.ge_label = "No Shield Signal 5% (01/15)";
  return ProcessLineCalDay(cfg, interactive);
}

DayResult ProcessDay_20260116(Bool_t interactive) {
  LineCalConfig cfg;
  cfg.date_label = "20260116";
  cfg.am_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.am_lo = 55;
  cfg.am_hi = 70;
  cfg.ba_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.day_datasets = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  cfg.postcal_bkg =
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116;
  cfg.postcal_sig = Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116;
  cfg.ge_label = "Graphite Castle Signal 10% (01/16)";
  cfg.rs_lo = 62.0;
  cfg.rs_hi = 71.5;
  return ProcessLineCalDay(cfg, interactive);
}

void PrintDayCals(const DayResult &d) {
  for (size_t i = 0; i < d.cal_funcs.size(); i++) {
    TF1 *c = d.cal_funcs[i];
    if (!c)
      continue;
    std::cout << "Cal pol1 " << d.cal_labels[i] << ": p0 = " << std::fixed
              << std::setprecision(6) << c->GetParameter(0) << " +/- "
              << c->GetParError(0) << ", p1 = " << c->GetParameter(1) << " +/- "
              << c->GetParError(1) << std::endl;
  }
}

void CalibrationLow() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/calibration_low_" +
                                    RESULT_TAG,
                                project_root + "/root_files");
  Bool_t interactive = kTRUE;

  Float_t ref_p0 = 0;
  PairCalResult ref_pair;
  DayResult d13 = ProcessDay_20260113(interactive, ref_p0, ref_pair);
  DayResult d14 = ProcessDay_20260114(interactive, ref_p0, ref_pair);
  DayResult d15 = ProcessDay_20260115(interactive);
  DayResult d16 = ProcessDay_20260116(interactive);

  PrintDayCals(d13);
  PrintDayCals(d14);
  PrintDayCals(d15);
  PrintDayCals(d16);

  // 15th/16th: line-cal (Am + Ba-133) days measured by rate subtraction.
  // Independent of the 13th/14th in both calibration source and background
  // treatment, which is what makes the method systematic meaningful.
  // One .result row per run; no budget math here.
  // 01/14 is the only cross-day gain transfer (no Am on the 14th) and has
  // the worst chi2 of any run; its in-situ row fails the combiner's quality
  // cut, but its rate-sub row is kept.
  const Bool_t INCLUDE_20260114 = kTRUE;

  std::vector<DayResult *> day_vec = {&d13};
  if (INCLUDE_20260114)
    day_vec.push_back(&d14);
  day_vec.push_back(&d15);
  day_vec.push_back(&d16);
  DayResult **days = day_vec.data();
  const Int_t n_days = (Int_t)day_vec.size();

  const TString results_dir = project_root + "/results";
  gSystem->mkdir(results_dir, kTRUE);
  const TString result_path =
      results_dir + "/ge_" + RESULT_TAG + RESULT_VARIANT + ".result";
  std::ofstream out(result_path.Data());
  out << "# Ge-73m gamma energy results   tag=" << RESULT_TAG << std::endl;
  out << "# method  ge_mu  fit_err  cal_err  bkg_err  chi2  gain_err  cal_off  "
         "label"
      << std::endl;
  out << std::fixed << std::setprecision(6);

  auto sanitize = [](TString s) {
    s.ReplaceAll(" ", "_");
    return s;
  };

  std::cout << std::endl;
  std::cout << "Method 1 -- in-situ simultaneous fit (Ge Peak mu, post-cal):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->ge_labels.size(); i++) {
      Float_t ce =
          i < days[d]->ge_cal_errs.size() ? days[d]->ge_cal_errs[i] : 0;
      Float_t be =
          i < days[d]->ge_bkg_errs.size() ? days[d]->ge_bkg_errs[i] : 0;
      Float_t gz =
          i < days[d]->ge_gain_errs.size() ? days[d]->ge_gain_errs[i] : 0;
      Float_t co =
          i < days[d]->ge_cal_offs.size() ? days[d]->ge_cal_offs[i] : 0;
      out << "insitu  " << days[d]->ge_mus[i] << "  " << days[d]->ge_errs[i]
          << "  " << ce << "  " << be << "  " << days[d]->ge_chi2s[i] << "  "
          << gz << "  " << co << "  " << sanitize(days[d]->ge_labels[i])
          << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->ge_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->ge_mus[i]
                << " +/- " << days[d]->ge_errs[i] << " (fit) +/- " << ce
                << " (cal) keV (chi2/ndf = " << std::setprecision(3)
                << days[d]->ge_chi2s[i] << ")" << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "Method 2 -- live-time rate subtraction (Ge Peak mu):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->rs_labels.size(); i++) {
      Float_t ce =
          i < days[d]->rs_cal_errs.size() ? days[d]->rs_cal_errs[i] : 0;
      out << "ratesub  " << days[d]->rs_mus[i] << "  " << days[d]->rs_errs[i]
          << "  " << ce << "  " << 0.0 << "  " << days[d]->rs_chi2s[i] << "  "
          << (i < days[d]->rs_gain_errs.size() ? days[d]->rs_gain_errs[i]
                                               : 0.0f)
          << "  "
          << (i < days[d]->rs_cal_offs.size() ? days[d]->rs_cal_offs[i] : 0.0f)
          << "  " << sanitize(days[d]->rs_labels[i]) << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->rs_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->rs_mus[i]
                << " +/- " << days[d]->rs_errs[i] << " (fit) +/- " << ce
                << " (cal) keV" << std::endl;
    }
  }
  out.close();
  std::cout << std::endl;
  std::cout << "Wrote " << result_path << std::endl;
  std::cout << "Run CombineGeResult.cpp for the budget." << std::endl;
}

void RunFitCSVExport() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/calibration_low_" +
                                    RESULT_TAG,
                                project_root + "/root_files");
  gSystem->mkdir(project_root + "/results", kTRUE);

  G_CSV_PREFIX = project_root + "/results/spectrum_CuShield10_20260113";
  Float_t ref_p0 = 0;
  ProcessPbAnchoredPair(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                        Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                        "20260113", "CuShield10", kTRUE,
                        Constants::POSTREACTOR_AM241_20260113, kTRUE, ref_p0);
  G_CSV_PREFIX = "";
  std::cout << std::endl;
  std::cout << "CSV export done. The sig file is the Pb-Ka + Ge spectrum; "
               "the bkg file is the Pb-only background channel."
            << std::endl;
}
