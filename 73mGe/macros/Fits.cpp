#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iomanip>
#include <vector>

const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;
const Float_t E_GE_73M = 68.752;
const Float_t E_BKG_73KEV = 73.0;

std::vector<Double_t> LoadEvents(const TString input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("calibrated/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open calibrated/" << input_name << ".root"
              << std::endl;
    return events;
  }
  TTree *tree = static_cast<TTree *>(file->Get("calibrated_tree"));
  if (!tree) {
    std::cerr << "ERROR: Cannot find calibrated_tree in calibrated/"
              << input_name << ".root" << std::endl;
    file->Close();
    delete file;
    return events;
  }
  events = RooFitUtils::LoadEventsFromTree(tree, "depositedEnergykeV");
  file->Close();
  delete file;
  return events;
}

struct BkgGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid;
};

BkgGeSimResult FitBkgGeSimultaneous(
    const TString bkg_input, const TString sig_input, const Float_t bkg_lo,
    const Float_t bkg_hi, const Float_t sig_lo, const Float_t sig_hi,
    const std::vector<Double_t> &bkg_peak_mus, const Double_t ge_mu_init,
    const Bool_t bkg_flat, const Bool_t sig_flat, const Bool_t interactive) {
  BkgGeSimResult out;
  out.valid = kFALSE;

  Int_t n_bkg_peaks = (Int_t)bkg_peak_mus.size();
  if (n_bkg_peaks < 1 || n_bkg_peaks > 2) {
    std::cerr << "ERROR: bkg_peak_mus must have 1 or 2 entries" << std::endl;
    return out;
  }

  std::vector<Double_t> bkg_events = LoadEvents(bkg_input);
  std::vector<Double_t> sig_events = LoadEvents(sig_input);
  if (bkg_events.empty() || sig_events.empty())
    return out;

  Bool_t use_step = kFALSE;
  Bool_t use_low_exp = kTRUE;
  Bool_t use_low_lin = kTRUE;
  Bool_t use_high_exp = kTRUE;

  RooFitUtils bkg_fitter(bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV,
                         bkg_flat, use_step, use_low_exp, use_low_lin,
                         use_high_exp);
  if (interactive)
    bkg_fitter.SetInteractive();

  FitResult seed;
  TString seed_label = "BkgSeed";
  if (n_bkg_peaks == 1) {
    seed = bkg_fitter.FitSinglePeak(bkg_input, seed_label);
  } else {
    seed = bkg_fitter.FitDoublePeak(bkg_input, seed_label, bkg_peak_mus[0],
                                    bkg_peak_mus[1]);
  }

  if (!seed.valid) {
    std::cerr << "ERROR: bkg-only seed fit failed for " << bkg_input
              << std::endl;
    return out;
  }

  RooFitUtils sim;
  if (interactive)
    sim.SetInteractive();

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back(ge_mu_init);

  sim.AddChannel("bkg", bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg_peaks, bkg_peak_mus, bkg_flat, use_step, use_low_exp,
                 use_low_lin, use_high_exp);
  sim.AddChannel("sig", sig_events, sig_lo, sig_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg_peaks + 1, sig_mus, sig_flat, use_step, use_low_exp,
                 use_low_lin, use_high_exp);

  for (Int_t i = 0; i < n_bkg_peaks; i++) {
    sim.LinkPeakShape("sig", i, "bkg", i);
  }
  sim.SeedChannel("bkg", seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_input, "BkgGe_sim");

  if (results.size() < 2)
    return out;

  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;
  return out;
}

void Fits() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/raw",
                                project_root + "/root_files");

  Bool_t interactive = kTRUE;

  std::vector<TString> run_names;
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> reduced_chi2;

  std::vector<Double_t> pb_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  std::vector<Double_t> one_bkg_mus = {(Double_t)E_BKG_73KEV};

  BkgGeSimResult cd_10 =
      FitBkgGeSimultaneous(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                           Constants::CDSHIELDSIGNAL_10PERCENT_20260113, 65, 81,
                           64, 80, pb_mus, E_GE_73M, kTRUE, kTRUE, interactive);
  BkgGeSimResult cd_25 =
      FitBkgGeSimultaneous(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                           Constants::CDSHIELDSIGNAL_25PERCENT_20260113, 66, 81,
                           65, 81, pb_mus, E_GE_73M, kTRUE, kTRUE, interactive);
  BkgGeSimResult cu_0113 = FitBkgGeSimultaneous(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113, 65, 82, 62, 80, pb_mus,
      E_GE_73M, kFALSE, kTRUE, interactive);
  BkgGeSimResult cu_0114 =
      FitBkgGeSimultaneous(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                           Constants::CUSHIELDSIGNAL_10PERCENT_20260114, 66, 82,
                           63, 80, pb_mus, E_GE_73M, kTRUE, kTRUE, interactive);

  BkgGeSimResult noshield_15 = FitBkgGeSimultaneous(
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115, 67, 77, 64, 77, one_bkg_mus,
      E_GE_73M, kTRUE, kTRUE, interactive);
  BkgGeSimResult graphite_16 = FitBkgGeSimultaneous(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116, 67, 80, 60,
      77, one_bkg_mus, E_GE_73M, kTRUE, kTRUE, interactive);

  BkgGeSimResult pairs[6] = {cd_10,   cd_25,       cu_0113,
                             cu_0114, noshield_15, graphite_16};
  TString labels[6] = {"Cd Shield Signal 10% (01/13)",
                       "Cd Shield Signal 25% (01/13)",
                       "Cu Shield Signal 10% (01/13)",
                       "Cu Shield Signal 10% (01/14)",
                       "No Shield Signal 5% (01/15)",
                       "No Shield Graphite Castle Signal 10% (01/16)"};

  for (Int_t i = 0; i < 6; i++) {
    if (!pairs[i].valid)
      continue;
    const PeakFitResult &ge = pairs[i].sig_channel.peaks.back();
    run_names.push_back(labels[i]);
    mu.push_back(ge.mu);
    mu_errors.push_back(ge.mu_error);
    reduced_chi2.push_back(pairs[i].sig_channel.reduced_chi2);
  }

  std::cout << std::endl;
  std::cout << "Individual Run Results (Ge Peak mu):" << std::endl;
  for (size_t i = 0; i < mu.size(); ++i) {
    std::cout << std::left << std::setw(50) << run_names[i] << ": "
              << std::fixed << std::setprecision(4) << mu[i] << " +/- "
              << mu_errors[i] << " keV"
              << " (chi2/ndf = " << std::setprecision(3) << reduced_chi2[i]
              << ")" << std::endl;
  }

  Float_t sum_weights = 0.0;
  Float_t weighted_sum = 0.0;
  for (size_t i = 0; i < mu.size(); ++i) {
    if (mu_errors[i] > 0) {
      Float_t w = 1.0 / (mu_errors[i] * mu_errors[i]);
      weighted_sum += mu[i] * w;
      sum_weights += w;
    }
  }
  if (sum_weights > 0) {
    Float_t combined_mu = weighted_sum / sum_weights;
    Float_t combined_error = std::sqrt(1.0 / sum_weights);
    std::cout << std::endl;
    std::cout << "Combined Ge mu: " << std::fixed << std::setprecision(4)
              << combined_mu << " +/- " << combined_error << " keV"
              << std::endl;
  }
}
