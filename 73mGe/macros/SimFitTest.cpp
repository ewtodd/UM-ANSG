#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <vector>

const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;
const Float_t E_GE_73M = 68.752;

std::vector<Double_t> LoadEvents(const TString input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    return events;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TString tree_name = Form("crystal%d_filtered_tree", c);
    TTree *tree = static_cast<TTree *>(file->Get(tree_name));
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

void PrintPeak(const TString tag, const TString label,
               const PeakFitResult &peak, Float_t reduced_chi2) {
  std::cout << "  [" << tag << "] " << label << ": mu = " << std::fixed
            << std::setprecision(4) << peak.mu << " +/- " << peak.mu_error
            << ", sigma = " << peak.sigma << " +/- " << peak.sigma_error
            << ", chi2/ndf = " << std::setprecision(3) << reduced_chi2
            << std::endl;
}

void RunCdShieldSimFit(const TString bkg_input, const TString sig_input,
                       Float_t bkg_lo, Float_t bkg_hi, Float_t sig_lo,
                       Float_t sig_hi) {
  std::cout << "\n=== Sim fit: " << bkg_input << " + " << sig_input
            << " ===" << std::endl;

  std::vector<Double_t> bkg_events = LoadEvents(bkg_input);
  std::vector<Double_t> sig_events = LoadEvents(sig_input);
  if (bkg_events.empty() || sig_events.empty())
    return;

  Bool_t use_flat_bkg = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp = kTRUE;
  Bool_t use_low_lin = kTRUE;
  Bool_t use_high_exp = kTRUE;

  RooFitUtils bkg_fitter(bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV,
                         use_flat_bkg, use_step, use_low_exp, use_low_lin,
                         use_high_exp);
  bkg_fitter.SetInteractive();
  FitResult pb_seed =
      bkg_fitter.FitDoublePeak(bkg_input, "Pb_KAlpha_seed", E_PB_KA2, E_PB_KA1);

  if (!pb_seed.valid) {
    std::cerr << "ERROR: bkg-only seed fit failed; aborting sim fit"
              << std::endl;
    return;
  }

  std::cout << "\nSeed acquired. Pb K-alpha shape: mu1=" << pb_seed.peaks[0].mu
            << " sigma1=" << pb_seed.peaks[0].sigma
            << " mu2=" << pb_seed.peaks[1].mu
            << " sigma2=" << pb_seed.peaks[1].sigma << std::endl;

  RooFitUtils sim;
  sim.SetInteractive();
  std::vector<Double_t> bkg_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  std::vector<Double_t> sig_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1,
                                   (Double_t)E_GE_73M};
  sim.AddChannel("bkg", bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV, 2,
                 bkg_mus, use_flat_bkg, use_step, use_low_exp, use_low_lin,
                 use_high_exp);
  sim.AddChannel("sig", sig_events, sig_lo, sig_hi, Constants::BIN_WIDTH_KEV, 3,
                 sig_mus, use_flat_bkg, use_step, use_low_exp, use_low_lin,
                 use_high_exp);

  sim.LinkPeakShape("sig", 0, "bkg", 0);
  sim.LinkPeakShape("sig", 1, "bkg", 1);

  sim.SeedChannel("bkg", pb_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(bkg_input, "Ge_sim");

  std::cout << "\nResults:" << std::endl;
  if (results.size() >= 1 && results[0].valid) {
    PrintPeak("sim/bkg", "Pb_KA1", results[0].peaks[0],
              results[0].reduced_chi2);
    PrintPeak("sim/bkg", "Pb_KA2", results[0].peaks[1],
              results[0].reduced_chi2);
  }
  if (results.size() >= 2 && results[1].valid) {
    PrintPeak("sim/sig", "Pb_KA1", results[1].peaks[0],
              results[1].reduced_chi2);
    PrintPeak("sim/sig", "Pb_KA2", results[1].peaks[1],
              results[1].reduced_chi2);
    PrintPeak("sim/sig", "Ge_73m", results[1].peaks[2],
              results[1].reduced_chi2);
  }
}

void SimFitTest() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/sim_fit_test",
                                project_root + "/root_files");

  RunCdShieldSimFit(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                    Constants::CDSHIELDSIGNAL_25PERCENT_20260113, 66, 81, 65,
                    81);
}
