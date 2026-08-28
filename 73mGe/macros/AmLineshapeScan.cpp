// Standalone: how far does the Am-241 centroid move between lineshape models?
//
// Am is the single anchor common to EVERY calibration in CalibrationLow.cpp
// (in-situ Am+Pb and line-cal Am+Ba alike) and dE_AM241 weights it at 0.1 eV,
// so a bias in this centroid is common-mode: it shifts every method together
// and no cross-check between them is sensitive to it.
//
// Runs non-interactively so no saved .roofits is loaded -- each variant is an
// independent automated fit of the same spectrum.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "RooFitUtils.hpp"
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
struct Variant {
  TString name;
  Bool_t flat_bkg, step, low_exp, low_lin, high_exp;
  Float_t lo, hi;
};
} // namespace

void AmLineshapeScan() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/am_lineshape_scan",
                                project_root + "/root_files");

  const TString run = Constants::POSTREACTOR_AM241_20260113;
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << run << std::endl;
    return;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *t =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!t)
      continue;
    std::vector<Double_t> ch = RooFitUtils::LoadEventsFromTree(t, "energykeV");
    events.insert(events.end(), ch.begin(), ch.end());
  }
  file->Close();
  std::cout << "Am-241 events loaded: " << events.size() << std::endl;

  // Nominal is the CalibrationLow.cpp FitCalPeak configuration:
  // flat bkg, no step, low-exp + low-lin tails, high-exp tail, range 51-71.
  std::vector<Variant> vs = {
      {"nominal (as in FitCalPeak)", kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, 51,
       71},
      {"no high-exp tail", kTRUE, kFALSE, kTRUE, kTRUE, kFALSE, 51, 71},
      {"no low-lin tail", kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, 51, 71},
      {"no low-exp tail", kTRUE, kFALSE, kFALSE, kTRUE, kTRUE, 51, 71},
      {"low-exp only", kTRUE, kFALSE, kTRUE, kFALSE, kFALSE, 51, 71},
      {"step enabled", kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, 51, 71},
      {"linear background", kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, 51, 71},
      {"narrow range 55-65", kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, 55, 65},
      {"wide range 48-71", kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, 48, 71},
  };

  std::cout << std::endl
            << std::left << std::setw(30) << "variant" << std::right
            << std::setw(12) << "mu [keV]" << std::setw(10) << "err[eV]"
            << std::setw(10) << "sigma" << std::setw(10) << "chi2/ndf"
            << std::endl;
  std::cout << std::string(72, '-') << std::endl;

  for (size_t i = 0; i < vs.size(); i++) {
    const Variant &v = vs[i];
    RooFitUtils fitter(events, v.lo, v.hi, Constants::BIN_WIDTH_KEV, v.flat_bkg,
                       v.step, v.low_exp, v.low_lin, v.high_exp);
    FitResult r = fitter.FitSinglePeak(run, Form("AmScan%zu", i));
    if (!r.valid || r.peaks.empty()) {
      std::cout << std::left << std::setw(30) << v.name << "   FIT FAILED"
                << std::endl;
      continue;
    }
    std::cout << std::left << std::setw(30) << v.name << std::right
              << std::fixed << std::setprecision(4) << std::setw(12)
              << r.peaks.at(0).mu << std::setprecision(1) << std::setw(10)
              << r.peaks.at(0).mu_error * 1000.0 << std::setprecision(4)
              << std::setw(10) << r.peaks.at(0).sigma << std::setprecision(2)
              << std::setw(10) << r.reduced_chi2 << std::endl;
  }
  std::cout << std::endl
            << "Saved interactive value in use by CalibrationLow.cpp: "
               "59.7075 +/- 0.0004 keV"
            << std::endl;
  std::cout << "Sensitivity: ~57 eV lower Am  ->  ~17 eV higher Ge energy."
            << std::endl;
}
