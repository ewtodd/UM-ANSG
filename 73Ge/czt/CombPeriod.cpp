// The period of the comb in the summed CZT energies.
//
// Every pixel's energies sit on their own fine grid; summed over the pixels
// they show a comb, and the chain bins its plots, chi2 and the 95.902 keV
// line's fits one period wide (COMB_PERIOD_RAW_KEV in GainTransferFit.hpp) so
// that each bin holds one tooth. This measures the period in a run's filtered
// energies (after the pixel gains, before the day map and the calibration):
// the spectrum over 50-110 keV in 1 eV bins, divided by its running mean over
// 2 keV, and the period between 0.20 and 0.35 keV with the most power in the
// remaining ripple. One run of each analysed day, since the constant claims
// the same comb on 01/13 and 01/15. Writes results/comb_period.result.
//
// Usage: root -b -q 'CombPeriod.cpp+'
#include "Constants.hpp"
#include "GainTransferFit.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "ProjectPaths.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

const Double_t SPECTRUM_LO_KEV = 50.0;
const Double_t SPECTRUM_HI_KEV = 110.0;
const Double_t BIN_KEV = 0.001;
// The running mean that takes out the spectrum's shape: several periods wide.
const Int_t SMOOTHING_BINS = 2001;
const Double_t PERIOD_LO_KEV = 0.200;
const Double_t PERIOD_HI_KEV = 0.350;
const Double_t PERIOD_STEP_KEV = 0.00005;

struct CombMeasurement {
  Double_t period_kev = 0;
  // The peak's power over the median power of the scanned periods.
  Double_t contrast = 0;
  Bool_t valid = kFALSE;
};

TH1D *FilteredSpectrum(const TString &input_name) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot open filtered/" << input_name << ".root"
              << std::endl;
    delete file;
    return nullptr;
  }
  const Int_t n_bins =
      TMath::Nint((SPECTRUM_HI_KEV - SPECTRUM_LO_KEV) / BIN_KEV);
  TH1D *spectrum = new TH1D("comb_spectrum_" + input_name, "", n_bins,
                            SPECTRUM_LO_KEV, SPECTRUM_HI_KEV);
  spectrum->SetDirectory(nullptr);
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t energy = 0;
    tree->SetBranchAddress("energykeV", &energy);
    Long64_t n_entries = tree->GetEntries();
    for (Long64_t j = 0; j < n_entries; j++) {
      tree->GetEntry(j);
      spectrum->Fill(energy);
    }
    tree->ResetBranchAddresses();
  }
  file->Close();
  delete file;
  return spectrum;
}

// The power of the ripple (the spectrum over its running mean, minus one) at
// each scanned period, and the period with the most, refined by a parabola
// through its neighbours.
CombMeasurement MeasureComb(const TH1D *spectrum) {
  CombMeasurement out;
  const Int_t n = spectrum->GetNbinsX();
  const Int_t half = SMOOTHING_BINS / 2;
  std::vector<Double_t> cumulative(n + 1, 0.0);
  for (Int_t b = 1; b <= n; b++)
    cumulative[b] = cumulative[b - 1] + spectrum->GetBinContent(b);
  std::vector<Double_t> energies, ripple;
  for (Int_t b = 1 + half; b <= n - half; b++) {
    Double_t mean =
        (cumulative[b + half] - cumulative[b - half - 1]) / SMOOTHING_BINS;
    if (mean <= 0)
      continue;
    energies.push_back(spectrum->GetBinCenter(b));
    ripple.push_back(spectrum->GetBinContent(b) / mean - 1.0);
  }
  const Int_t n_periods =
      TMath::Nint((PERIOD_HI_KEV - PERIOD_LO_KEV) / PERIOD_STEP_KEV) + 1;
  std::vector<Double_t> periods, powers;
  for (Int_t p = 0; p < n_periods; p++) {
    Double_t period = PERIOD_LO_KEV + p * PERIOD_STEP_KEV;
    Double_t sum_cos = 0, sum_sin = 0;
    for (size_t k = 0; k < energies.size(); k++) {
      Double_t phase = 2.0 * TMath::Pi() * energies[k] / period;
      sum_cos += ripple[k] * TMath::Cos(phase);
      sum_sin += ripple[k] * TMath::Sin(phase);
    }
    periods.push_back(period);
    powers.push_back(sum_cos * sum_cos + sum_sin * sum_sin);
  }
  Int_t best =
      (Int_t)(std::max_element(powers.begin(), powers.end()) - powers.begin());
  if (best == 0 || best == n_periods - 1) {
    std::cerr << "ERROR: the comb's period lies at the edge of the scan"
              << std::endl;
    return out;
  }
  Double_t below = powers[best - 1], at = powers[best],
           above = powers[best + 1];
  Double_t offset = 0.5 * (below - above) / (below - 2.0 * at + above);
  out.period_kev = periods[best] + offset * PERIOD_STEP_KEV;
  std::vector<Double_t> sorted = powers;
  std::sort(sorted.begin(), sorted.end());
  out.contrast = at / sorted[sorted.size() / 2];
  out.valid = kTRUE;
  return out;
}

} // namespace

void CombPeriod() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/czt/chain",
                                project_root + "/root_files");
  const Int_t n_runs = 2;
  const TString runs[n_runs] = {Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                                Constants::NOSHIELDSIGNAL_5PERCENT_20260115};
  const char *keys[n_runs] = {"cu10_20260113", "noshield_20260115"};

  const TString result_path = project_root + "/results/comb_period.result";
  std::ofstream result(result_path.Data());
  result << "# Written by czt/CombPeriod.cpp: the comb's period in each run's "
            "filtered energies, keV."
         << std::endl;
  result << std::fixed << std::setprecision(4);
  for (Int_t r = 0; r < n_runs; r++) {
    TH1D *spectrum = FilteredSpectrum(runs[r]);
    if (!spectrum)
      continue;
    CombMeasurement comb = MeasureComb(spectrum);
    delete spectrum;
    if (!comb.valid)
      continue;
    std::cout << "COMB-PERIOD " << runs[r] << ": " << std::fixed
              << std::setprecision(4) << comb.period_kev << " keV (power "
              << std::setprecision(0) << comb.contrast
              << " times the scan's median); the chain uses "
              << std::setprecision(4) << COMB_PERIOD_RAW_KEV << std::endl;
    result << keys[r] << " " << comb.period_kev << std::endl;
  }
  std::cout << "Wrote " << result_path << std::endl;
}
