#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include <TArrayF.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <vector>

// Light-output window matches LO_LOWER/LO_UPPER in python/noise_study.py.
const Float_t LO_MIN = 750.0;
const Float_t LO_MAX = 1500.0;

// Stored waveforms are already baseline-subtracted and polarity-inverted, so
// the first N samples have mean ~ 0 by construction and RMS reduces to
// sqrt(sum(s^2)/(N-1)) -- the same quantity upstream computes in
// WaveformProcessingUtils::SubtractBaseline.
Float_t ComputeMeanRMS(const TString output_name) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TFile *file = IO::OpenForReading(output_name + ".root");
  TTree *tree = static_cast<TTree *>(file->Get("features"));

  TArrayF *samples = nullptr;
  Float_t light_output;
  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("light_output", &light_output);

  const Int_t N = Constants::DEFAULT_PROCESSING_CONFIG.num_samples_baseline;

  Double_t sum_rms = 0.0;
  Long64_t n_accepted = 0;
  Long64_t n_entries = tree->GetEntries();

  for (Long64_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    if (light_output < LO_MIN || light_output > LO_MAX)
      continue;
    if (samples->GetSize() < N)
      continue;

    Double_t sum_sq = 0.0;
    for (Int_t j = 0; j < N; j++) {
      Float_t v = samples->GetAt(j);
      sum_sq += v * v;
    }
    sum_rms += TMath::Sqrt(sum_sq / (N - 1));
    n_accepted++;
  }

  file->Close();
  delete file;

  Float_t mean_rms = n_accepted > 0 ? sum_rms / n_accepted : -1.0;
  std::cout << std::fixed << std::setprecision(4) << output_name
            << ": mean baseline RMS = " << mean_rms << " ADC (over "
            << n_accepted << " events in LO [" << LO_MIN << ", " << LO_MAX
            << "] keVee, N=" << N << " baseline samples)" << std::endl;
  return mean_rms;
}

void NoiseRMS() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");

  ComputeMeanRMS(Constants::AM241);
  ComputeMeanRMS(Constants::NA22);
}
