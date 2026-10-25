#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>

void InitialPlots(
    const std::vector<TString> output_names,
    const std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {

  Int_t n_files = output_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString output_name = output_names.at(i);
    TString output_filepath = "root_files/" + output_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t long_integral;
    features_tree->SetBranchAddress("long_integral", &long_integral);
    Float_t pulse_height;
    features_tree->SetBranchAddress("pulse_height", &pulse_height);

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);

    TH1F *long_integral_hist =
        new TH1F("",
                 Form("; Pulse Integral [a.u.]; Counts / %d a.u.",
                      Constants::PI_BIN_WIDTH),
                 Constants::PI_HIST_NBINS, Constants::PI_HIST_XMIN,
                 Constants::PI_HIST_XMAX);

    TH1F *pulse_height_hist = new TH1F(
        "",
        Form("; Pulse Height [ADC]; Counts / %d ADC", Constants::PH_BIN_WIDTH),
        Constants::PH_HIST_NBINS, Constants::PH_HIST_XMIN,
        Constants::PH_HIST_XMAX);

    Int_t num_entries = features_tree->GetEntries();

    for (Int_t i = 0; i < num_entries; i++) {
      features_tree->GetEntry(i);
      long_integral_hist->Fill(long_integral);
      pulse_height_hist->Fill(pulse_height);
    }

    Int_t color = colors.at(i);

    PlottingUtils::ConfigureAndDrawHistogram(long_integral_hist, color);
    PlottingUtils::SaveFigure(canvas, output_name + "_long_integral");

    PlottingUtils::ConfigureAndDrawHistogram(pulse_height_hist, color);
    canvas->Update();

    PlottingUtils::SaveFigure(canvas, output_name + "_pulse_height");

    output->cd();
    long_integral_hist->Write("long_integral", TObject::kOverwrite);
    pulse_height_hist->Write("pulse_height", TObject::kOverwrite);
    output->Close();

    delete canvas;
  }
}

void InitialProcessing() {
  InitUtils::SetROOTPreferences();

  Int_t max_threads = 16;
  Bool_t reprocess = kTRUE;

  if (reprocess) {
    PulseIntegralProcessor processor(
        Constants::DEFAULT_PROCESSING_CONFIG,
        Constants::DEFAULT_PULSE_INTEGRAL_SETTINGS);
    processor.ProcessFilesParallel(Constants::ALL_FILEPATHS,
                                   Constants::ALL_OUTPUT_NAMES, max_threads);
  }

  InitialPlots(Constants::ALL_OUTPUT_NAMES);
}
