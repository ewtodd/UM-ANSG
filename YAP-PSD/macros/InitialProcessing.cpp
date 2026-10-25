#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>

void InitialPlots(
    const std::vector<TString> output_names,
    const std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString output_name = output_names.at(i);
    TString output_filepath =
        project_root + "/root_files/" + output_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t long_integral;
    features_tree->SetBranchAddress("long_integral", &long_integral);

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);

    TH1F *long_integral_hist =
        new TH1F("",
                 Form("; Pulse Integral [a.u.]; Counts / %d a.u.",
                      Constants::PI_BIN_WIDTH),
                 Constants::PI_HIST_NBINS, Constants::PI_HIST_XMIN,
                 Constants::PI_HIST_XMAX);

    Int_t num_entries = features_tree->GetEntries();

    for (Int_t i = 0; i < num_entries; i++) {
      features_tree->GetEntry(i);
      long_integral_hist->Fill(long_integral);
    }

    Int_t color = colors.at(i);
    PlottingUtils::ConfigureAndDrawHistogram(long_integral_hist, color);
    PlottingUtils::SaveFigure(canvas, output_name + "_long_integral", "",
                              PlotSaveOptions::kLOG);

    output->cd();
    long_integral_hist->Write("long_integral", TObject::kOverwrite);
    output->Close();

    delete canvas;
  }
}

void InitialProcessing() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");

  Int_t max_threads = 16;
  Bool_t reprocess = kTRUE;

  if (reprocess) {
    const TString input_dir = project_root + "/root_files/";
    std::vector<TString> input_filepaths;
    Int_t n_inputs = Int_t(Constants::ALL_FILENAMES.size());
    for (Int_t i = 0; i < n_inputs; i++) {
      input_filepaths.push_back(input_dir + Constants::ALL_FILENAMES.at(i));
    }
    PulseIntegralProcessor processor(
        Constants::DEFAULT_PROCESSING_CONFIG,
        Constants::DEFAULT_PULSE_INTEGRAL_SETTINGS);
    processor.ProcessFilesParallel(input_filepaths, Constants::ALL_OUTPUT_NAMES,
                                   max_threads);
  }

  InitialPlots(Constants::ALL_OUTPUT_NAMES);
}
