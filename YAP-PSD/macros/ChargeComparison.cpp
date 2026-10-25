#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <vector>

void CalculateChargeComparison(const std::vector<TString> output_names) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();

  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);
    TString filepath = project_root + "/root_files/" + output_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("features"));

    TArrayF *samples = nullptr;
    Int_t trigger_position;

    tree->SetBranchAddress("Samples", &samples);
    tree->SetBranchAddress("trigger_position", &trigger_position);

    Float_t charge_comparison;

    TBranch *cc_branch = tree->GetBranch("charge_comparison");
    if (cc_branch) {
      tree->SetBranchAddress("charge_comparison", &charge_comparison);
      cc_branch->Reset();
    } else {
      cc_branch = tree->Branch("charge_comparison", &charge_comparison,
                               "charge_comparison/F");
    }

    Int_t n_entries = tree->GetEntries();

    Float_t short_calc, long_calc;

    for (Int_t i = 0; i < n_entries; i++) {
      short_calc = 0;
      long_calc = 0;
      tree->GetEntry(i);

      Int_t n_samples = samples->GetSize();
      Int_t start = std::max<Int_t>(
          trigger_position -
              Constants::DEFAULT_PULSE_INTEGRAL_SETTINGS.pre_gate,
          0);
      Int_t short_end =
          std::min(start + Constants::OPTIMAL_SHORT_GATE, n_samples);
      Int_t long_end =
          std::min(start + Constants::OPTIMAL_LONG_GATE, n_samples);

      for (Int_t j = 0; j < n_samples; j++) {
        Float_t value = samples->GetAt(j);
        if (j < short_end)
          short_calc += value;
        long_calc += value;
      }

      if (long_calc > short_calc) {
        charge_comparison = 1 - (short_calc / long_calc);
      } else {
        charge_comparison = -1;
      }

      cc_branch->Fill();
    }

    tree->Write("features", TObject::kOverwrite);
    file->Close();
    delete file;
    std::cout << "Calculated charge comparison for " << output_name
              << std::endl;
  }
}

void PlotChargeComparison(const std::vector<TString> output_names) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();

  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);
    TString filepath = project_root + "/root_files/" + output_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("features"));

    Float_t charge_comparison;
    Float_t light_output;
    tree->SetBranchAddress("charge_comparison", &charge_comparison);
    tree->SetBranchAddress("light_output", &light_output);

    Int_t n_entries = tree->GetEntries();

    tree->LoadBaskets();

    TH2F *charge_comparison_vs_LO =
        new TH2F(PlottingUtils::GetRandomName(), "", Constants::LO_HIST_NBINS,
                 Constants::LO_HIST_XMIN, Constants::LO_HIST_XMAX,
                 Constants::CC_HIST_NBINS, Constants::CC_HIST_XMIN,
                 Constants::CC_HIST_XMAX);

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      charge_comparison_vs_LO->Fill(light_output, charge_comparison);
    }

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureAndDraw2DHistogram(
        charge_comparison_vs_LO, canvas, ";Light Output [keVee]; PSP_{CC}");
    PlottingUtils::SaveFigure(canvas, "cc_vs_lo_" + output_name, "",
                              PlotSaveOptions::kLINEAR);

    charge_comparison_vs_LO->Write("charge_comparison_vs_LO",
                                   TObject::kOverwrite);

    delete canvas;
    file->Close();
    delete file;
    std::cout << "Plotted charge comparison for " << output_name << std::endl;
  }
}

void ChargeComparison() {
  Bool_t recalculate_cc = kTRUE;

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");

  if (recalculate_cc)
    CalculateChargeComparison(Constants::ALL_OUTPUT_NAMES);

  PlotChargeComparison(Constants::ALL_OUTPUT_NAMES);
}
