#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TColor.h>
#include <TROOT.h>
#include <TSystem.h>
#include <vector>

const Int_t FINE_BIN_WIDTH = 75;
const Float_t LO_WF_MIN = 0;
const Float_t LO_WF_MAX = 1750;
const Int_t N_FINE_BINS =
    static_cast<Int_t>((LO_WF_MAX - LO_WF_MIN) / FINE_BIN_WIDTH);
const Int_t TARGET_EVENTS = 5000;
const Int_t MIN_EVENTS = 2500;
const Int_t DISPLAY_SAMPLES = 100;

Int_t AlphaColor(Float_t t) {
  Int_t r = 255 - static_cast<Int_t>(75 * t);
  Int_t g = 180 - static_cast<Int_t>(180 * t);
  Int_t b = 120 - static_cast<Int_t>(120 * t);
  return TColor::GetColor(r, g, b);
}

Int_t GammaColor(Float_t t) {
  Int_t r = 120 - static_cast<Int_t>(120 * t);
  Int_t g = 180 - static_cast<Int_t>(180 * t);
  Int_t b = 255 - static_cast<Int_t>(75 * t);
  return TColor::GetColor(r, g, b);
}

void ComputeAndPlot(const TString output_name, const TString label,
                    TGraph *alpha_template, TGraph *gamma_template,
                    TGraph *compare_template, Bool_t compare_is_alpha) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString filepath = project_root + "/root_files/" + output_name + ".root";
  TFile *file = TFile::Open(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << filepath << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree" << std::endl;
    file->Close();
    delete file;
    return;
  }

  TArrayF *samples = nullptr;
  Float_t light_output;
  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("light_output", &light_output);

  tree->GetEntry(0);
  Int_t wavelength = samples->GetSize();

  std::vector<std::vector<Double_t>> fine_sums(
      N_FINE_BINS, std::vector<Double_t>(wavelength, 0.0));
  std::vector<Int_t> fine_counts(N_FINE_BINS, 0);

  Int_t n_entries = tree->GetEntries();

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    if (light_output < LO_WF_MIN || light_output >= LO_WF_MAX)
      continue;

    Int_t bin = static_cast<Int_t>((light_output - LO_WF_MIN) / FINE_BIN_WIDTH);
    if (bin < 0 || bin >= N_FINE_BINS)
      continue;

    for (Int_t j = 0; j < wavelength; j++) {
      fine_sums[bin][j] += samples->At(j);
    }
    fine_counts[bin]++;
  }

  file->Close();
  delete file;

  std::vector<TGraph *> graphs;
  std::vector<TString> bin_labels;
  std::vector<Int_t> colors;
  std::vector<Bool_t> is_alpha;

  std::vector<Double_t> merged_sum(wavelength, 0.0);
  Int_t merged_count = 0;
  Float_t merged_lo_min = LO_WF_MIN;

  Int_t first_adaptive_bin = 200 / FINE_BIN_WIDTH;

  for (Int_t b = 0; b < N_FINE_BINS; b++) {
    for (Int_t j = 0; j < wavelength; j++) {
      merged_sum[j] += fine_sums[b][j];
    }
    merged_count += fine_counts[b];

    Float_t bin_upper = LO_WF_MIN + (b + 1) * FINE_BIN_WIDTH;
    Bool_t is_last = (b == N_FINE_BINS - 1);
    Bool_t below_floor = (b < first_adaptive_bin - 1);

    if (!below_floor && (merged_count >= TARGET_EVENTS || is_last) &&
        merged_count >= MIN_EVENTS) {
      Double_t peak = 0;
      for (Int_t j = 0; j < wavelength; j++) {
        merged_sum[j] /= merged_count;
        if (merged_sum[j] > peak)
          peak = merged_sum[j];
      }

      if (peak > 0) {
        Int_t n_display = std::min(DISPLAY_SAMPLES, wavelength);
        TGraph *g = new TGraph(n_display);
        for (Int_t j = 0; j < n_display; j++) {
          g->SetPoint(j, j, merged_sum[j] / peak);
        }

        Float_t energy_frac = (merged_lo_min + bin_upper) / (2.0 * LO_WF_MAX);
        if (energy_frac > 1.0)
          energy_frac = 1.0;

        Bool_t is_alpha_like = kFALSE;
        if (compare_template) {
          Double_t *cy = compare_template->GetY();
          Int_t cn = compare_template->GetN();
          Int_t last = n_display - 1;
          if (last >= cn)
            last = cn - 1;
          Double_t bin_last = merged_sum[last] / peak;
          Double_t comp_last = cy[last];
          Double_t tolerance = 0.0001;
          if (compare_is_alpha)
            is_alpha_like = !(bin_last > comp_last + tolerance);
          else
            is_alpha_like = (bin_last < comp_last - tolerance);
        }

        graphs.push_back(g);
        bin_labels.push_back(Form("%d-%d", static_cast<Int_t>(merged_lo_min),
                                  static_cast<Int_t>(bin_upper)));
        is_alpha.push_back(is_alpha_like);
        colors.push_back(is_alpha_like ? AlphaColor(energy_frac)
                                       : GammaColor(energy_frac));

        std::cout << output_name << " " << bin_labels.back()
                  << " keVee: " << merged_count << " events" << std::endl;
      }

      merged_count = 0;
      merged_lo_min = bin_upper;
      for (Int_t j = 0; j < wavelength; j++) {
        merged_sum[j] = 0.0;
      }
    }
  }

  Int_t n_graphs = static_cast<Int_t>(graphs.size());
  if (n_graphs == 0) {
    std::cout << "No bins with sufficient statistics for " << output_name
              << std::endl;
    return;
  }

  TCanvas *canvas = new TCanvas(PlottingUtils::GetRandomName(), "", 1200, 800);

  TPad *plot_pad = new TPad("plot_pad", "", 0.0, 0.0, 0.85, 1.0);
  plot_pad->Draw();

  TPad *legend_pad = new TPad("legend_pad", "", 0.82, 0.0, 1.0, 1.0);
  legend_pad->Draw();

  plot_pad->cd();

  for (Int_t i = 0; i < n_graphs; i++) {
    PlottingUtils::ConfigureGraph(graphs[i], colors[i],
                                  ";Sample [2 ns];Amplitude [a.u.]");
    if (i == 0) {
      graphs[i]->Draw("AL");
      graphs[i]->GetYaxis()->SetRangeUser(-0.05, 1.15);
    } else {
      graphs[i]->Draw("L SAME");
    }
  }

  if (alpha_template) {
    alpha_template->Set(std::min(DISPLAY_SAMPLES, alpha_template->GetN()));
    alpha_template->SetLineColor(kRed + 2);
    alpha_template->SetLineWidth(PlottingUtils::GetLineWidth() + 1);
    alpha_template->SetLineStyle(2);
    alpha_template->Draw("L SAME");
  }

  if (gamma_template) {
    gamma_template->Set(std::min(DISPLAY_SAMPLES, gamma_template->GetN()));
    gamma_template->SetLineColor(kBlue + 2);
    gamma_template->SetLineWidth(PlottingUtils::GetLineWidth() + 1);
    gamma_template->SetLineStyle(2);
    gamma_template->Draw("L SAME");
  }

  TString final_label =
      label == Constants::AM241_LABEL ? "(a) " + label : "(b) " + label;
  TLatex *plot_label = PlottingUtils::AddText(final_label, 0.82, 0.82);
  plot_label->SetTextSize(35);
  plot_label->Draw();
  legend_pad->cd();

  TLegend *leg = PlottingUtils::AddLegend(0.0, 0.92, 0.05, 0.92);
  for (Int_t i = 0; i < n_graphs; i++) {
    leg->AddEntry(graphs[i], bin_labels[i], "l");
  }
  if (alpha_template)
    leg->AddEntry(alpha_template, "#alpha template", "l");
  if (gamma_template)
    leg->AddEntry(gamma_template, "#gamma template", "l");
  leg->Draw();

  plot_pad->cd();
  plot_pad->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(canvas, "avg_waveforms_" + output_name, "",
                            PlotSaveOptions::kLOG);

  delete canvas;
  for (Int_t i = 0; i < n_graphs; i++) {
    delete graphs[i];
  }
}

TGraph *LoadTemplate(const TString output_name) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString filepath = project_root + "/root_files/" + output_name + ".root";
  TFile *file = TFile::Open(filepath, "READ");
  if (!file || file->IsZombie())
    return nullptr;

  TGraph *g = static_cast<TGraph *>(file->Get("average_waveform"));
  if (g)
    g = static_cast<TGraph *>(g->Clone());

  file->Close();
  delete file;
  return g;
}

void AverageWaveforms() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");

  TGraph *alpha_template = LoadTemplate(Constants::AM241);
  TGraph *gamma_template = LoadTemplate(Constants::NA22);

  ComputeAndPlot(Constants::AM241, Constants::AM241_LABEL, alpha_template,
                 gamma_template, gamma_template, kFALSE);
  ComputeAndPlot(Constants::NA22, Constants::NA22_LABEL, alpha_template,
                 gamma_template, alpha_template, kTRUE);

  delete alpha_template;
  delete gamma_template;
}
