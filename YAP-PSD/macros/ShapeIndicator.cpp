#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <vector>

struct SpectralCuts {
  Float_t min_light_output = 900.0;
  Float_t max_light_output = 1200.0;
  Int_t max_events_per_source = 5000;
};

const SpectralCuts CUTS;

void CalculateAverageWaveforms(const std::vector<TString> output_names,
                               const SpectralCuts &cuts = CUTS) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();
  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);

    TString filename = project_root + "/root_files/" + output_name + ".root";

    TFile *wf_file = TFile::Open(filename, "UPDATE");
    if (!wf_file || wf_file->IsZombie()) {
      std::cout << "Error: Could not load " << output_name << std::endl;
      return;
    }

    TTree *tree = static_cast<TTree *>(wf_file->Get("features"));
    if (!tree) {
      std::cout << "Error: Could not find features tree!" << std::endl;
      wf_file->Close();
      return;
    }

    TArrayF *samples = nullptr;
    Float_t light_output;

    tree->SetBranchAddress("Samples", &samples);
    tree->SetBranchAddress("light_output", &light_output);

    std::vector<std::vector<Float_t>> waveforms;
    Int_t event_counts = 0;

    Float_t min_light_output = output_name == Constants::AM241 ? 1300 : 500;
    Float_t max_light_output = 1550;
    Int_t max_events_per_source = CUTS.max_events_per_source;

    Long64_t n_entries = tree->GetEntries();

    tree->GetEntry(0);
    Int_t wavelength = samples->GetSize();

    for (Long64_t i = 0; i < n_entries; ++i) {
      if (tree->GetEntry(i) <= 0)
        continue;

      if (event_counts >= max_events_per_source)
        continue;

      if (light_output < min_light_output || light_output > max_light_output)
        continue;

      std::vector<Float_t> waveform;

      for (Int_t j = 0; j < wavelength; j++) {
        waveform.push_back(samples->At(j));
      }

      waveforms.push_back(waveform);
      event_counts++;
    }

    std::vector<Float_t> average_waveform(wavelength, 0.0);

    Int_t n_waveforms = waveforms.size();
    for (Int_t wave = 0; wave < n_waveforms; wave++) {
      std::vector<Float_t> wf = waveforms.at(wave);
      for (Int_t size = 0; size < wavelength; size++) {
        average_waveform[size] += wf[size];
      }
    }

    for (Int_t i = 0; i < wavelength; i++) {
      average_waveform[i] /= Float_t(waveforms.size());
    }

    Float_t peak =
        *std::max_element(average_waveform.begin(), average_waveform.end());
    if (peak > 0) {
      for (Int_t i = 0; i < wavelength; i++) {
        average_waveform[i] /= peak;
      }
    }

    TGraph *g = new TGraph(wavelength);
    for (Int_t i = 0; i < wavelength; i++) {
      g->SetPoint(i, i, average_waveform[i]);
    }
    g->SetName("average_waveform");
    PlottingUtils::ConfigureGraph(g,
                                  PlottingUtils::GetDefaultColors().at(entry),
                                  ";Sample [2 ns];Amplitude [a.u.]");

    wf_file->cd();
    g->Write("average_waveform", TObject::kOverwrite);
    wf_file->Close();
  }
}

void CalculateRawWeightingFunction(const TString alpha_output_name,
                                   const TString gamma_output_name) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString alpha_filename =
      project_root + "/root_files/" + alpha_output_name + ".root";
  TFile *alpha_file = TFile::Open(alpha_filename, "UPDATE");
  TGraph *alpha_average =
      static_cast<TGraph *>(alpha_file->Get("average_waveform"));

  TString gamma_filename =
      project_root + "/root_files/" + gamma_output_name + ".root";
  TFile *gamma_file = TFile::Open(gamma_filename, "UPDATE");
  TGraph *gamma_average =
      static_cast<TGraph *>(gamma_file->Get("average_waveform"));

  Int_t wavelength = alpha_average->GetN();
  Double_t *alpha_y = alpha_average->GetY();
  Double_t *gamma_y = gamma_average->GetY();

  TGraph *wf = new TGraph(wavelength);
  Double_t *wf_x = wf->GetX();
  Double_t *wf_y = wf->GetY();

  for (Int_t i = 0; i < wavelength; i++) {
    wf_x[i] = i;
    wf_y[i] = (alpha_y[i] - gamma_y[i]) / (alpha_y[i] + gamma_y[i]);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawGraph(wf, kAzure,
                                       ";Sample [2 ns];Amplitude [a.u.]");
  PlottingUtils::SaveFigure(canvas, "raw_weighting_function", "",
                            PlotSaveOptions::kLINEAR);

  alpha_file->cd();
  wf->Write("raw_weighting_function", TObject::kOverwrite);
  gamma_file->cd();
  wf->Write("raw_weighting_function", TObject::kOverwrite);
  alpha_file->Close();
  gamma_file->Close();
}

void CalculateCleanWeightingFunction(const TString alpha_output_name,
                                     const TString gamma_output_name) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString alpha_filename =
      project_root + "/root_files/" + alpha_output_name + ".root";
  TFile *alpha_file = TFile::Open(alpha_filename, "UPDATE");
  TGraph *alpha_average =
      static_cast<TGraph *>(alpha_file->Get("average_waveform"));

  TString gamma_filename =
      project_root + "/root_files/" + gamma_output_name + ".root";
  TFile *gamma_file = TFile::Open(gamma_filename, "UPDATE");
  TGraph *gamma_average =
      static_cast<TGraph *>(gamma_file->Get("average_waveform"));

  Int_t wavelength = alpha_average->GetN();
  Double_t *alpha_y = alpha_average->GetY();
  Double_t *gamma_y = gamma_average->GetY();

  TGraph *wf = new TGraph(wavelength);
  Double_t *wf_x = wf->GetX();
  Double_t *wf_y = wf->GetY();

  Float_t numerator = 0, denominator = 0;

  for (Int_t i = 0; i < wavelength; i++) {
    numerator = 0;
    denominator = 0;
    wf_x[i] = i;

    if (i < 20) {
      wf_y[i] = 0;
    } else {
      numerator = alpha_y[i] - gamma_y[i];
      denominator = alpha_y[i] + gamma_y[i];
      if (i < 50)
        wf_y[i] = denominator > 5e-2 ? numerator / denominator : 0;
      else
        wf_y[i] = numerator / denominator;
    }
  }

  Double_t peak = 0;
  for (Int_t i = 0; i < wavelength; i++) {
    if (std::abs(wf_y[i]) > std::abs(peak))
      peak = wf_y[i];
  }

  for (Int_t i = 0; i < wavelength; i++) {
    wf_y[i] /= std::abs(peak);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawGraph(alpha_average, kRed + 2,
                                       ";Sample [2 ns];Amplitude [a.u.]");
  PlottingUtils::ConfigureGraph(gamma_average, kBlue + 2,
                                ";Sample [2 ns];Amplitude [a.u.]");
  gamma_average->Draw("SAME");
  PlottingUtils::ConfigureGraph(wf, kGreen + 1,
                                ";Sample [2 ns];Amplitude [a.u.]");
  wf->Draw("SAME");
  alpha_average->GetYaxis()->SetRangeUser(-1.1, 1.1);

  TLegend *leg = PlottingUtils::AddLegend(0.5, 0.85, 0.65, 0.85);
  leg->AddEntry(alpha_average, "f_{#alpha}(t) (Am-241)", "l");
  leg->AddEntry(gamma_average, "f_{#gamma}(t) (Na-22)", "l");
  leg->AddEntry(wf, "Weighting function p(t)", "l");

  PlottingUtils::SaveFigure(canvas, "clean_weighting_function", "",
                            PlotSaveOptions::kLINEAR);

  alpha_file->cd();
  wf->Write("clean_weighting_function", TObject::kOverwrite);
  gamma_file->cd();
  wf->Write("clean_weighting_function", TObject::kOverwrite);
  alpha_file->Close();
  gamma_file->Close();
}

void CalculateShapeIndicator(const std::vector<TString> output_names) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();
  TGraph *raw_weighting_function = nullptr;
  TGraph *clean_weighting_function = nullptr;

  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);
    TString filepath = project_root + "/root_files/" + output_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("features"));

    if (output_name == Constants::AM241) {
      raw_weighting_function =
          static_cast<TGraph *>(file->Get("raw_weighting_function"));
      clean_weighting_function =
          static_cast<TGraph *>(file->Get("clean_weighting_function"));
    }

    TArrayF *samples = nullptr;
    Int_t trigger_position;

    tree->SetBranchAddress("Samples", &samples);
    tree->SetBranchAddress("trigger_position", &trigger_position);

    Float_t raw_shape_indicator;
    Float_t clean_shape_indicator;

    TBranch *raw_branch = tree->GetBranch("raw_shape_indicator");
    TBranch *clean_branch = tree->GetBranch("clean_shape_indicator");

    if (raw_branch) {
      tree->SetBranchAddress("raw_shape_indicator", &raw_shape_indicator);
      raw_branch->Reset();
    } else {
      raw_branch = tree->Branch("raw_shape_indicator", &raw_shape_indicator,
                                "raw_shape_indicator/F");
    }

    if (clean_branch) {
      tree->SetBranchAddress("clean_shape_indicator", &clean_shape_indicator);
      clean_branch->Reset();
    } else {
      clean_branch =
          tree->Branch("clean_shape_indicator", &clean_shape_indicator,
                       "clean_shape_indicator/F");
    }

    Int_t n_entries = tree->GetEntries();

    Float_t raw_num_calc, raw_denom_calc;
    Float_t clean_num_calc, clean_denom_calc;

    for (Int_t i = 0; i < n_entries; i++) {
      raw_num_calc = 0;
      raw_denom_calc = 0;
      clean_num_calc = 0;
      clean_denom_calc = 0;

      tree->GetEntry(i);

      Int_t n_samples = samples->GetSize();
      Int_t start = std::max<Int_t>(
          trigger_position - Constants::DEFAULT_PROCESSING_CONFIG.pre_gate, 0);
      Int_t long_end =
          std::min(start + Constants::OPTIMAL_LONG_GATE, n_samples);

      Double_t *raw_wf_y = raw_weighting_function->GetY();
      Double_t *clean_wf_y = clean_weighting_function->GetY();
      for (Int_t j = 0; j < n_samples; j++) {
        Float_t value = samples->GetAt(j);
        raw_num_calc += raw_wf_y[j] * value;
        raw_denom_calc += value;
        clean_num_calc += clean_wf_y[j] * value;
        clean_denom_calc += value;
      }

      if (raw_denom_calc > 1e-3) {
        raw_shape_indicator = raw_num_calc / raw_denom_calc;
      } else {
        raw_shape_indicator = -1;
      }

      if (clean_denom_calc > 1e-3) {
        clean_shape_indicator = clean_num_calc / clean_denom_calc;
      } else {
        clean_shape_indicator = -1;
      }

      raw_branch->Fill();
      clean_branch->Fill();
    }

    tree->Write("features", TObject::kOverwrite);
    file->Close();
    delete file;
    std::cout << "Calculated shape indicator for " << output_name << std::endl;
  }
}

void PlotShapeIndicator(const std::vector<TString> output_names) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  Int_t n_files = output_names.size();

  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);
    TString filepath = project_root + "/root_files/" + output_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("features"));

    Float_t raw_shape_indicator, clean_shape_indicator;
    Float_t light_output;
    tree->SetBranchAddress("raw_shape_indicator", &raw_shape_indicator);
    tree->SetBranchAddress("clean_shape_indicator", &clean_shape_indicator);
    tree->SetBranchAddress("light_output", &light_output);

    Int_t n_entries = tree->GetEntries();

    tree->LoadBaskets();

    TH2F *clean_shape_indicator_vs_LO =
        new TH2F(PlottingUtils::GetRandomName(), "", Constants::LO_HIST_NBINS,
                 Constants::LO_HIST_XMIN, Constants::LO_HIST_XMAX,
                 Constants::SI_HIST_NBINS, Constants::SI_HIST_XMIN,
                 Constants::SI_HIST_XMAX);

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      clean_shape_indicator_vs_LO->Fill(light_output, clean_shape_indicator);
    }

    TCanvas *clean_canvas = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureAndDraw2DHistogram(
        clean_shape_indicator_vs_LO, clean_canvas,
        ";Light Output [keVee]; PSP_{SI}");
    PlottingUtils::SaveFigure(clean_canvas, "clean_si_vs_lo_" + output_name, "",
                              PlotSaveOptions::kLINEAR);

    clean_shape_indicator_vs_LO->Write("clean_shape_indicator_vs_LO",
                                       TObject::kOverwrite);

    delete clean_canvas;
    file->Close();
    delete file;
    std::cout << "Plotted shape indicator for " << output_name << std::endl;
  }
}

std::vector<Float_t> LoadSIValues(const TString &output_name,
                                  const SpectralCuts &cuts) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString filepath = project_root + "/root_files/" + output_name + ".root";
  TFile *file = new TFile(filepath, "READ");
  TTree *tree = static_cast<TTree *>(file->Get("features"));

  Float_t clean_shape_indicator, light_output;
  tree->SetBranchAddress("clean_shape_indicator", &clean_shape_indicator);
  tree->SetBranchAddress("light_output", &light_output);
  tree->LoadBaskets();

  std::vector<Float_t> values;
  Int_t n_entries = tree->GetEntries();
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    if (light_output >= cuts.min_light_output &&
        light_output <= cuts.max_light_output) {
      values.push_back(clean_shape_indicator);
    }
  }

  file->Close();
  delete file;
  return values;
}

TF1 *FitSIGaussian(TH1F *hist) {
  if (!hist || hist->GetEntries() == 0)
    return nullptr;

  Int_t first_bin = hist->FindFirstBinAbove(hist->GetMaximum() * 0.01);
  Int_t last_bin = hist->FindLastBinAbove(hist->GetMaximum() * 0.01);
  Double_t fit_min = hist->GetBinCenter(first_bin);
  Double_t fit_max = hist->GetBinCenter(last_bin);

  TF1 *fit = new TF1(PlottingUtils::GetRandomName(), "gaus", fit_min, fit_max);
  fit->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
  fit->SetParLimits(0, hist->GetMaximum() * 0.1, hist->GetMaximum() * 2.0);
  fit->SetParLimits(1, Constants::SI_HIST_XMIN, Constants::SI_HIST_XMAX);
  fit->SetParLimits(2, 0, 1);

  Int_t status = hist->Fit(fit, "RQSN");
  if (status != 0) {
    delete fit;
    return nullptr;
  }
  return fit;
}

void PlotBestSIFOM(const SpectralCuts &cuts = CUTS) {
  TString alpha_output = Constants::AM241;
  TString gamma_output = Constants::NA22;
  TString alpha_label = Constants::AM241_LABEL;
  TString gamma_label = Constants::NA22_LABEL;

  std::vector<Float_t> alpha_si = LoadSIValues(alpha_output, cuts);
  std::vector<Float_t> gamma_si = LoadSIValues(gamma_output, cuts);

  std::cout << "SI FOM plot: " << alpha_si.size() << " alpha events, "
            << gamma_si.size() << " gamma events in " << cuts.min_light_output
            << "-" << cuts.max_light_output << " keVee" << std::endl;

  TH1F *hist_alpha =
      new TH1F(PlottingUtils::GetRandomName(), "", Constants::SI_HIST_NBINS,
               Constants::SI_HIST_XMIN, 0.3);
  TH1F *hist_gamma =
      new TH1F(PlottingUtils::GetRandomName(), "", Constants::SI_HIST_NBINS,
               Constants::SI_HIST_XMIN, 0.3);
  hist_alpha->SetDirectory(0);
  hist_gamma->SetDirectory(0);

  for (size_t i = 0; i < alpha_si.size(); i++)
    hist_alpha->Fill(alpha_si[i]);
  for (size_t i = 0; i < gamma_si.size(); i++)
    hist_gamma->Fill(gamma_si[i]);

  TF1 *fit_alpha = FitSIGaussian(hist_alpha);
  TF1 *fit_gamma = FitSIGaussian(hist_gamma);

  if (!fit_alpha || !fit_gamma) {
    std::cout << "Error: Gaussian fit failed for SI FOM plot" << std::endl;
    delete hist_alpha;
    delete hist_gamma;
    delete fit_alpha;
    delete fit_gamma;
    return;
  }

  Double_t mean_alpha = fit_alpha->GetParameter(1);
  Double_t sigma_alpha = std::abs(fit_alpha->GetParameter(2));
  Double_t mean_gamma = fit_gamma->GetParameter(1);
  Double_t sigma_gamma = std::abs(fit_gamma->GetParameter(2));

  Double_t fwhm_alpha = 2.355 * sigma_alpha;
  Double_t fwhm_gamma = 2.355 * sigma_gamma;
  Double_t separation = std::abs(mean_alpha - mean_gamma);
  Double_t fom = separation / (fwhm_alpha + fwhm_gamma);

  std::cout << "SI FOM = " << fom << std::endl;

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);

  PlottingUtils::ConfigureHistogram(hist_alpha, kRed + 2, "");
  PlottingUtils::ConfigureHistogram(hist_gamma, kBlue + 2, "");

  hist_alpha->SetFillColorAlpha(kRed + 2, 0.3);
  hist_gamma->SetFillColorAlpha(kBlue + 2, 0.3);

  hist_alpha->GetXaxis()->SetTitle("PSP_{SI}");
  hist_alpha->GetYaxis()->SetTitle("Counts");

  Double_t y_max =
      std::max(hist_alpha->GetMaximum(), hist_gamma->GetMaximum()) * 1.2;
  hist_alpha->SetMaximum(y_max);
  hist_alpha->Draw("HIST");
  hist_gamma->Draw("HIST SAME");

  Double_t draw_min_alpha = mean_alpha - 3.0 * sigma_alpha;
  Double_t draw_max_alpha = mean_alpha + 3.0 * sigma_alpha;
  TF1 *draw_fit_alpha = new TF1(PlottingUtils::GetRandomName(), "gaus",
                                draw_min_alpha, draw_max_alpha);
  draw_fit_alpha->SetParameters(fit_alpha->GetParameter(0), mean_alpha,
                                sigma_alpha);
  draw_fit_alpha->SetLineColor(kRed);
  draw_fit_alpha->SetLineWidth(2);
  draw_fit_alpha->SetLineStyle(1);
  draw_fit_alpha->Draw("SAME");

  Double_t draw_min_gamma = mean_gamma - 3.0 * sigma_gamma;
  Double_t draw_max_gamma = mean_gamma + 3.0 * sigma_gamma;
  TF1 *draw_fit_gamma = new TF1(PlottingUtils::GetRandomName(), "gaus",
                                draw_min_gamma, draw_max_gamma);
  draw_fit_gamma->SetParameters(fit_gamma->GetParameter(0), mean_gamma,
                                sigma_gamma);
  draw_fit_gamma->SetLineColor(kBlue);
  draw_fit_gamma->SetLineWidth(2);
  draw_fit_gamma->SetLineStyle(1);
  draw_fit_gamma->Draw("SAME");

  TLegend *leg = PlottingUtils::AddLegend(0.55, 0.75, 0.55, 0.83);
  leg->AddEntry(hist_alpha, Form("%s (#alpha)", alpha_label.Data()), "f");
  leg->AddEntry(hist_gamma, Form("%s (#gamma)", gamma_label.Data()), "f");
  leg->AddEntry(draw_fit_alpha, "#alpha fit", "l");
  leg->AddEntry(draw_fit_gamma, "#gamma fit", "l");

  char fom_text[100];
  sprintf(fom_text, "FOM: %.3f", fom);
  leg->AddEntry((TObject *)0, fom_text, "");
  leg->Draw();

  PlottingUtils::SaveFigure(canvas, "best_si_fom_" + gamma_label, "",
                            PlotSaveOptions::kLOG);

  delete fit_alpha;
  delete fit_gamma;
  delete canvas;
  delete hist_alpha;
  delete hist_gamma;
}

void ShapeIndicator() {
  Bool_t recalculate_average = kTRUE;
  Bool_t recalculate_si = kTRUE;
  Bool_t replot_si = kTRUE;

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");

  if (recalculate_average)
    CalculateAverageWaveforms(Constants::SINGLE_OUTPUT_NAMES);

  CalculateRawWeightingFunction(Constants::AM241, Constants::NA22);
  CalculateCleanWeightingFunction(Constants::AM241, Constants::NA22);

  if (recalculate_si)
    CalculateShapeIndicator(Constants::ALL_OUTPUT_NAMES);

  if (replot_si)
    PlotShapeIndicator(Constants::ALL_OUTPUT_NAMES);

  PlotBestSIFOM();
}
