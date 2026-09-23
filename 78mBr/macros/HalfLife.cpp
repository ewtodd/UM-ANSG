#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TSystem.h>
#include <TTree.h>
#include <vector>

void CalculatePSPvsLO(std::vector<TString> input_names,
                      Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t short_integral, long_integral, light_output_keVee, psp;
    features_tree->SetBranchAddress("short_integral", &short_integral);
    features_tree->SetBranchAddress("long_integral", &long_integral);
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->Branch("psp", &psp, "psp/F");

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      psp = long_integral > 1e-3
                ? (long_integral - short_integral) / long_integral
                : 0;
      features_tree->GetBranch("psp")->Fill();
    }

    features_tree->Write("", TObject::kOverwrite);
    output->Close();
  }
}

void PlotPSPvsLO(std::vector<TString> input_names) {
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    const Int_t psp_nbins = 100;
    const Float_t psp_min = 0, psp_max = 1;
    const Float_t psp_bin_width = (psp_max - psp_min) / psp_nbins;
    TH2F *PSPvsLO =
        new TH2F("",
                 Form("; Light Output [keVee]; PSP; Counts / %d keV / %.2f",
                      Constants::LO_BIN_WIDTH, psp_bin_width),
                 Constants::LO_HIST_NBINS, Constants::LO_HIST_XMIN,
                 Constants::LO_HIST_XMAX, psp_nbins, psp_min, psp_max);

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();

    TString output_filepath = "root_files/" + input_name + ".root";
    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));
    Float_t light_output_keVee, psp;
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->SetBranchAddress("psp", &psp);

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      PSPvsLO->Fill(light_output_keVee, psp);
    }

    PlottingUtils::ConfigureAndDraw2DHistogram(PSPvsLO, canvas);
    PlottingUtils::SaveFigure(canvas, input_name + "_psp_vs_light_output", "",
                              PlotSaveOptions::kLINEAR);

    PSPvsLO->Write("PSP vs. Light Output", TObject::kOverwrite);
    output->Close();
    delete canvas;
    delete PSPvsLO;
  }
}

const Float_t CANDIDATE_LO_MIN_KEVEE = 125;
const Float_t CANDIDATE_LO_MAX_KEVEE = 175;

// A candidate is a double only if TSpectrum resolves two peaks; their positions
// seed the template fit. Shared with the Monte Carlo so simulated waveforms
// face the same selection as real ones. SearchHighRes overwrites its input, so
// it works on a copy.
Bool_t FindDoublePeaks(const std::vector<Double_t> &waveform,
                       Int_t &peak1_position, Int_t &peak2_position) {
  std::vector<Double_t> work(waveform);
  Int_t old_ignore_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;
  TSpectrum spectrum(10);
  Int_t npeaks = spectrum.SearchHighRes(
      work.data(), work.data(), (Int_t)work.size(), 2, 5, kFALSE, 3, kTRUE, 3);
  gErrorIgnoreLevel = old_ignore_level;
  if (npeaks < 2)
    return kFALSE;
  Double_t *positions = spectrum.GetPositionX();
  peak1_position = (Int_t)positions[0];
  peak2_position = (Int_t)positions[1];
  return kTRUE;
}

void GetCandidateWaveforms(std::vector<TString> input_names,
                           Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t entries = input_names.size();

  TString candidates_filepath = "root_files/candidates.root";

  TFile *output = new TFile(candidates_filepath, "RECREATE");
  TTree *output_tree = new TTree("features", "Candidate waveform features.");
  TArrayF *output_samples = nullptr;
  Float_t output_psp, output_light_output_keVee;

  output_tree->Branch("light_output", &output_light_output_keVee,
                      "light_output/F");
  output_tree->Branch("psp", &output_psp, "output_psp/F");
  output_tree->Branch("Samples", &output_samples);

  for (Int_t i = 0; i < entries; i++) {

    TString input_name = input_names[i];

    TString input_filepath = "root_files/" + input_name + ".root";

    TFile *input = new TFile(input_filepath, "READ");
    TTree *features_tree = static_cast<TTree *>(input->Get("features"));

    Float_t light_output_keVee, psp;
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->SetBranchAddress("psp", &psp);

    TArrayF *samples = nullptr;
    features_tree->SetBranchAddress("Samples", &samples);

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      output_psp = psp;
      output_light_output_keVee = light_output_keVee;
      if (0.05 < psp && psp < 0.275 &&
          output_light_output_keVee < CANDIDATE_LO_MAX_KEVEE &&
          output_light_output_keVee > CANDIDATE_LO_MIN_KEVEE) {
        output_samples = new TArrayF(*samples);
        output_tree->Fill();
      }
    }
    input->Close();
  }

  output->cd();
  output_tree->Write("", TObject::kOverwrite);
  output->Close();
}

void CreateTemplateWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString template_148keV_filepath = "root_files/calibration_Eu152.root";
  TFile *input_148keV = new TFile(template_148keV_filepath, "READ");
  TTree *features_tree_148keV =
      static_cast<TTree *>(input_148keV->Get("features"));

  TArrayF *samples_148keV = nullptr;
  Float_t pulse_height_148keV;
  features_tree_148keV->SetBranchAddress("Samples", &samples_148keV);
  features_tree_148keV->SetBranchAddress("pulse_height", &pulse_height_148keV);

  Int_t num_entries_148keV = features_tree_148keV->GetEntries();
  std::vector<std::vector<Double_t>> selected_waveforms_148keV;
  Int_t nsamples_148keV = 0;

  for (Int_t i = 0; i < num_entries_148keV; i++) {
    features_tree_148keV->GetEntry(i);

    if (pulse_height_148keV > 1500 && pulse_height_148keV < 1800) {
      nsamples_148keV = samples_148keV->GetSize();
      std::vector<Double_t> waveform(nsamples_148keV);
      for (Int_t j = 0; j < nsamples_148keV; j++) {
        waveform[j] = (Double_t)samples_148keV->At(j);
      }
      selected_waveforms_148keV.push_back(waveform);
    }
  }

  input_148keV->Close();

  std::cout << "Found " << selected_waveforms_148keV.size()
            << " waveforms for template 148 keV" << std::endl;

  std::vector<Double_t> template_waveform_148keV(nsamples_148keV, 0.0);
  for (const auto &wf : selected_waveforms_148keV) {
    for (Int_t j = 0; j < nsamples_148keV; j++) {
      template_waveform_148keV[j] += wf[j];
    }
  }

  for (Int_t j = 0; j < nsamples_148keV; j++) {
    template_waveform_148keV[j] /= selected_waveforms_148keV.size();
  }

  std::vector<Double_t> x_values(nsamples_148keV);
  for (Int_t j = 0; j < nsamples_148keV; j++) {
    x_values[j] = 2 * j;
  }

  TString output_filepath = "root_files/template_waveforms.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TGraph *template_graph_148keV = new TGraph(nsamples_148keV, x_values.data(),
                                             template_waveform_148keV.data());
  template_graph_148keV->SetName("template_148keV");
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      template_graph_148keV, kRed,
      "122 keV 152Eu Template Waveform;Time [ns];ADC");
  template_graph_148keV->SetLineWidth(2);
  template_graph_148keV->Draw();
  PlottingUtils::SaveFigure(canvas, "template_148keV", "",
                            PlotSaveOptions::kLINEAR);
  template_graph_148keV->Write("", TObject::kOverwrite);

  TString template_32keV_filepath = "root_files/calibration_Eu152.root";
  TFile *input_32keV = new TFile(template_32keV_filepath, "READ");
  TTree *features_tree_32keV =
      static_cast<TTree *>(input_32keV->Get("features"));

  TArrayF *samples_32keV = nullptr;
  Float_t pulse_height_32keV;
  features_tree_32keV->SetBranchAddress("Samples", &samples_32keV);
  features_tree_32keV->SetBranchAddress("pulse_height", &pulse_height_32keV);

  Int_t num_entries_32keV = features_tree_32keV->GetEntries();
  std::vector<std::vector<Double_t>> selected_waveforms_32keV;
  Int_t nsamples_32keV = 0;

  for (Int_t i = 0; i < num_entries_32keV; i++) {
    features_tree_32keV->GetEntry(i);

    if (pulse_height_32keV > 450 && pulse_height_32keV < 600) {
      nsamples_32keV = samples_32keV->GetSize();
      std::vector<Double_t> waveform(nsamples_32keV);
      for (Int_t j = 0; j < nsamples_32keV; j++) {
        waveform[j] = (Double_t)samples_32keV->At(j);
      }
      selected_waveforms_32keV.push_back(waveform);
    }
  }

  input_32keV->Close();

  std::cout << "Found " << selected_waveforms_32keV.size()
            << " waveforms for template 32 keV" << std::endl;

  std::vector<Double_t> template_waveform_32keV(nsamples_32keV, 0.0);
  for (const auto &wf : selected_waveforms_32keV) {
    for (Int_t j = 0; j < nsamples_32keV; j++) {
      template_waveform_32keV[j] += wf[j];
    }
  }

  for (Int_t j = 0; j < nsamples_32keV; j++) {
    template_waveform_32keV[j] /= selected_waveforms_32keV.size();
  }

  TGraph *template_graph_32keV = new TGraph(nsamples_32keV, x_values.data(),
                                            template_waveform_32keV.data());
  template_graph_32keV->SetName("template_32keV");
  delete canvas;
  canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      template_graph_32keV, kBlue,
      "33.4 keV La K-alpha Template Waveform;Time [ns];ADC");
  template_graph_32keV->SetLineWidth(2);
  template_graph_32keV->Draw();
  PlottingUtils::SaveFigure(canvas, "template_32keV", "",
                            PlotSaveOptions::kLINEAR);
  output->cd();
  template_graph_32keV->Write("", TObject::kOverwrite);

  output->Close();
}

void AnalyzeDoubleWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess) {
    return;
  }

  TString input_filepath = "root_files/candidates.root";

  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayF *samples = nullptr;
  Float_t light_output, psp;
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("light_output", &light_output);
  features_tree->SetBranchAddress("psp", &psp);

  TString output_filepath = "root_files/double_waveforms.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree = new TTree("features", "Double waveform features.");
  TArrayF *output_samples = nullptr;
  Float_t output_psp, output_light_output;
  Int_t peak1_position, peak2_position;

  output_tree->Branch("light_output", &output_light_output, "light_output/F");
  output_tree->Branch("psp", &output_psp, "psp/F");
  output_tree->Branch("peak1_position", &peak1_position, "peak1_position/I");
  output_tree->Branch("peak2_position", &peak2_position, "peak2_position/I");
  output_tree->Branch("Samples", &output_samples);

  Int_t num_entries = features_tree->GetEntries();

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    Int_t nsamples = samples->GetSize();
    std::vector<Double_t> waveform_array(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      waveform_array[j] = (Double_t)samples->At(j);
    }

    if (!FindDoublePeaks(waveform_array, peak1_position, peak2_position))
      continue;

    output_psp = psp;
    output_light_output = light_output;
    if (output_samples)
      delete output_samples;
    output_samples = new TArrayF(*samples);
    output_tree->Fill();
  }

  input->Close();
  output->cd();
  output_tree->Write("", TObject::kOverwrite);
  output->Close();
}

struct FitData {
  const std::vector<Double_t> *waveform;
  const Double_t *template_148_y;
  const Double_t *template_32_y;
  Int_t template_npoints;
  Int_t nsamples;
  Int_t candidate_peak_pos;
  Int_t template_148_peak_pos;
  Int_t template_32_peak_pos;
  Double_t calibration_slope;
  Double_t template_148_peak_height;
  Double_t template_32_peak_height;
};

Double_t CalculateChi2(const Double_t *params, const FitData &data) {
  Double_t energy_148 = params[0];
  Double_t energy_32 = params[1];
  Double_t time_diff = params[2];

  Double_t ph_148 = energy_148 / data.calibration_slope;
  Double_t scale_148 = ph_148 / data.template_148_peak_height;

  Double_t ph_32 = energy_32 / data.calibration_slope;
  Double_t scale_32 = ph_32 / data.template_32_peak_height;

  Double_t time_diff_samples_exact = time_diff / 2.0;

  Double_t chi2 = 0.0;

  for (Int_t j = 0; j < data.nsamples; j++) {
    Double_t template_val = 0.0;

    Double_t exact_index_148 =
        j - data.candidate_peak_pos + data.template_148_peak_pos;

    if (exact_index_148 >= 0 && exact_index_148 < data.template_npoints - 1) {
      Int_t idx_low = (Int_t)exact_index_148;
      Int_t idx_high = idx_low + 1;
      Double_t frac = exact_index_148 - idx_low;

      Double_t interp_val = (1.0 - frac) * data.template_148_y[idx_low] +
                            frac * data.template_148_y[idx_high];
      template_val += interp_val * scale_148;
    } else if (exact_index_148 >= 0 &&
               exact_index_148 < data.template_npoints) {
      template_val += data.template_148_y[(Int_t)exact_index_148] * scale_148;
    }

    Double_t peak_32_pos_exact =
        data.candidate_peak_pos + time_diff_samples_exact;
    Double_t exact_index_32 = j - peak_32_pos_exact + data.template_32_peak_pos;

    if (exact_index_32 >= 0 && exact_index_32 < data.template_npoints - 1) {
      Int_t idx_low = (Int_t)exact_index_32;
      Int_t idx_high = idx_low + 1;
      Double_t frac = exact_index_32 - idx_low;

      Double_t interp_val = (1.0 - frac) * data.template_32_y[idx_low] +
                            frac * data.template_32_y[idx_high];
      template_val += interp_val * scale_32;
    } else if (exact_index_32 >= 0 && exact_index_32 < data.template_npoints) {
      template_val += data.template_32_y[(Int_t)exact_index_32] * scale_32;
    }

    Double_t diff = (*data.waveform)[j] - template_val;
    chi2 += diff * diff;
  }

  return chi2;
}

void FitDoubleWaveforms(const TString &input_name = "double_waveforms.root",
                        const TString &output_name = "fitted_doubles.root",
                        Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString calibration_filepath = "root_files/calibration_function.root";
  TFile *calib_file = new TFile(calibration_filepath, "READ");
  TF1 *calibration_function =
      static_cast<TF1 *>(calib_file->Get("calibration"));

  Double_t calibration_slope = calibration_function->GetParameter(1);

  TString template_filepath = "root_files/template_waveforms.root";
  TFile *template_file = new TFile(template_filepath, "READ");
  TGraph *template_graph_148 =
      static_cast<TGraph *>(template_file->Get("template_148keV"));
  TGraph *template_graph_32 =
      static_cast<TGraph *>(template_file->Get("template_32keV"));

  Double_t *template_148_y = template_graph_148->GetY();
  Double_t *template_32_y = template_graph_32->GetY();
  Int_t template_npoints = template_graph_148->GetN();

  Int_t template_148_peak_pos = TMath::LocMax(template_npoints, template_148_y);
  Int_t template_32_peak_pos = TMath::LocMax(template_npoints, template_32_y);
  Double_t template_148_peak_height = template_148_y[template_148_peak_pos];
  Double_t template_32_peak_height = template_32_y[template_32_peak_pos];

  TString input_filepath = "root_files/" + input_name;
  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayF *samples = nullptr;
  Float_t light_output, psp;
  Int_t peak1_position, peak2_position; // Add these branch addresses
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("light_output", &light_output);
  features_tree->SetBranchAddress("psp", &psp);
  features_tree->SetBranchAddress("peak1_position", &peak1_position);
  features_tree->SetBranchAddress("peak2_position", &peak2_position);

  TString output_filepath = "root_files/" + output_name;
  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree =
      new TTree("features", "Fitted double waveform features.");
  Float_t peak1_light_output, peak2_light_output, time_difference, chi2;

  output_tree->Branch("peak1_light_output", &peak1_light_output,
                      "peak1_light_output/F");
  output_tree->Branch("peak2_light_output", &peak2_light_output,
                      "peak2_light_output/F");
  output_tree->Branch("time_difference", &time_difference, "time_difference/F");
  output_tree->Branch("chi2", &chi2, "chi2/F");

  Int_t feature_tree_entries = features_tree->GetEntries();
  Int_t num_entries =
      TMath::Min(1000000, feature_tree_entries); // Back to processing all
  Int_t plot_count = TMath::Min(100, num_entries);

  std::cout << "Fitting " << num_entries << " candidate waveforms..."
            << std::endl;

  ROOT::Math::Minimizer *minimizer =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  minimizer->SetMaxFunctionCalls(1000);
  minimizer->SetMaxIterations(1000);
  minimizer->SetTolerance(0.01);
  minimizer->SetPrintLevel(-1);

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    if (i % 10000 == 0) { // More frequent updates
      std::cout << "Processing waveform " << i << "/" << num_entries
                << std::endl;
    }

    Int_t nsamples = samples->GetSize();
    std::vector<Double_t> waveform(nsamples);
    std::vector<Double_t> x_values(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      waveform[j] = (Double_t)samples->At(j);
      x_values[j] = j * 2.0;
    }

    Double_t *waveform_ptr = waveform.data();
    Int_t candidate_peak_pos = TMath::LocMax(nsamples, waveform_ptr);
    Double_t candidate_peak_height = waveform[candidate_peak_pos];

    Double_t initial_energy_148 =
        calibration_function->Eval(candidate_peak_height);

    Double_t initial_time_diff =
        TMath::Abs(peak2_position - peak1_position) * 2.0;

    Double_t second_peak_height = 0.0;
    Int_t second_peak_sample = (Int_t)peak2_position;
    if (second_peak_sample >= 0 && second_peak_sample < nsamples) {
      second_peak_height = waveform[second_peak_sample];
    }

    Double_t initial_energy_32 = calibration_function->Eval(second_peak_height);

    FitData fitData;
    fitData.waveform = &waveform;
    fitData.template_148_y = template_148_y;
    fitData.template_32_y = template_32_y;
    fitData.template_npoints = template_npoints;
    fitData.nsamples = nsamples;
    fitData.candidate_peak_pos = candidate_peak_pos;
    fitData.template_148_peak_pos = template_148_peak_pos;
    fitData.template_32_peak_pos = template_32_peak_pos;
    fitData.calibration_slope = calibration_slope;
    fitData.template_148_peak_height = template_148_peak_height;
    fitData.template_32_peak_height = template_32_peak_height;

    ROOT::Math::Functor functor(
        [&fitData](const Double_t *params) {
          return CalculateChi2(params, fitData);
        },
        3);

    minimizer->SetFunction(functor);

    minimizer->SetVariable(0, "energy_148", initial_energy_148, 1);
    minimizer->SetVariableLimits(0, initial_energy_148 * 0.8,
                                 initial_energy_148 * 1.2);

    minimizer->SetVariable(1, "energy_32", initial_energy_32, 0.01);
    minimizer->SetVariableLimits(1, 0, initial_energy_32 * 1.4);

    Double_t time_step = 0.01;
    minimizer->SetVariable(2, "time_diff", initial_time_diff, time_step);
    minimizer->SetVariableLimits(2, initial_time_diff * 0.95,
                                 initial_time_diff * 1.05);

    Bool_t fit_success = minimizer->Minimize();

    Double_t best_energy_148, best_energy_32, best_time_diff, best_chi2;

    if (fit_success) {
      const Double_t *fit_results = minimizer->X();
      best_energy_148 = fit_results[0];
      best_energy_32 = fit_results[1];
      best_time_diff = fit_results[2];
      best_chi2 = minimizer->MinValue();
    }

    peak1_light_output = best_energy_148;
    peak2_light_output = best_energy_32;
    time_difference = best_time_diff;
    chi2 = best_chi2;

    if (fit_success)
      output_tree->Fill();

    if (i < plot_count) {
      Double_t ph_148_best = best_energy_148 / calibration_slope;
      Double_t scale_148_best = ph_148_best / template_148_peak_height;

      Double_t ph_32_best = best_energy_32 / calibration_slope;
      Double_t scale_32_best = ph_32_best / template_32_peak_height;

      Int_t time_diff_samples_best = TMath::Nint(best_time_diff / 2.0);

      std::vector<Double_t> best_template(nsamples, 0.0);
      std::vector<Double_t> residuals(nsamples, 0.0);

      for (Int_t j = 0; j < nsamples; j++) {
        Int_t template_index = j - candidate_peak_pos + template_148_peak_pos;
        if (template_index >= 0 && template_index < template_npoints) {
          best_template[j] += template_148_y[template_index] * scale_148_best;
        }
      }

      Int_t peak_32_pos_best = candidate_peak_pos + time_diff_samples_best;
      for (Int_t j = 0; j < nsamples; j++) {
        Int_t template_index = j - peak_32_pos_best + template_32_peak_pos;
        if (template_index >= 0 && template_index < template_npoints) {
          best_template[j] += template_32_y[template_index] * scale_32_best;
        }
      }

      for (Int_t j = 0; j < nsamples; j++) {
        residuals[j] = waveform[j] - best_template[j];
      }

      TGraph *original = new TGraph(nsamples, x_values.data(), waveform.data());
      TGraph *fitted =
          new TGraph(nsamples, x_values.data(), best_template.data());
      TGraph *residual_graph =
          new TGraph(nsamples, x_values.data(), residuals.data());

      TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
      canvas->Divide(1, 2);

      canvas->cd(1);
      PlottingUtils::ConfigureAndDrawGraph(original, kBlue + 1,
                                           ";Time [ns];ADC");
      fitted->SetLineColor(kRed);
      fitted->SetLineWidth(2);
      fitted->Draw("L SAME");

      Double_t peak1_time = peak1_position * 2.0;
      Double_t peak2_time = peak2_position * 2.0;

      Double_t peak1_height = 0.0;
      Double_t peak2_height = 0.0;
      if (peak1_position >= 0 && peak1_position < nsamples) {
        peak1_height = waveform[peak1_position];
      }
      if (peak2_position >= 0 && peak2_position < nsamples) {
        peak2_height = waveform[peak2_position];
      }

      TMarker *marker_peak1 = new TMarker(peak1_time, peak1_height, 22);
      marker_peak1->SetMarkerColor(kGreen + 2);
      marker_peak1->SetMarkerSize(2);
      marker_peak1->Draw();

      TMarker *marker_peak2 = new TMarker(peak2_time, peak2_height, 22);
      marker_peak2->SetMarkerColor(kGreen + 2);
      marker_peak2->SetMarkerSize(2);
      marker_peak2->Draw();

      TLegend *leg1 = new TLegend(0.75, 0.50, 0.88, 0.88);
      leg1->SetBorderSize(1);
      leg1->SetFillColor(kWhite);
      leg1->SetTextSize(0.04);
      leg1->SetMargin(0.1);
      leg1->SetTextFont(132);
      leg1->AddEntry(original, "Original", "lp");
      leg1->AddEntry(fitted, "Best Fit", "l");
      leg1->AddEntry(marker_peak1, "TSpectrum Peak", "p");
      leg1->AddEntry((TObject *)0, Form("E_{1} = %.1f keVee", best_energy_148),
                     "");
      leg1->AddEntry((TObject *)0, Form("E_{2} = %.1f keVee", best_energy_32),
                     "");
      leg1->AddEntry((TObject *)0,
                     Form("#Deltat_{init} = %.0f ns", initial_time_diff), "");
      leg1->AddEntry((TObject *)0,
                     Form("#Deltat_{fit} = %.0f ns", best_time_diff), "");
      leg1->AddEntry((TObject *)0, Form("#chi^{2} = %.0f", chi2), "");
      leg1->Draw();

      canvas->cd(2);
      PlottingUtils::ConfigureAndDrawGraph(residual_graph, kGreen + 2,
                                           ";Time [ns];Residual");

      TLine *zero_line = new TLine(x_values[0], 0, x_values[nsamples - 1], 0);
      zero_line->SetLineStyle(2);
      zero_line->SetLineColor(kBlack);
      zero_line->Draw();

      PlottingUtils::SaveFigure(canvas, Form("double_waveform_fit_%d", i),
                                "double_waveform_fits",
                                PlotSaveOptions::kLINEAR);

      delete canvas;
      delete original;
      delete fitted;
      delete residual_graph;
      delete marker_peak1;
      delete marker_peak2;
    }
  }

  delete minimizer;

  input->Close();
  template_file->Close();
  calib_file->Close();

  output->cd();
  output_tree->Write();
  output->Close();
}

void PlotDoublePeaks(const TString &fitted_name = "fitted_doubles.root",
                     const TString &plot_prefix = "",
                     Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString output_filepath = "root_files/" + fitted_name;
  TFile *output = new TFile(output_filepath, "UPDATE");
  TTree *features_tree = static_cast<TTree *>(output->Get("features"));

  Float_t peak1_light_output, peak2_light_output;

  features_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  features_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);

  const Int_t peak1_nbins = 150, peak1_min = 100, peak1_max = 200;
  const Int_t peak2_nbins = 150, peak2_min = 0, peak2_max = 50;
  const Float_t peak1_bin_width =
      (Float_t)(peak1_max - peak1_min) / peak1_nbins;
  const Float_t peak2_bin_width =
      (Float_t)(peak2_max - peak2_min) / peak2_nbins;
  TH2F *LO1vsLO2 = new TH2F(
      "",
      Form("; Peak 1 Light Output [keVee]; Peak 2 Light Output [keVee]; "
           "Counts / %.2f keV / %.2f keV",
           peak1_bin_width, peak2_bin_width),
      peak1_nbins, peak1_min, peak1_max, peak2_nbins, peak2_min, peak2_max);
  TH1F *peak1_hist = new TH1F(
      "",
      Form("; Peak 1 Light Output [keVee]; Counts / %.2f keV", peak1_bin_width),
      peak1_nbins, peak1_min, peak1_max);
  TH1F *peak2_hist = new TH1F(
      "",
      Form("; Peak 2 Light Output [keVee]; Counts / %.2f keV", peak2_bin_width),
      peak2_nbins, peak2_min, peak2_max);

  Int_t num_entries = features_tree->GetEntries();
  for (Int_t j = 0; j < num_entries; j++) {
    features_tree->GetEntry(j);
    LO1vsLO2->Fill(peak1_light_output, peak2_light_output);
    peak1_hist->Fill(peak1_light_output);
    peak2_hist->Fill(peak2_light_output);
  }

  TCanvas *canvas2d = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(LO1vsLO2, canvas2d);
  PlottingUtils::SaveFigure(canvas2d,
                            plot_prefix + "peak1_vs_peak2_light_output", "",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvas1 = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(peak1_hist, kBlue + 1);
  PlottingUtils::SaveFigure(canvas1, plot_prefix + "peak1_light_output");

  TCanvas *canvas2 = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(peak2_hist, kRed + 1);
  PlottingUtils::SaveFigure(canvas2, plot_prefix + "peak2_light_output");

  output->cd();
  LO1vsLO2->Write("Peak1 vs Peak2 Light Output", TObject::kOverwrite);
  peak1_hist->Write("Peak1 Light Output", TObject::kOverwrite);
  peak2_hist->Write("Peak2 Light Output", TObject::kOverwrite);

  delete canvas2d;
  delete canvas1;
  delete canvas2;
  delete LO1vsLO2;
  delete peak1_hist;
  delete peak2_hist;

  output->Close();
  delete output;
}

struct HalfLifeResult {
  Double_t value = 0, error = 0;
  Bool_t valid = kFALSE;
};

struct ClusterCut {
  Double_t mean_x = 0, sigma_x = 1, mean_y = 0, sigma_y = 1;
  TF2 *gaussian2d = nullptr;
};

const Int_t TIME_LOWER_NS = 0, TIME_UPPER_NS = 300, TIME_BIN_WIDTH_NS = 2;
// 100 ns, not 65: on simulated doubles the waveform fit reshapes the
// distribution below ~100 ns and the exponential comes out 0.2 ns long there;
// from 100 ns the input half-life is recovered. Above ~200 ns the second pulse
// runs off the end of the recorded window and the acceptance falls.
const Int_t FIT_LOWER_NS = 100, FIT_UPPER_NS = 200;
const Double_t NOMINAL_CUT_SIGMA = 2.0;

ClusterCut FitCluster(TH2F *LO1vsLO2) {
  ClusterCut cut;
  cut.gaussian2d =
      new TF2("gaussian2d", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",
              0, 200, 0, 200);
  cut.gaussian2d->SetParameters(LO1vsLO2->GetMaximum(), LO1vsLO2->GetMean(1),
                                LO1vsLO2->GetRMS(1), LO1vsLO2->GetMean(2),
                                LO1vsLO2->GetRMS(2));
  LO1vsLO2->Fit(cut.gaussian2d, "Q");
  cut.mean_x = cut.gaussian2d->GetParameter(1);
  cut.sigma_x = TMath::Abs(cut.gaussian2d->GetParameter(2));
  cut.mean_y = cut.gaussian2d->GetParameter(3);
  cut.sigma_y = TMath::Abs(cut.gaussian2d->GetParameter(4));
  std::cout << "2D Gaussian fit: " << std::endl;
  std::cout << "Peak 1 mean: " << cut.mean_x << " +/- " << cut.sigma_x
            << std::endl;
  std::cout << "Peak 2 mean: " << cut.mean_y << " +/- " << cut.sigma_y
            << std::endl;
  return cut;
}

// Time differences of the events within cut_sigma of the cluster centre, in
// units of the fitted widths.
TH1F *SelectTimeDifferences(TTree *features_tree, const ClusterCut &cut,
                            Double_t cut_sigma) {
  Float_t peak1_light_output, peak2_light_output, time_difference;
  features_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  features_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);
  features_tree->SetBranchAddress("time_difference", &time_difference);

  TH1F *time_diff_hist = new TH1F(
      "", Form("; Time Difference [ns]; Counts / %d ns", TIME_BIN_WIDTH_NS),
      (TIME_UPPER_NS - TIME_LOWER_NS) / TIME_BIN_WIDTH_NS, TIME_LOWER_NS,
      TIME_UPPER_NS);
  time_diff_hist->Sumw2();

  Int_t num_entries = features_tree->GetEntries();
  for (Int_t j = 0; j < num_entries; j++) {
    features_tree->GetEntry(j);
    Double_t dist_x = (peak1_light_output - cut.mean_x) / cut.sigma_x;
    Double_t dist_y = (peak2_light_output - cut.mean_y) / cut.sigma_y;
    if (TMath::Sqrt(dist_x * dist_x + dist_y * dist_y) < cut_sigma)
      time_diff_hist->Fill(time_difference);
  }
  return time_diff_hist;
}

// Likelihood, not chi2: the tail bins hold tens of counts, and a chi2 fit
// with errors taken from the bin contents overweights downward fluctuations
// there. On simulated doubles that shortened the recovered half-life by
// 0.3 ns; the likelihood fit recovers the input.
HalfLifeResult FitExponential(TH1F *time_diff_hist, TF1 *exponential) {
  HalfLifeResult result;
  exponential->SetRange(FIT_LOWER_NS, FIT_UPPER_NS);
  exponential->SetParameters(time_diff_hist->GetMaximum(), 100);
  exponential->SetParNames("Amplitude", "Lifetime");
  TFitResultPtr fit = time_diff_hist->Fit(exponential, "RLQS");
  if (!fit.Get() || !fit->IsValid())
    return result;
  result.value = exponential->GetParameter(1) * TMath::Log(2);
  result.error = exponential->GetParError(1) * TMath::Log(2);
  // Minuit reports success on fits that ran away on a starved histogram;
  // only a half-life in the physical range with a finite, smaller error
  // counts.
  result.valid = TMath::Finite(result.value) && TMath::Finite(result.error) &&
                 result.value > 1.0 && result.value < 100.0 &&
                 result.error > 0 && result.error < result.value;
  return result;
}

HalfLifeResult
FitAndExtractHalfLife(const TString &fitted_name = "fitted_doubles.root",
                      const TString &results_name = "half_life_results.root",
                      const TString &plot_prefix = "",
                      Bool_t reprocess = kFALSE,
                      Double_t cut_sigma = NOMINAL_CUT_SIGMA) {
  HalfLifeResult result;
  if (!reprocess)
    return result;
  TString input_filepath = "root_files/" + fitted_name;
  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));
  TH2F *LO1vsLO2 =
      static_cast<TH2F *>(input->Get("Peak1 vs Peak2 Light Output"));

  ClusterCut cut = FitCluster(LO1vsLO2);
  TH1F *time_diff_hist = SelectTimeDifferences(features_tree, cut, cut_sigma);
  TF1 *exponential = new TF1("exponential", "[0]*TMath::Exp(-x/[1])",
                             FIT_LOWER_NS, FIT_UPPER_NS);
  result = FitExponential(time_diff_hist, exponential);
  std::cout << "Lifetime: " << exponential->GetParameter(1) << " +/- "
            << exponential->GetParError(1) << " ns" << std::endl;
  std::cout << "Half-life: " << result.value << " +/- " << result.error
            << " ns   (cut " << cut_sigma << " sigma)" << std::endl;

  TCanvas *canvas_2d = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(LO1vsLO2, canvas_2d);
  cut.gaussian2d->Draw("CONT3 SAME");
  PlottingUtils::SaveFigure(canvas_2d, plot_prefix + "fitted_2d_gaussian", "",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvas_time = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(time_diff_hist, kBlue + 1);
  time_diff_hist->SetMarkerSize(2);
  time_diff_hist->Draw("E SAME");
  TMarker *marker_lower =
      new TMarker(FIT_LOWER_NS, exponential->Eval(FIT_LOWER_NS), 21);
  marker_lower->SetMarkerColor(kRed);
  marker_lower->SetMarkerSize(1);
  marker_lower->Draw();

  Float_t marker_upper_y =
      exponential->Eval(FIT_UPPER_NS) > time_diff_hist->GetMinimum()
          ? exponential->Eval(FIT_UPPER_NS)
          : time_diff_hist->GetMinimum();
  TMarker *marker_upper = new TMarker(FIT_UPPER_NS, marker_upper_y, 21);
  marker_upper->SetMarkerColor(kRed);
  marker_upper->SetMarkerSize(1);
  marker_upper->Draw();
  exponential->SetRange(TIME_LOWER_NS, TIME_UPPER_NS);
  exponential->Draw("SAME");

  TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  leg->AddEntry((TObject *)0,
                Form("t_{1/2} = %.1f #pm %.1f ns", result.value, result.error),
                "");
  leg->AddEntry(marker_upper, "Fit region");
  leg->SetMargin(0.1);
  leg->Draw();

  PlottingUtils::SaveFigure(canvas_time, plot_prefix + "time_difference_fit");

  TString output_filepath = "root_files/" + results_name;
  TFile *output = new TFile(output_filepath, "RECREATE");
  LO1vsLO2->Write("2D_Gaussian_Fit");
  time_diff_hist->Write("Time_Difference");
  cut.gaussian2d->Write("gaussian2d_fit");
  exponential->Write("exponential_fit");
  output->Close();

  delete canvas_2d;
  delete canvas_time;
  delete time_diff_hist;
  delete cut.gaussian2d;
  delete exponential;
  delete output;

  input->Close();
  delete input;
  return result;
}

// The quoted uncertainty is the sensitivity of the half-life to where the
// fiducial cut is drawn, not the parabolic error of one fit: refit with the
// cut radius stepped from 0.5 to 5 sigma (below that the fit is noise-limited,
// above it the selection no longer changes) and take the spread. This is the
// perturbation analysis of Fondement et al., Nucl. Phys. A 1075 (2026) 123485.
HalfLifeResult ScanCutWidth(const TString &fitted_name = "fitted_doubles.root",
                            const TString &results_name = "cut_scan.root",
                            const TString &plot_prefix = "",
                            Bool_t reprocess = kFALSE) {
  HalfLifeResult summary;
  if (!reprocess)
    return summary;
  TString input_filepath = "root_files/" + fitted_name;
  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));
  TH2F *LO1vsLO2 =
      static_cast<TH2F *>(input->Get("Peak1 vs Peak2 Light Output"));
  ClusterCut cut = FitCluster(LO1vsLO2);

  const Double_t k_min = 0.5, k_max = 5.0, k_step = 0.1;
  std::vector<Double_t> ks, values, errors;
  TF1 *exponential = new TF1("exponential_scan", "[0]*TMath::Exp(-x/[1])",
                             FIT_LOWER_NS, FIT_UPPER_NS);
  std::cout << "Cut-width scan:" << std::endl;
  for (Double_t k = k_min; k <= k_max + 1e-9; k += k_step) {
    TH1F *time_diff_hist = SelectTimeDifferences(features_tree, cut, k);
    HalfLifeResult r = FitExponential(time_diff_hist, exponential);
    if (r.valid) {
      ks.push_back(k);
      values.push_back(r.value);
      errors.push_back(r.error);
      std::cout << "  k = " << std::fixed << std::setprecision(1) << k
                << "   t1/2 = " << std::setprecision(3) << r.value << " +/- "
                << r.error << " ns   (" << (Int_t)time_diff_hist->GetEntries()
                << " events)" << std::endl;
    }
    delete time_diff_hist;
  }
  delete exponential;

  Int_t n = (Int_t)ks.size();
  if (n == 0) {
    input->Close();
    delete input;
    return summary;
  }
  Double_t sw = 0, swx = 0, sx = 0, sxx = 0;
  for (Int_t i = 0; i < n; i++) {
    Double_t w = 1.0 / (errors[i] * errors[i]);
    sw += w;
    swx += w * values[i];
    sx += values[i];
    sxx += values[i] * values[i];
  }
  Double_t mean = swx / sw;
  Double_t spread = TMath::Sqrt(TMath::Max(0.0, sxx / n - (sx / n) * (sx / n)));
  Double_t lo = *std::min_element(values.begin(), values.end());
  Double_t hi = *std::max_element(values.begin(), values.end());
  summary.value = mean;
  summary.error = spread;
  summary.valid = kTRUE;
  std::cout << "Cut-width scan summary (" << n << " cuts from " << k_min
            << " to " << k_max << " sigma):" << std::endl;
  std::cout << "  weighted mean t1/2 = " << std::setprecision(3) << mean
            << " ns   spread (rms) = " << spread << " ns   range " << lo
            << " to " << hi << " ns" << std::endl;

  TGraphErrors *graph =
      new TGraphErrors(n, ks.data(), values.data(), nullptr, errors.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(graph, kBlue + 1,
                                "; Cut radius [#sigma]; t_{1/2} [ns]");
  graph->Draw("APE");
  TLine *mean_line = new TLine(k_min, mean, k_max, mean);
  mean_line->SetLineColor(kRed);
  mean_line->SetLineStyle(2);
  mean_line->Draw();
  PlottingUtils::SaveFigure(canvas, plot_prefix + "half_life_vs_cut_width", "",
                            PlotSaveOptions::kLINEAR);

  TString output_filepath = "root_files/" + results_name;
  TFile *output = new TFile(output_filepath, "RECREATE");
  graph->Write("half_life_vs_cut_width");
  output->Close();

  delete canvas;
  delete mean_line;
  delete graph;
  delete cut.gaussian2d;
  delete output;
  input->Close();
  delete input;
  return summary;
}

void HalfLife() {
  InitUtils::SetROOTPreferences();

  Bool_t reprocess_initial = kTRUE;
  Bool_t reprocess_candidate = kTRUE;
  Bool_t reprocess_waveform = kTRUE;
  Bool_t reprocess_halflife = kTRUE;

  std::vector<TString> input_names = Constants::ALL_OUTPUT_NAMES;

  CalculatePSPvsLO(input_names, reprocess_initial);
  if (reprocess_initial)
    PlotPSPvsLO(input_names);

  std::vector<TString> irradiation_names = Constants::IRRADIATION_DATASETS;

  GetCandidateWaveforms(irradiation_names, reprocess_candidate);
  AnalyzeDoubleWaveforms(reprocess_waveform);
  if (reprocess_waveform)
    PlotPSPvsLO({"double_waveforms"});
  CreateTemplateWaveforms(reprocess_waveform);
  FitDoubleWaveforms("double_waveforms.root", "fitted_doubles.root",
                     reprocess_waveform);
  PlotDoublePeaks("fitted_doubles.root", "", reprocess_halflife);
  FitAndExtractHalfLife("fitted_doubles.root", "half_life_results.root", "",
                        reprocess_halflife);
  ScanCutWidth("fitted_doubles.root", "cut_scan.root", "", reprocess_halflife);
}
