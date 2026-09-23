#include "HalfLife.cpp"
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
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TSystem.h>
#include <TTree.h>
#include <vector>

// Simulate the half-life the data gave, read back from the extraction, so the
// closure asks whether this measurement reproduces itself.
Double_t MeasuredHalfLife() {
  TFile file("root_files/half_life_results.root", "READ");
  TF1 *exponential = static_cast<TF1 *>(file.Get("exponential_fit"));
  if (!exponential) {
    std::cerr << "No exponential_fit in half_life_results.root; run HalfLife "
                 "first"
              << std::endl;
    return 0;
  }
  return exponential->GetParameter(1) * TMath::Log(2);
}

void GenerateMCData(Int_t num_events, Double_t half_life,
                    Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  std::cout << "Generating Monte Carlo data..." << std::endl;

  TString calibration_filepath = "root_files/calibration_function.root";
  TFile *calib_file = new TFile(calibration_filepath, "READ");
  TF1 *calibration_function =
      static_cast<TF1 *>(calib_file->Get("calibration"));

  if (!calibration_function) {
    std::cerr << "Error: Could not load calibration function!" << std::endl;
    return;
  }

  Double_t calibration_slope = calibration_function->GetParameter(1);
  Double_t calibration_intercept = calibration_function->GetParameter(0);

  TString template_filepath = "root_files/template_waveforms.root";
  TFile *template_file = new TFile(template_filepath, "READ");
  TGraph *template_graph_148 =
      static_cast<TGraph *>(template_file->Get("template_148keV"));
  TGraph *template_graph_32 =
      static_cast<TGraph *>(template_file->Get("template_32keV"));

  if (!template_graph_148 || !template_graph_32) {
    std::cerr << "Error: Could not load template waveforms!" << std::endl;
    return;
  }

  Double_t *template_148_y = template_graph_148->GetY();
  Double_t *template_32_y = template_graph_32->GetY();
  Int_t template_npoints = template_graph_148->GetN();

  Int_t template_148_peak_pos = TMath::LocMax(template_npoints, template_148_y);
  Int_t template_32_peak_pos = TMath::LocMax(template_npoints, template_32_y);
  Double_t template_148_peak_height = template_148_y[template_148_peak_pos];
  Double_t template_32_peak_height = template_32_y[template_32_peak_pos];

  Double_t true_light_output_148 = 148.55;
  Double_t true_light_output_32 = 32.3;
  Double_t sigma_adc_148 = 51.2089;
  Double_t sigma_adc_32 = 34.3313;
  Double_t lifetime = half_life / TMath::Log(2);

  Double_t sigma_keV_148 = calibration_slope * sigma_adc_148;
  Double_t sigma_keV_32 = calibration_slope * sigma_adc_32;

  std::cout << "MC Parameters:" << std::endl;
  std::cout << "  Peak 1: " << true_light_output_148
            << " keVee, sigma = " << sigma_keV_148 << " keVee" << std::endl;
  std::cout << "  Peak 2: " << true_light_output_32
            << " keVee, sigma = " << sigma_keV_32 << " keVee" << std::endl;
  std::cout << "  Half-life: " << half_life << " ns (lifetime = " << lifetime
            << " ns)" << std::endl;
  std::cout << "  Number of events: " << num_events << std::endl;

  TString output_filepath = "root_files/mc_double_waveforms.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree = new TTree("features", "MC double waveform features.");

  TArrayF *output_samples = nullptr;
  Float_t output_light_output;
  Int_t peak1_position, peak2_position;
  Float_t true_peak1_light_output, true_peak2_light_output, true_time_diff;

  output_tree->Branch("light_output", &output_light_output, "light_output/F");
  output_tree->Branch("peak1_position", &peak1_position, "peak1_position/I");
  output_tree->Branch("peak2_position", &peak2_position, "peak2_position/I");
  output_tree->Branch("Samples", &output_samples);
  output_tree->Branch("true_peak1_light_output", &true_peak1_light_output,
                      "true_peak1_light_output/F");
  output_tree->Branch("true_peak2_light_output", &true_peak2_light_output,
                      "true_peak2_light_output/F");
  output_tree->Branch("true_time_diff", &true_time_diff, "true_time_diff/F");

  TRandom3 *rng = new TRandom3(0);
  Int_t rejected_light_output = 0, rejected_peak_finding = 0;

  for (Int_t i = 0; i < num_events; i++) {
    if (i % 10000 == 0) {
      std::cout << "Generating event " << i << "/" << num_events << std::endl;
    }

    Double_t sampled_light_output_148 =
        rng->Gaus(true_light_output_148, sigma_keV_148);
    Double_t sampled_light_output_32 =
        rng->Gaus(true_light_output_32, sigma_keV_32);

    if (sampled_light_output_148 < 0)
      sampled_light_output_148 = 0.1;
    if (sampled_light_output_32 < 0)
      sampled_light_output_32 = 0.1;

    Double_t sampled_time_diff = rng->Exp(lifetime);

    if (sampled_time_diff > 300.0) {
      i--;
      continue;
    }

    true_peak1_light_output = sampled_light_output_148;
    true_peak2_light_output = sampled_light_output_32;
    true_time_diff = sampled_time_diff;

    Double_t ph_148 = sampled_light_output_148 / calibration_slope;
    Double_t ph_32 = sampled_light_output_32 / calibration_slope;

    Double_t scale_148 = ph_148 / template_148_peak_height;
    Double_t scale_32 = ph_32 / template_32_peak_height;

    Int_t nsamples = template_npoints;
    std::vector<Float_t> waveform(nsamples, 0);

    Int_t first_peak_sample_pos = template_148_peak_pos;
    // The second pulse lands at a continuous time, not on the sample grid:
    // place it by the same linear interpolation the fit model uses. Snapping
    // to the nearest sample put every fitted time on an even number of ns and
    // produced a comb in the time-difference histogram.
    Double_t second_peak_pos_exact =
        first_peak_sample_pos + sampled_time_diff / 2.0;

    for (Int_t j = 0; j < nsamples; j++) {
      Int_t template_index = j - first_peak_sample_pos + template_148_peak_pos;
      if (template_index >= 0 && template_index < template_npoints) {
        waveform[j] += (Float_t)(template_148_y[template_index] * scale_148);
      }
    }

    for (Int_t j = 0; j < nsamples; j++) {
      Double_t exact_index = j - second_peak_pos_exact + template_32_peak_pos;
      if (exact_index < 0 || exact_index >= template_npoints - 1)
        continue;
      Int_t idx_low = (Int_t)exact_index;
      Double_t frac = exact_index - idx_low;
      Double_t value = (1.0 - frac) * template_32_y[idx_low] +
                       frac * template_32_y[idx_low + 1];
      waveform[j] += (Float_t)(value * scale_32);
    }

    Double_t max_value = 0;
    for (Int_t j = 0; j < nsamples; j++) {
      if (waveform[j] > max_value) {
        max_value = waveform[j];
      }
    }

    output_light_output = calibration_function->Eval(max_value);

    // Same selection as the real candidates: the light-output window of
    // GetCandidateWaveforms, then TSpectrum must resolve two peaks, whose
    // positions seed the fit exactly as they do for data. Overlapping pulses
    // the finder cannot separate never reach the fit in data, so they must not
    // reach it here either.
    if (output_light_output < CANDIDATE_LO_MIN_KEVEE ||
        output_light_output > CANDIDATE_LO_MAX_KEVEE) {
      rejected_light_output++;
      continue;
    }
    std::vector<Double_t> waveform_array(waveform.begin(), waveform.end());
    if (!FindDoublePeaks(waveform_array, peak1_position, peak2_position)) {
      rejected_peak_finding++;
      continue;
    }

    if (output_samples)
      delete output_samples;
    output_samples = new TArrayF(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      output_samples->SetAt(waveform[j], j);
    }

    output_tree->Fill();
  }

  std::cout << "Generated " << num_events
            << " MC events: " << output_tree->GetEntries()
            << " pass the candidate selection (" << rejected_light_output
            << " outside the light-output window, " << rejected_peak_finding
            << " without two resolved peaks)" << std::endl;

  output->cd();
  output_tree->Write("", TObject::kOverwrite);
  output->Close();

  template_file->Close();
  calib_file->Close();

  delete rng;
  delete output;
  delete template_file;
  delete calib_file;

  std::cout << "MC data saved to " << output_filepath << std::endl;
}

void AnalyzeMCEfficiency(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  std::cout << "Analyzing MC efficiency..." << std::endl;

  TString mc_input_filepath = "root_files/mc_double_waveforms.root";
  TFile *mc_input = new TFile(mc_input_filepath, "READ");
  TTree *mc_input_tree = static_cast<TTree *>(mc_input->Get("features"));

  Int_t total_generated = mc_input_tree->GetEntries();
  std::cout << "Total MC events generated: " << total_generated << std::endl;

  mc_input->Close();
  delete mc_input;

  FitDoubleWaveforms("mc_double_waveforms.root", "mc_fitted_doubles.root",
                     reprocess);

  TString mc_fitted_filepath = "root_files/mc_fitted_doubles.root";
  TFile *mc_fitted = new TFile(mc_fitted_filepath, "READ");
  TTree *mc_fitted_tree = static_cast<TTree *>(mc_fitted->Get("features"));

  Int_t successfully_fitted = mc_fitted_tree->GetEntries();
  std::cout << "Successfully fitted events: " << successfully_fitted
            << std::endl;

  Float_t peak1_light_output, peak2_light_output;
  mc_fitted_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  mc_fitted_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);

  TH2F *LO1vsLO2 =
      new TH2F("", "; Peak 1 Light Output [keVee]; Peak 2 Light Output [keVee]",
               150, 100, 200, 150, 0, 50);

  for (Int_t i = 0; i < successfully_fitted; i++) {
    mc_fitted_tree->GetEntry(i);
    LO1vsLO2->Fill(peak1_light_output, peak2_light_output);
  }

  TF2 *gaussian2d =
      new TF2("gaussian2d", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",
              0, 200, 0, 200);
  gaussian2d->SetParameters(LO1vsLO2->GetMaximum(), LO1vsLO2->GetMean(1),
                            LO1vsLO2->GetRMS(1), LO1vsLO2->GetMean(2),
                            LO1vsLO2->GetRMS(2));

  LO1vsLO2->Fit(gaussian2d, "Q");

  Double_t mean_x = gaussian2d->GetParameter(1);
  Double_t sigma_x = TMath::Abs(gaussian2d->GetParameter(2));
  Double_t mean_y = gaussian2d->GetParameter(3);
  Double_t sigma_y = TMath::Abs(gaussian2d->GetParameter(4));

  std::cout << "2D Gaussian fit: " << std::endl;
  std::cout << "Peak 1 mean: " << mean_x << " +/- " << sigma_x << std::endl;
  std::cout << "Peak 2 mean: " << mean_y << " +/- " << sigma_y << std::endl;

  Int_t events_after_2d_cut = 0;
  mc_fitted_tree->GetEntry(0);
  for (Int_t i = 0; i < successfully_fitted; i++) {
    mc_fitted_tree->GetEntry(i);

    Double_t dist_x = (peak1_light_output - mean_x) / sigma_x;
    Double_t dist_y = (peak2_light_output - mean_y) / sigma_y;
    Double_t distance = TMath::Sqrt(dist_x * dist_x + dist_y * dist_y);

    if (distance < 2.0) {
      events_after_2d_cut++;
    }
  }

  std::cout << "Events passing 2D Gaussian cut (2 sigma): "
            << events_after_2d_cut << std::endl;

  mc_fitted->Close();
  delete mc_fitted;

  Double_t fitting_efficiency =
      (Double_t)successfully_fitted / (Double_t)total_generated;
  Double_t cut_efficiency =
      (Double_t)events_after_2d_cut / (Double_t)successfully_fitted;
  Double_t total_efficiency =
      (Double_t)events_after_2d_cut / (Double_t)total_generated;

  std::cout << "\n=== Efficiency Results ===" << std::endl;
  std::cout << "Fitting efficiency: " << fitting_efficiency * 100 << "% ("
            << successfully_fitted << "/" << total_generated << ")"
            << std::endl;
  std::cout << "2D Gaussian cut efficiency: " << cut_efficiency * 100 << "% ("
            << events_after_2d_cut << "/" << successfully_fitted << ")"
            << std::endl;
  std::cout << "Total efficiency: " << total_efficiency * 100 << "% ("
            << events_after_2d_cut << "/" << total_generated << ")"
            << std::endl;

  TString efficiency_filepath = "root_files/mc_efficiency.root";
  TFile *efficiency_file = new TFile(efficiency_filepath, "RECREATE");

  TNamed *total_gen =
      new TNamed("total_generated", Form("%d", total_generated));
  TNamed *fitted =
      new TNamed("successfully_fitted", Form("%d", successfully_fitted));
  TNamed *after_cut =
      new TNamed("events_after_2d_cut", Form("%d", events_after_2d_cut));
  TNamed *eff_fit =
      new TNamed("fitting_efficiency", Form("%.6f", fitting_efficiency));
  TNamed *eff_cut = new TNamed("cut_efficiency", Form("%.6f", cut_efficiency));
  TNamed *eff_total =
      new TNamed("total_efficiency", Form("%.6f", total_efficiency));

  total_gen->Write();
  fitted->Write();
  after_cut->Write();
  eff_fit->Write();
  eff_cut->Write();
  eff_total->Write();

  LO1vsLO2->Write("2D_Light_Output");
  gaussian2d->Write("gaussian2d_fit");

  efficiency_file->Close();
  delete efficiency_file;
  delete LO1vsLO2;
  delete gaussian2d;

  std::cout << "Efficiency results saved to " << efficiency_filepath
            << std::endl;
}

// The point of the Monte Carlo: run the SAME 2D cut and exponential fit that
// produce the measured half-life on simulated doubles with a known half-life,
// and see whether it comes back.
void CheckHalfLifeClosure(Double_t half_life, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;
  PlotDoublePeaks("mc_fitted_doubles.root", "mc_", reprocess);
  HalfLifeResult recovered = FitAndExtractHalfLife(
      "mc_fitted_doubles.root", "mc_half_life_results.root", "mc_", reprocess);
  if (!recovered.valid) {
    std::cerr << "MC half-life extraction failed" << std::endl;
    return;
  }
  Double_t pull = (recovered.error > 0)
                      ? (recovered.value - half_life) / recovered.error
                      : 0.0;
  std::cout << std::endl << "=== MC half-life closure ===" << std::endl;
  std::cout << "  simulated: " << half_life << " ns" << std::endl;
  std::cout << "  recovered: " << recovered.value << " +/- " << recovered.error
            << " ns   (" << recovered.value - half_life << " ns, " << pull
            << " sigma)" << std::endl;
  HalfLifeResult scan = ScanCutWidth("mc_fitted_doubles.root",
                                     "mc_cut_scan.root", "mc_", reprocess);
  if (scan.valid)
    std::cout << "  cut-width scan: " << scan.value << " +/- " << scan.error
              << " ns (weighted mean +/- spread), input " << half_life
              << " ns is " << (scan.value - half_life) / scan.error
              << " spreads away" << std::endl;
}

void PlotTrueMC(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  std::cout << "Plotting MC comparison..." << std::endl;

  TString mc_input_filepath = "root_files/mc_double_waveforms.root";
  TFile *mc_input = new TFile(mc_input_filepath, "READ");
  TTree *mc_input_tree = static_cast<TTree *>(mc_input->Get("features"));

  Float_t true_peak1_light_output, true_peak2_light_output, true_time_diff;
  mc_input_tree->SetBranchAddress("true_peak1_light_output",
                                  &true_peak1_light_output);
  mc_input_tree->SetBranchAddress("true_peak2_light_output",
                                  &true_peak2_light_output);
  mc_input_tree->SetBranchAddress("true_time_diff", &true_time_diff);

  TH2F *true_LO1vsLO2 = new TH2F(
      "",
      "; True Peak 1 Light Output [keVee]; True Peak 2 Light Output [keVee]",
      150, 100, 200, 150, 0, 50);
  TH1F *true_peak1_hist =
      new TH1F("", "; True Peak 1 Light Output [keVee]; Counts", 150, 100, 200);
  TH1F *true_peak2_hist =
      new TH1F("", "; True Peak 2 Light Output [keVee]; Counts", 150, 0, 50);
  TH1F *true_time_diff_hist =
      new TH1F("", "; True Time Difference [ns]; Counts / 2 ns", 150, 0, 300);

  Int_t num_entries = mc_input_tree->GetEntries();
  for (Int_t i = 0; i < num_entries; i++) {
    mc_input_tree->GetEntry(i);
    true_LO1vsLO2->Fill(true_peak1_light_output, true_peak2_light_output);
    true_peak1_hist->Fill(true_peak1_light_output);
    true_peak2_hist->Fill(true_peak2_light_output);
    true_time_diff_hist->Fill(true_time_diff);
  }

  TCanvas *canvas_2d_true = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(true_LO1vsLO2, canvas_2d_true);
  PlottingUtils::SaveFigure(canvas_2d_true,
                            "mc_true_peak1_vs_peak2_light_output", "",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvas_p1_true = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(true_peak1_hist, kBlue + 1);
  PlottingUtils::SaveFigure(canvas_p1_true, "mc_true_peak1_light_output");

  TCanvas *canvas_p2_true = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(true_peak2_hist, kRed + 1);
  PlottingUtils::SaveFigure(canvas_p2_true, "mc_true_peak2_light_output");

  TCanvas *canvas_time_true = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(true_time_diff_hist, kGreen + 1);
  PlottingUtils::SaveFigure(canvas_time_true, "mc_true_time_difference");

  mc_input->Close();

  TString mc_fitted_filepath = "root_files/mc_fitted_doubles.root";
  TFile *mc_fitted = new TFile(mc_fitted_filepath, "READ");
  TTree *mc_fitted_tree = static_cast<TTree *>(mc_fitted->Get("features"));

  Float_t peak1_light_output, peak2_light_output, time_difference;
  mc_fitted_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  mc_fitted_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);
  mc_fitted_tree->SetBranchAddress("time_difference", &time_difference);

  TH2F *fitted_LO1vsLO2 = new TH2F("",
                                   "; Fitted Peak 1 Light Output [keVee]; "
                                   "Fitted Peak 2 Light Output [keVee]",
                                   150, 100, 200, 150, 0, 50);
  TH1F *fitted_peak1_hist = new TH1F(
      "", "; Fitted Peak 1 Light Output [keVee]; Counts", 150, 100, 200);
  TH1F *fitted_peak2_hist =
      new TH1F("", "; Fitted Peak 2 Light Output [keVee]; Counts", 150, 0, 50);
  TH1F *fitted_time_diff_hist =
      new TH1F("", "; Fitted Time Difference [ns]; Counts / 2 ns", 150, 0, 300);

  Int_t num_fitted_entries = mc_fitted_tree->GetEntries();
  for (Int_t i = 0; i < num_fitted_entries; i++) {
    mc_fitted_tree->GetEntry(i);
    fitted_LO1vsLO2->Fill(peak1_light_output, peak2_light_output);
    fitted_peak1_hist->Fill(peak1_light_output);
    fitted_peak2_hist->Fill(peak2_light_output);
    fitted_time_diff_hist->Fill(time_difference);
  }

  TCanvas *canvas_2d_fitted = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(fitted_LO1vsLO2, canvas_2d_fitted);
  PlottingUtils::SaveFigure(canvas_2d_fitted,
                            "mc_fitted_peak1_vs_peak2_light_output", "",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvas_p1_fitted = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(fitted_peak1_hist, kBlue + 1);
  PlottingUtils::SaveFigure(canvas_p1_fitted, "mc_fitted_peak1_light_output");

  TCanvas *canvas_p2_fitted = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(fitted_peak2_hist, kRed + 1);
  PlottingUtils::SaveFigure(canvas_p2_fitted, "mc_fitted_peak2_light_output");

  TCanvas *canvas_time_fitted = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(fitted_time_diff_hist, kGreen + 1);
  PlottingUtils::SaveFigure(canvas_time_fitted, "mc_fitted_time_difference");

  mc_fitted->Close();

  delete canvas_2d_true;
  delete canvas_p1_true;
  delete canvas_p2_true;
  delete canvas_time_true;
  delete canvas_2d_fitted;
  delete canvas_p1_fitted;
  delete canvas_p2_fitted;
  delete canvas_time_fitted;

  delete true_LO1vsLO2;
  delete true_peak1_hist;
  delete true_peak2_hist;
  delete true_time_diff_hist;
  delete fitted_LO1vsLO2;
  delete fitted_peak1_hist;
  delete fitted_peak2_hist;
  delete fitted_time_diff_hist;

  delete mc_input;
  delete mc_fitted;

  std::cout << "MC comparison plots saved." << std::endl;
}
void MonteCarloEfficiency() {
  InitUtils::SetROOTPreferences();
  Double_t half_life = MeasuredHalfLife();
  if (!(half_life > 0))
    return;
  GenerateMCData(5000000, half_life, kTRUE);
  AnalyzeMCEfficiency(kTRUE);
  CheckHalfLifeClosure(half_life, kTRUE);
  PlotTrueMC(kTRUE);
}
