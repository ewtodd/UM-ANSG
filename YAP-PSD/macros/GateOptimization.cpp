#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

struct OptimizationParameters {
  const Float_t min_light_output = 900.0;
  const Float_t max_light_output = 1200.0;

  const Int_t min_short_gate = 5;
  const Int_t max_short_gate = 40;
  const Int_t short_gate_step = 1;
  const Int_t min_long_gate = 45;
  const Int_t max_long_gate = 220;
  const Int_t long_gate_step = 5;
};

struct OptimizationResult {
  Int_t short_gate;
  Int_t long_gate;
  Double_t fom;
  Double_t mean_alpha;
  Double_t sigma_alpha;
  Double_t amp_alpha;
  Double_t mean_gamma;
  Double_t sigma_gamma;
  Double_t amp_gamma;
};

struct Event {
  std::vector<Float_t> samples;
  Float_t light_output;
  Int_t trigger_position;
};

struct GateCombination {
  Int_t short_gate;
  Int_t long_gate;
};

struct ComputeResult {
  Int_t short_gate;
  Int_t long_gate;
  std::vector<Float_t> alpha_cc;
  std::vector<Float_t> gamma_cc;
};

const OptimizationParameters OPTIMIZATION_PARAMETERS;

std::vector<Event> LoadEvents(const TString &output_name) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString filepath = project_root + "/root_files/" + output_name + ".root";
  TFile *file = new TFile(filepath, "READ");
  TTree *tree = static_cast<TTree *>(file->Get("features"));

  TArrayF *samples = nullptr;
  Float_t light_output_keVee;
  Int_t trigger_position;

  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("light_output", &light_output_keVee);
  tree->SetBranchAddress("trigger_position", &trigger_position);
  tree->LoadBaskets();

  Int_t n_entries = tree->GetEntries();
  std::vector<Event> events;
  events.reserve(n_entries);

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    if (light_output_keVee <= OPTIMIZATION_PARAMETERS.min_light_output ||
        light_output_keVee >= OPTIMIZATION_PARAMETERS.max_light_output)
      continue;

    Event evt;
    evt.light_output = light_output_keVee;
    evt.trigger_position = trigger_position;

    Int_t n_samples = samples->GetSize();
    evt.samples.resize(n_samples);
    for (Int_t j = 0; j < n_samples; j++) {
      evt.samples[j] = samples->At(j);
    }

    events.push_back(std::move(evt));
  }

  file->Close();
  delete file;

  std::cout << "Loaded " << events.size() << " events from " << output_name
            << " (within light output window)" << std::endl;

  return events;
}

std::vector<Float_t> ComputeCCValues(const std::vector<Event> &events,
                                     Int_t short_gate, Int_t long_gate,
                                     Int_t pre_gate) {
  std::vector<Float_t> cc_values;
  cc_values.reserve(events.size());

  for (size_t i = 0; i < events.size(); i++) {
    Int_t n_samples = static_cast<Int_t>(events[i].samples.size());
    Int_t start = std::max<Int_t>(events[i].trigger_position - pre_gate, 0);
    Int_t short_end = std::min(start + short_gate, n_samples);
    Int_t long_end = std::min(start + long_gate, n_samples);

    Float_t short_sum = 0.0f;
    Float_t long_sum = 0.0f;

    for (Int_t j = start; j < long_end; j++) {
      Float_t value = events[i].samples[j];
      if (j < short_end)
        short_sum += value;
      long_sum += value;
    }

    if (long_sum > short_sum) {
      cc_values.push_back(1.0f - (short_sum / long_sum));
    }
  }

  return cc_values;
}

TH1F *BuildHistogram(const std::vector<Float_t> &cc_values) {
  TH1F *hist = new TH1F(PlottingUtils::GetRandomName(), ";CC PSD;Counts",
                        Constants::CC_HIST_NBINS, Constants::CC_HIST_XMIN,
                        Constants::CC_HIST_XMAX);
  hist->SetDirectory(0);

  for (size_t i = 0; i < cc_values.size(); i++) {
    hist->Fill(cc_values[i]);
  }

  return hist;
}

TF1 *FitGaussian(TH1F *hist) {
  if (!hist || hist->GetEntries() == 0)
    return nullptr;

  Int_t first_bin = hist->FindFirstBinAbove(hist->GetMaximum() * 0.01);
  Int_t last_bin = hist->FindLastBinAbove(hist->GetMaximum() * 0.01);
  Double_t fit_min = hist->GetBinCenter(first_bin);
  Double_t fit_max = hist->GetBinCenter(last_bin);

  TF1 *fit = new TF1(PlottingUtils::GetRandomName(), "gaus", fit_min, fit_max);

  Double_t amp = hist->GetMaximum();
  Double_t mean = hist->GetMean();
  Double_t rms = hist->GetRMS();

  fit->SetParameters(amp, mean, rms);
  fit->SetParLimits(0, amp * 0.1, amp * 2.0);
  fit->SetParLimits(1, 0, 1);
  fit->SetParLimits(2, 0, 1);

  Int_t status = hist->Fit(fit, "LRQSN");
  if (status != 0) {
    delete fit;
    return nullptr;
  }

  return fit;
}

Double_t
CalculateFigureOfMerit(TH1F *hist_alpha, TH1F *hist_gamma, Double_t &separation,
                       Double_t &fwhm_alpha, Double_t &fwhm_gamma,
                       Double_t &out_amp_alpha, Double_t &out_mean_alpha,
                       Double_t &out_sigma_alpha, Double_t &out_amp_gamma,
                       Double_t &out_mean_gamma, Double_t &out_sigma_gamma) {
  if (!hist_alpha || !hist_gamma)
    return -1;

  if (hist_alpha->GetEntries() == 0 || hist_gamma->GetEntries() == 0)
    return -1;

  TF1 *fit_alpha = FitGaussian(hist_alpha);
  TF1 *fit_gamma = FitGaussian(hist_gamma);

  if (!fit_alpha || !fit_gamma) {
    delete fit_alpha;
    delete fit_gamma;
    return -1;
  }

  Double_t mean_alpha = fit_alpha->GetParameter(1);
  Double_t sigma_alpha = std::abs(fit_alpha->GetParameter(2));
  Double_t mean_gamma = fit_gamma->GetParameter(1);
  Double_t sigma_gamma = std::abs(fit_gamma->GetParameter(2));

  out_amp_alpha = fit_alpha->GetParameter(0);
  out_mean_alpha = mean_alpha;
  out_sigma_alpha = sigma_alpha;
  out_amp_gamma = fit_gamma->GetParameter(0);
  out_mean_gamma = mean_gamma;
  out_sigma_gamma = sigma_gamma;

  delete fit_alpha;
  delete fit_gamma;

  if (sigma_alpha <= 0 || sigma_gamma <= 0)
    return -1;

  fwhm_alpha = 2.355 * sigma_alpha;
  fwhm_gamma = 2.355 * sigma_gamma;
  separation = std::abs(mean_alpha - mean_gamma);

  Double_t fom = separation / (fwhm_alpha + fwhm_gamma);
  return fom;
}

std::vector<OptimizationResult>
OptimizeGates(const std::vector<Event> &alpha_events,
              const std::vector<Event> &gamma_events) {

  std::cout << "Starting gate optimization..." << std::endl;
  std::cout << "  Light output range: "
            << OPTIMIZATION_PARAMETERS.min_light_output << "-"
            << OPTIMIZATION_PARAMETERS.max_light_output << " keVee"
            << std::endl;
  std::cout << "  Short gate range: " << OPTIMIZATION_PARAMETERS.min_short_gate
            << "-" << OPTIMIZATION_PARAMETERS.max_short_gate << " (step "
            << OPTIMIZATION_PARAMETERS.short_gate_step << ")" << std::endl;
  std::cout << "  Long gate range: " << OPTIMIZATION_PARAMETERS.min_long_gate
            << "-" << OPTIMIZATION_PARAMETERS.max_long_gate << " (step "
            << OPTIMIZATION_PARAMETERS.long_gate_step << ")" << std::endl;

  std::vector<GateCombination> combinations;
  for (Int_t sg = OPTIMIZATION_PARAMETERS.min_short_gate;
       sg <= OPTIMIZATION_PARAMETERS.max_short_gate;
       sg += OPTIMIZATION_PARAMETERS.short_gate_step) {
    for (Int_t lg = OPTIMIZATION_PARAMETERS.min_long_gate;
         lg <= OPTIMIZATION_PARAMETERS.max_long_gate;
         lg += OPTIMIZATION_PARAMETERS.long_gate_step) {
      if (sg < lg) {
        combinations.push_back({sg, lg});
      }
    }
  }

  std::cout << "  Total valid combinations to test: " << combinations.size()
            << std::endl;

  Int_t pre_gate = Constants::DEFAULT_PROCESSING_CONFIG.pre_gate;

  unsigned int n_threads =
      std::min<unsigned int>(std::thread::hardware_concurrency(),
                             static_cast<unsigned int>(combinations.size()));
  if (n_threads == 0)
    n_threads = 1;
  std::cout << "  Using " << n_threads << " threads" << std::endl;

  std::mutex print_mutex;
  std::vector<std::future<ComputeResult>> futures;
  futures.reserve(combinations.size());

  for (size_t i = 0; i < combinations.size(); i++) {
    futures.push_back(std::async(
        std::launch::async,
        [&alpha_events, &gamma_events, &print_mutex,
         pre_gate](Int_t sg, Int_t lg) -> ComputeResult {
          {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "  Computing gates " << sg << "/" << lg << "..."
                      << std::endl;
          }

          ComputeResult result;
          result.short_gate = sg;
          result.long_gate = lg;
          result.alpha_cc = ComputeCCValues(alpha_events, sg, lg, pre_gate);
          result.gamma_cc = ComputeCCValues(gamma_events, sg, lg, pre_gate);

          return result;
        },
        combinations[i].short_gate, combinations[i].long_gate));
  }

  std::vector<OptimizationResult> results;

  for (size_t i = 0; i < futures.size(); i++) {
    ComputeResult cr = futures[i].get();

    if (static_cast<Int_t>(cr.alpha_cc.size()) < 100 ||
        static_cast<Int_t>(cr.gamma_cc.size()) < 100) {
      std::cout << "  Insufficient statistics for gates " << cr.short_gate
                << "/" << cr.long_gate << " (alpha:" << cr.alpha_cc.size()
                << ", gamma:" << cr.gamma_cc.size() << ")" << std::endl;
      continue;
    }

    TH1F *hist_alpha = BuildHistogram(cr.alpha_cc);
    TH1F *hist_gamma = BuildHistogram(cr.gamma_cc);

    Double_t separation, fwhm_alpha, fwhm_gamma;
    Double_t amp_alpha, mean_alpha, sigma_alpha;
    Double_t amp_gamma, mean_gamma, sigma_gamma;
    Double_t fom = CalculateFigureOfMerit(
        hist_alpha, hist_gamma, separation, fwhm_alpha, fwhm_gamma, amp_alpha,
        mean_alpha, sigma_alpha, amp_gamma, mean_gamma, sigma_gamma);

    if (fom > 0) {
      OptimizationResult opt;
      opt.short_gate = cr.short_gate;
      opt.long_gate = cr.long_gate;
      opt.fom = fom;
      opt.amp_alpha = amp_alpha;
      opt.mean_alpha = mean_alpha;
      opt.sigma_alpha = sigma_alpha;
      opt.amp_gamma = amp_gamma;
      opt.mean_gamma = mean_gamma;
      opt.sigma_gamma = sigma_gamma;
      results.push_back(opt);

      std::cout << "  Gates " << cr.short_gate << "/" << cr.long_gate
                << " -> FOM = " << fom << std::endl;
    } else {
      std::cout << "  Gates " << cr.short_gate << "/" << cr.long_gate
                << " -> FOM calculation failed" << std::endl;
    }

    delete hist_alpha;
    delete hist_gamma;
  }

  std::sort(results.begin(), results.end(),
            [](const OptimizationResult &a, const OptimizationResult &b) {
              return a.fom > b.fom;
            });

  std::cout << "Optimization completed. Found " << results.size()
            << " valid combinations." << std::endl;

  return results;
}

void PlotPSDWithFOM(TH1F *hist_alpha, TH1F *hist_gamma,
                    const TString alpha_label, const TString gamma_label,
                    const OptimizationResult &opt, const TString psd_type,
                    const TString output_name) {
  if (!hist_alpha || !hist_gamma) {
    std::cout << "Error: Invalid histograms for PSD+FOM plot" << std::endl;
    return;
  }

  if (hist_alpha->GetEntries() == 0 || hist_gamma->GetEntries() == 0) {
    std::cout << "Error: Empty histograms for PSD+FOM plot" << std::endl;
    return;
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);

  PlottingUtils::ConfigureHistogram(hist_alpha, kRed + 2, "");
  PlottingUtils::ConfigureHistogram(hist_gamma, kBlue + 2, "");

  hist_alpha->SetFillColorAlpha(kRed + 2, 0.3);
  hist_gamma->SetFillColorAlpha(kBlue + 2, 0.3);

  hist_alpha->GetXaxis()->SetTitle(psd_type);
  hist_alpha->GetYaxis()->SetTitle("Counts");

  Double_t max_alpha = hist_alpha->GetMaximum();
  Double_t max_gamma = hist_gamma->GetMaximum();
  Double_t y_max = std::max(max_alpha, max_gamma) * 1.2;

  hist_alpha->SetMaximum(y_max);
  hist_alpha->Draw("HIST");
  hist_gamma->Draw("HIST SAME");

  // Reconstruct fit curves from stored parameters (no refitting)
  Double_t draw_min_alpha = opt.mean_alpha - 3.0 * opt.sigma_alpha;
  Double_t draw_max_alpha = opt.mean_alpha + 3.0 * opt.sigma_alpha;
  TF1 *fit_alpha = new TF1(PlottingUtils::GetRandomName(), "gaus",
                           draw_min_alpha, draw_max_alpha);
  fit_alpha->SetParameters(opt.amp_alpha, opt.mean_alpha, opt.sigma_alpha);
  fit_alpha->SetLineColor(kRed);
  fit_alpha->SetLineWidth(2);
  fit_alpha->SetLineStyle(1);
  fit_alpha->Draw("SAME");

  Double_t draw_min_gamma = opt.mean_gamma - 3.0 * opt.sigma_gamma;
  Double_t draw_max_gamma = opt.mean_gamma + 3.0 * opt.sigma_gamma;
  TF1 *fit_gamma = new TF1(PlottingUtils::GetRandomName(), "gaus",
                           draw_min_gamma, draw_max_gamma);
  fit_gamma->SetParameters(opt.amp_gamma, opt.mean_gamma, opt.sigma_gamma);
  fit_gamma->SetLineColor(kBlue);
  fit_gamma->SetLineWidth(2);
  fit_gamma->SetLineStyle(1);
  fit_gamma->Draw("SAME");

  TLegend *leg = PlottingUtils::AddLegend(0.2, 0.45, 0.55, 0.83);
  leg->AddEntry(hist_alpha, Form("%s (#alpha)", alpha_label.Data()), "f");
  leg->AddEntry(hist_gamma, Form("%s (#gamma)", gamma_label.Data()), "f");
  leg->AddEntry(fit_alpha, "#alpha fit", "l");
  leg->AddEntry(fit_gamma, "#gamma fit", "l");

  char fom_text[100];
  sprintf(fom_text, "FOM: %.3f", opt.fom);
  leg->AddEntry((TObject *)0, fom_text, "");
  leg->Draw();

  PlottingUtils::SaveFigure(canvas, output_name + "_" + gamma_label + "", "",
                            PlotSaveOptions::kLOG);

  delete canvas;
}

void GateOptimization() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT, project_root + "/plots",
                                project_root + "/root_files");
  ROOT::EnableThreadSafety();

  TString alpha_output = Constants::AM241;
  TString gamma_output = Constants::NA22;
  TString alpha_label = Constants::AM241_LABEL;
  TString gamma_label = Constants::NA22_LABEL;

  std::vector<Event> alpha_events = LoadEvents(alpha_output);
  std::vector<Event> gamma_events = LoadEvents(gamma_output);

  std::vector<OptimizationResult> opt_results =
      OptimizeGates(alpha_events, gamma_events);
  std::cout << "Top gate combinations:" << std::endl;
  std::cout << "Rank | Short | Long | FOM" << std::endl;

  Int_t max_results = std::min(10, (Int_t)opt_results.size());
  for (Int_t i = 0; i < max_results; ++i) {
    OptimizationResult result = opt_results.at(i);
    printf("%4d | %5d | %4d | %6.3f \n", i + 1, result.short_gate,
           result.long_gate, result.fom);
  }

  if (opt_results.empty()) {
    std::cout << "No valid results found." << std::endl;
    return;
  }

  std::cout << "Creating PSD plot with optimal gates..." << std::endl;
  OptimizationResult best = opt_results.at(0);

  Int_t pre_gate = Constants::DEFAULT_PROCESSING_CONFIG.pre_gate;
  std::vector<Float_t> best_alpha_cc =
      ComputeCCValues(alpha_events, best.short_gate, best.long_gate, pre_gate);
  std::vector<Float_t> best_gamma_cc =
      ComputeCCValues(gamma_events, best.short_gate, best.long_gate, pre_gate);

  TH1F *best_alpha = BuildHistogram(best_alpha_cc);
  TH1F *best_gamma = BuildHistogram(best_gamma_cc);

  if (best_alpha && best_gamma) {
    std::string output_name = "optimal_gates_" +
                              std::to_string(best.short_gate) + "_" +
                              std::to_string(best.long_gate);
    PlotPSDWithFOM(best_alpha, best_gamma, alpha_label, gamma_label, best,
                   "PSP_{CC}", output_name);
  }

  delete best_alpha;
  delete best_gamma;
}
