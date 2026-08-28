#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iomanip>
#include <vector>

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_79 = 79.6142;
const Float_t E_BA133_81 = 80.9979;
const Float_t E_BA133_276 = 276.3989;
const Float_t E_BA133_303 = 302.8508;
const Float_t E_BA133_356 = 356.0129;
const Float_t E_BA133_384 = 383.8485;
const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;
const Float_t E_GE_73M = 68.752;
const Float_t E_CD114M = 95.9023;
const Float_t E_ANNIHILATION = 510.999;

const Float_t CAL_RANGE_LOW = -10;
const Float_t CAL_RANGE_HIGH = 600;

struct CalibrationData {
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;
};

struct PbGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid;
};

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label);

std::vector<Double_t> LoadEvents(const TString input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    return events;
  }

  TTree *tree = static_cast<TTree *>(file->Get("total_energy_tree"));
  if (tree)
    events = RooFitUtils::LoadEventsFromTree(tree, "energykeV");
  else
    std::cerr << "ERROR: total_energy_tree not found in filtered/" << input_name
              << ".root" << std::endl;

  file->Close();
  delete file;
  return events;
}

FitResult FitCalibrationPeak(const TString input_name, const TString peak_name,
                             const Bool_t interactive) {
  std::vector<Double_t> events = LoadEvents(input_name);
  if (events.empty())
    return {};

  Bool_t use_flat_background = kFALSE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  Float_t fit_low = 0, fit_high = 0;

  if (peak_name == "Am_59.5keV") {
    if (input_name == Constants::POSTREACTOR_AM241_20260113) {
      fit_low = 50;
      fit_high = 70;
    } else if (input_name == Constants::POSTREACTOR_AM241_BA133_20260116) {
      fit_low = 55;
      fit_high = 70;
    } else {
      fit_low = 51;
      fit_high = 71;
    }
  } else if (peak_name == "Ba_53.16keV") {
    fit_low = 48;
    fit_high = 58;
  } else if (peak_name == "Ba_80.98keV") {
    fit_low = 75;
    fit_high = 90;
  } else if (peak_name == "Cd114m_95.9keV") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115) {
      fit_low = 91;
      fit_high = 100;
    } else {
      fit_low = 91;
      fit_high = 103;
    }
  } else {
    return {};
  }

  RooFitUtils fitter(events, fit_low, fit_high, Constants::BIN_WIDTH_KEV,
                     use_flat_background, use_step, use_low_exp_tail,
                     use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter.SetInteractive();
  FitResult result;
  if (peak_name == "Ba_80.98keV")
    result =
        fitter.FitDoublePeak(input_name, peak_name, E_BA133_79, E_BA133_81);
  else
    result = fitter.FitSinglePeak(input_name, peak_name);
  return result;
}

FitResult FitCalibrationPeakFull(const TString input_name,
                                 const TString peak_name,
                                 const Bool_t interactive) {
  std::vector<Double_t> events = LoadEvents(input_name);
  if (events.empty())
    return {};

  Bool_t use_flat_background = kFALSE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  Float_t fit_low = 0, fit_high = 0;
  if (peak_name == "Ba_276.40keV") {
    fit_low = 270;
    fit_high = 283;
  } else if (peak_name == "Ba_302.85keV") {
    fit_low = 296;
    fit_high = 310;
  } else if (peak_name == "Ba_356.01keV") {
    fit_low = 349;
    fit_high = 363;
  } else if (peak_name == "Ba_383.85keV") {
    fit_low = 377;
    fit_high = 391;
  } else if (peak_name == "Annihilation_511keV") {
    fit_low = 500;
    fit_high = 525;
  } else {
    return {};
  }

  RooFitUtils fitter(events, fit_low, fit_high, Constants::BIN_WIDTH_KEV,
                     use_flat_background, use_step, use_low_exp_tail,
                     use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter.SetInteractive();
  return fitter.FitSinglePeak(input_name, peak_name);
}

PbGeSimResult FitPbGeSimultaneous(const TString bkg_input,
                                  const TString sig_input, const Float_t bkg_lo,
                                  const Float_t bkg_hi, const Float_t sig_lo,
                                  const Float_t sig_hi, const Bool_t bkg_flat,
                                  const Bool_t sig_flat,
                                  const Bool_t interactive) {
  PbGeSimResult out;
  out.valid = kFALSE;

  std::vector<Double_t> bkg_events = LoadEvents(bkg_input);
  std::vector<Double_t> sig_events = LoadEvents(sig_input);
  if (bkg_events.empty() || sig_events.empty())
    return out;

  Bool_t use_step = kFALSE;
  Bool_t use_low_exp = kTRUE;
  Bool_t use_low_lin = kTRUE;
  Bool_t use_high_exp = kTRUE;

  RooFitUtils bkg_fitter(bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV,
                         bkg_flat, use_step, use_low_exp, use_low_lin,
                         use_high_exp);
  if (interactive)
    bkg_fitter.SetInteractive();
  FitResult seed =
      bkg_fitter.FitDoublePeak(bkg_input, "PbKa_seed", E_PB_KA2, E_PB_KA1);
  if (!seed.valid) {
    std::cerr << "ERROR: bkg-only seed fit failed for " << bkg_input
              << std::endl;
    return out;
  }

  RooFitUtils sim;
  if (interactive)
    sim.SetInteractive();
  std::vector<Double_t> bkg_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  std::vector<Double_t> sig_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1,
                                   (Double_t)E_GE_73M};
  sim.AddChannel("bkg", bkg_events, bkg_lo, bkg_hi, Constants::BIN_WIDTH_KEV, 2,
                 bkg_mus, bkg_flat, use_step, use_low_exp, use_low_lin,
                 use_high_exp);
  sim.AddChannel("sig", sig_events, sig_lo, sig_hi, Constants::BIN_WIDTH_KEV, 3,
                 sig_mus, sig_flat, use_step, use_low_exp, use_low_lin,
                 use_high_exp);
  sim.LinkPeakShape("sig", 0, "bkg", 0);
  sim.LinkPeakShape("sig", 1, "bkg", 1);
  sim.SeedChannel("bkg", seed);

  std::vector<FitResult> results =
      sim.FitSimultaneous(sig_input, "PbKa_Ge_sim");

  if (results.size() < 2)
    return out;

  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;
  return out;
}

TF1 *CreateAndSaveCalibration(
    const std::vector<Float_t> &mu,
    const std::vector<Float_t> &calibration_values_keV,
    const std::vector<Float_t> &mu_errors, const TString date_label) {

  Int_t size = calibration_values_keV.size();

  TGraph *calibration_curve =
      new TGraph(size, mu.data(), calibration_values_keV.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      calibration_curve, kBlue,
      "; Precalibrated Energy [keV]; Deposited Energy [keV]");

  calibration_curve->GetXaxis()->SetRangeUser(-5, CAL_RANGE_HIGH);
  calibration_curve->GetYaxis()->SetRangeUser(-5, CAL_RANGE_HIGH);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(5);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("AP");

  TF1 *calibration_fit =
      new TF1("cal_" + date_label, "pol2", CAL_RANGE_LOW, CAL_RANGE_HIGH);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetParameter(2, 0);
  calibration_fit->SetNpx(1000);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");

  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return calibration_fit;
}

void AddZeroPoint(CalibrationData &cal_data) {
  cal_data.run_names.push_back("Zero Point");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);
}

void AddCalibrationPoint(CalibrationData &cal_data, const TString run_name,
                         const Float_t mu, const Float_t mu_error,
                         const Float_t true_energy,
                         const Float_t reduced_chi2) {
  cal_data.run_names.push_back(run_name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_error);
  cal_data.calibration_values_keV.push_back(true_energy);
  cal_data.reduced_chi2.push_back(reduced_chi2);
}

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label) {
  std::cout << "Calibration Points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); ++i) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0) {
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
    }
    std::cout << std::endl;
  }
}

Double_t InvertMasterPol2(TF1 *master, Double_t e_true) {
  Double_t p0 = master->GetParameter(0);
  Double_t p1 = master->GetParameter(1);
  Double_t p2 = master->GetParameter(2);

  if (std::abs(p2) > 1e-12) {
    Double_t disc = p1 * p1 - 4.0 * p2 * (p0 - e_true);
    if (disc >= 0)
      return (-p1 + std::sqrt(disc)) / (2.0 * p2);
    return (e_true - p0) / p1;
  }
  return (p1 != 0) ? (e_true - p0) / p1 : e_true;
}

TF1 *ComposeDriftCalibration(TF1 *master, const TString date_label,
                             const CalibrationData &cal_data) {
  Int_t n = cal_data.mu.size();
  if (n < 2) {
    std::cerr << "WARNING: Only " << n << " peaks for drift on " << date_label
              << ", need >= 2" << std::endl;
    return nullptr;
  }

  std::vector<Float_t> master_gm(n);
  for (Int_t i = 0; i < n; i++)
    master_gm[i] = static_cast<Float_t>(
        InvertMasterPol2(master, cal_data.calibration_values_keV[i]));

  TGraph *drift_graph = new TGraph(n, cal_data.mu.data(), master_gm.data());
  TF1 *drift_fit =
      new TF1("drift_" + date_label, "pol1", CAL_RANGE_LOW, CAL_RANGE_HIGH);
  drift_fit->SetParameter(0, 0);
  drift_fit->SetParameter(1, 1);
  drift_graph->Fit(drift_fit, "RQ");

  Double_t alpha = drift_fit->GetParameter(0);
  Double_t beta = drift_fit->GetParameter(1);

  std::cout << "  Drift " << date_label << ": alpha = " << std::fixed
            << std::setprecision(6) << alpha << ", beta = " << beta << " (" << n
            << " peaks)" << std::endl;

  for (Int_t i = 0; i < n; i++) {
    Double_t mapped = alpha + beta * cal_data.mu[i];
    Double_t cal_true = master->Eval(mapped);
    std::cout << "    " << std::left << std::setw(35) << cal_data.run_names[i]
              << " E_meas = " << std::fixed << std::setprecision(4)
              << cal_data.mu[i]
              << ", E_true = " << cal_data.calibration_values_keV[i]
              << ", cal = " << cal_true << ", residual = "
              << (cal_true - cal_data.calibration_values_keV[i]) << " keV"
              << std::endl;
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      drift_graph, kBlue,
      "; Measured Energy [keV]; Master Gain-Matched Energy [keV]");
  drift_graph->SetMarkerStyle(5);
  drift_graph->SetMarkerSize(2);
  drift_graph->Draw("AP");
  drift_fit->Draw("SAME");
  PlottingUtils::SaveFigure(canvas, "drift_" + date_label, "",
                            PlotSaveOptions::kLINEAR);

  Double_t p0 = master->GetParameter(0);
  Double_t p1 = master->GetParameter(1);
  Double_t p2 = master->GetParameter(2);

  TString formula = Form("[0] + [1]*(%.12g + %.12g*x) + "
                         "[2]*(%.12g + %.12g*x)*(%.12g + %.12g*x)",
                         alpha, beta, alpha, beta, alpha, beta);
  TF1 *cal_func =
      new TF1("cal_" + date_label, formula, CAL_RANGE_LOW, CAL_RANGE_HIGH);
  cal_func->SetParameter(0, p0);
  cal_func->SetParameter(1, p1);
  cal_func->SetParameter(2, p2);
  cal_func->SetNpx(1000);

  delete drift_graph;
  delete drift_fit;

  return cal_func;
}

TF1 *BuildCalibration_20260115(const Bool_t interactive) {
  std::cout << "Building Master Calibration (Jan 15)" << std::endl;

  CalibrationData cal_data;

  FitResult am241_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_20260115, "Am_59.5keV", interactive);
  if (am241_result.valid)
    AddCalibrationPoint(
        cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
        am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba53_result = FitCalibrationPeak(
      Constants::POSTREACTOR_BA133_20260115, "Ba_53.16keV", interactive);
  if (ba53_result.valid)
    AddCalibrationPoint(
        cal_data, "Post-reactor Ba-133 53.16keV", ba53_result.peaks.at(0).mu,
        ba53_result.peaks.at(0).mu_error, E_BA133_53, ba53_result.reduced_chi2);

  FitResult ba81_result = FitCalibrationPeak(
      Constants::POSTREACTOR_BA133_20260115, "Ba_80.98keV", interactive);
  if (ba81_result.valid) {
    AddCalibrationPoint(
        cal_data, "Post-reactor Ba-133 79.61keV", ba81_result.peaks.at(0).mu,
        ba81_result.peaks.at(0).mu_error, E_BA133_79, ba81_result.reduced_chi2);
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 80.98keV",
                        ba81_result.peaks.at(1).mu,
                        ba81_result.peaks.at(1).mu_error, E_BA133_81, -1);
  }

  FitResult ba276_result = FitCalibrationPeakFull(
      Constants::POSTREACTOR_BA133_20260115, "Ba_276.40keV", interactive);
  if (ba276_result.valid)
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 276.40keV",
                        ba276_result.peaks.at(0).mu,
                        ba276_result.peaks.at(0).mu_error, E_BA133_276,
                        ba276_result.reduced_chi2);

  FitResult ba303_result = FitCalibrationPeakFull(
      Constants::POSTREACTOR_BA133_20260115, "Ba_302.85keV", interactive);
  if (ba303_result.valid)
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 302.85keV",
                        ba303_result.peaks.at(0).mu,
                        ba303_result.peaks.at(0).mu_error, E_BA133_303,
                        ba303_result.reduced_chi2);

  FitResult ba356_result = FitCalibrationPeakFull(
      Constants::POSTREACTOR_BA133_20260115, "Ba_356.01keV", interactive);
  if (ba356_result.valid)
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 356.01keV",
                        ba356_result.peaks.at(0).mu,
                        ba356_result.peaks.at(0).mu_error, E_BA133_356,
                        ba356_result.reduced_chi2);

  FitResult ba384_result = FitCalibrationPeakFull(
      Constants::POSTREACTOR_BA133_20260115, "Ba_383.85keV", interactive);
  if (ba384_result.valid)
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 383.85keV",
                        ba384_result.peaks.at(0).mu,
                        ba384_result.peaks.at(0).mu_error, E_BA133_384,
                        ba384_result.reduced_chi2);

  FitResult ann_result =
      FitCalibrationPeakFull(Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                             "Annihilation_511keV", interactive);
  if (ann_result.valid)
    AddCalibrationPoint(cal_data, "511 keV Annihilation",
                        ann_result.peaks.at(0).mu,
                        ann_result.peaks.at(0).mu_error, E_ANNIHILATION,
                        ann_result.reduced_chi2);

  PrintCalibrationSummary(cal_data, "Master_20260115");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "Master_20260115");
  return cal_func;
}

TF1 *BuildCalibration_20260113(TF1 *master, const Bool_t interactive) {
  std::cout << "Building Drift Calibration for Jan 13" << std::endl;

  CalibrationData cal_data;

  FitResult am241_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_20260113, "Am_59.5keV", interactive);
  if (am241_result.valid)
    AddCalibrationPoint(
        cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
        am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  PbGeSimResult cu_pair =
      FitPbGeSimultaneous(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                          Constants::CUSHIELDSIGNAL_10PERCENT_20260113, 65, 82,
                          62, 80, kFALSE, kTRUE, interactive);
  if (cu_pair.valid) {
    AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                        cu_pair.bkg_channel.peaks.at(0).mu,
                        cu_pair.bkg_channel.peaks.at(0).mu_error, E_PB_KA2,
                        cu_pair.bkg_channel.reduced_chi2);
    AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                        cu_pair.bkg_channel.peaks.at(1).mu,
                        cu_pair.bkg_channel.peaks.at(1).mu_error, E_PB_KA1, -1);
  }

  FitResult ann_result =
      FitCalibrationPeakFull(Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                             "Annihilation_511keV", interactive);
  if (ann_result.valid)
    AddCalibrationPoint(cal_data, "511 keV Annihilation",
                        ann_result.peaks.at(0).mu,
                        ann_result.peaks.at(0).mu_error, E_ANNIHILATION,
                        ann_result.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260113");

  return ComposeDriftCalibration(master, "20260113", cal_data);
}

TF1 *BuildCalibration_20260114(TF1 *master, const Bool_t interactive) {
  std::cout << "Building Drift Calibration for Jan 14" << std::endl;

  CalibrationData cal_data;

  PbGeSimResult cu_pair =
      FitPbGeSimultaneous(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                          Constants::CUSHIELDSIGNAL_10PERCENT_20260114, 66, 82,
                          63, 80, kTRUE, kTRUE, interactive);
  if (cu_pair.valid) {
    AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                        cu_pair.bkg_channel.peaks.at(0).mu,
                        cu_pair.bkg_channel.peaks.at(0).mu_error, E_PB_KA2,
                        cu_pair.bkg_channel.reduced_chi2);
    AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                        cu_pair.bkg_channel.peaks.at(1).mu,
                        cu_pair.bkg_channel.peaks.at(1).mu_error, E_PB_KA1, -1);
  }

  FitResult ann_result =
      FitCalibrationPeakFull(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                             "Annihilation_511keV", interactive);
  if (ann_result.valid)
    AddCalibrationPoint(cal_data, "511 keV Annihilation",
                        ann_result.peaks.at(0).mu,
                        ann_result.peaks.at(0).mu_error, E_ANNIHILATION,
                        ann_result.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260114");

  return ComposeDriftCalibration(master, "20260114", cal_data);
}

TF1 *BuildCalibration_20260116(TF1 *master, const Bool_t interactive) {
  std::cout << "Building Drift Calibration for Jan 16" << std::endl;

  CalibrationData cal_data;

  FitResult am241_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_BA133_20260116, "Am_59.5keV", interactive);
  if (am241_result.valid)
    AddCalibrationPoint(
        cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
        am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba53_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_BA133_20260116, "Ba_53.16keV", interactive);
  if (ba53_result.valid)
    AddCalibrationPoint(
        cal_data, "Post-reactor Ba-133 53.16keV", ba53_result.peaks.at(0).mu,
        ba53_result.peaks.at(0).mu_error, E_BA133_53, ba53_result.reduced_chi2);

  FitResult ba81_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_BA133_20260116, "Ba_80.98keV", interactive);
  if (ba81_result.valid) {
    AddCalibrationPoint(
        cal_data, "Post-reactor Ba-133 79.61keV", ba81_result.peaks.at(0).mu,
        ba81_result.peaks.at(0).mu_error, E_BA133_79, ba81_result.reduced_chi2);
    AddCalibrationPoint(cal_data, "Post-reactor Ba-133 80.98keV",
                        ba81_result.peaks.at(1).mu,
                        ba81_result.peaks.at(1).mu_error, E_BA133_81, -1);
  }

  FitResult ann_result = FitCalibrationPeakFull(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      "Annihilation_511keV", interactive);
  if (ann_result.valid)
    AddCalibrationPoint(cal_data, "511 keV Annihilation",
                        ann_result.peaks.at(0).mu,
                        ann_result.peaks.at(0).mu_error, E_ANNIHILATION,
                        ann_result.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260116");

  return ComposeDriftCalibration(master, "20260116", cal_data);
}

void PulseHeightToDepositedEnergy(const std::vector<TString> &input_names,
                                  TF1 *calibration_function,
                                  const TString date_label) {
  TFile *cal_file = IO::OpenForWriting("calibrated/calibration_function_" +
                                       date_label + ".root");
  calibration_function->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  Int_t n_files = input_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_name = input_names[i];
    TFile *in_file = IO::OpenForReading("filtered/" + input_name + ".root");
    TFile *file = IO::OpenForWriting("calibrated/" + input_name + ".root");

    Float_t deposited_energy = 0;
    TTree *cal_tree = new TTree("calibrated_tree", "Calibrated events");
    cal_tree->SetDirectory(file);
    cal_tree->Branch("depositedEnergykeV", &deposited_energy,
                     "depositedEnergykeV/F");

    TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                          Form("; Deposited Energy [keV]; Counts / %d eV",
                               Constants::BIN_WIDTH_EV),
                          Constants::HIST_NBINS, Constants::HIST_XMIN,
                          Constants::HIST_XMAX);
    hist->SetDirectory(0);

    TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                                Form("; Deposited Energy [keV]; Counts / %d eV",
                                     Constants::BIN_WIDTH_EV),
                                Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                                Constants::ZOOMED_XMAX);
    zoomedHist->SetDirectory(0);

    TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                              Form("; Deposited Energy [keV]; Counts / %d eV",
                                   Constants::BIN_WIDTH_EV),
                              Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                              Constants::PEAK_XMAX);
    peakHist->SetDirectory(0);

    TTree *tree = static_cast<TTree *>(in_file->Get("total_energy_tree"));
    if (tree) {
      Float_t energy = 0;
      tree->SetBranchAddress("energykeV", &energy);
      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        deposited_energy = calibration_function->Eval(energy);
        cal_tree->Fill();
        hist->Fill(deposited_energy);
        zoomedHist->Fill(deposited_energy);
        peakHist->Fill(deposited_energy);
      }
    } else {
      std::cerr << "ERROR: total_energy_tree not found in filtered/"
                << input_name << ".root" << std::endl;
    }

    file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);

    delete hist;
    delete zoomedHist;
    delete peakHist;

    std::cout << "Calibrated " << input_name << std::endl;
    file->Close();
    delete file;
    in_file->Close();
    delete in_file;
  }
}

void Calibration() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/raw",
                                project_root + "/root_files");

  Bool_t interactive = kTRUE;

  TF1 *master = BuildCalibration_20260115(interactive);
  TF1 *cal_jan13 = BuildCalibration_20260113(master, interactive);
  TF1 *cal_jan14 = BuildCalibration_20260114(master, interactive);
  TF1 *cal_jan15 = master;
  TF1 *cal_jan16 = BuildCalibration_20260116(master, interactive);

  std::cout << std::endl;
  std::cout << "Master Calibration (Jan 15, pol2):" << std::endl;
  std::cout << "  p0 = " << std::fixed << std::setprecision(6)
            << master->GetParameter(0) << " +/- " << master->GetParError(0)
            << std::endl;
  std::cout << "  p1 = " << master->GetParameter(1) << " +/- "
            << master->GetParError(1) << std::endl;
  std::cout << "  p2 = " << std::scientific << std::setprecision(4)
            << master->GetParameter(2) << " +/- " << master->GetParError(2)
            << std::endl;
  std::cout << std::endl;

  std::vector<TString> datasets_jan13 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
      Constants::POSTREACTOR_AM241_20260113,
      Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113,
      Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113};
  PulseHeightToDepositedEnergy(datasets_jan13, cal_jan13, "20260113");

  std::vector<TString> datasets_jan14 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_90PERCENT_20260114};
  PulseHeightToDepositedEnergy(datasets_jan14, cal_jan14, "20260114");

  std::vector<TString> datasets_jan15 = {
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
      Constants::POSTREACTOR_AM241_20260115,
      Constants::POSTREACTOR_BA133_20260115, Constants::SHUTTERCLOSED_20260115};
  PulseHeightToDepositedEnergy(datasets_jan15, cal_jan15, "Master_20260115");

  std::vector<TString> datasets_jan16 = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  PulseHeightToDepositedEnergy(datasets_jan16, cal_jan16, "20260116");
}
