#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "GainTransferFit.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <RtypesCore.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>

// calibration constants and the shared fit core live in GainTransferFit.hpp

struct CalibrationData {
  std::vector<Float_t> mu, mu_errors, calibration_values_keV, reduced_chi2;
  std::vector<Float_t> energy_errors;
  std::vector<TString> run_names;
};

struct DayResult {
  std::vector<TString> cal_labels;
  std::vector<TF1 *> cal_funcs;
  std::vector<TString> ge_labels;
  std::vector<Float_t> ge_mus, ge_errs, ge_chi2s;
  std::vector<Float_t> ge_cal_errs;
  // The OFFSET part of the reference calibration error: its value at the
  // best-constrained energy. The remainder, sqrt(ref^2 - off^2), is the
  // multiplicative term, and the gain transfer is a pure scale factor.
  std::vector<Float_t> ge_cal_offs;
  // The TRANSFER part of ge_cal_errs, parallel to ge_*. ge_cal_err is
  // hypot(reference_cal_err, gain_err), so the reference part -- COMMON to
  // every run on the day -- is sqrt(cal^2 - gain^2) and the transfer part is
  // independent per run. The combiner needs that split to build the covariance
  // correctly instead of assuming the whole calibration is correlated or none
  // of it is.
  std::vector<Float_t> ge_gain_errs;
  // With the signal-run Pb anchor: the error of that run's own offset, also
  // inside ge_cal_errs. It comes from the run's own Pb fit, so it is
  // independent between runs, unlike the transfer.
  std::vector<Float_t> ge_anchor_errs;
  std::vector<Float_t> ge_bkg_errs;
  std::vector<TString> rs_labels;
  std::vector<Float_t> rs_mus, rs_errs, rs_chi2s, rs_cal_errs, rs_gain_errs;
  std::vector<Float_t> rs_cal_offs, rs_sys_errs;
};

struct CalFit {
  TF1 *func = nullptr;
  Double_t p0 = 0, p1 = 1;
  Int_t npar = 2;
  Double_t var_p0 = 0, var_p1 = 0, cov_p0p1 = 0;
};

struct PairCalResult {
  TF1 *cal_func = nullptr;
  Double_t ge_cal_err = 0; // calibration sigma propagated to the Ge energy
  Double_t ge_cal_off = 0; // its value at the best-constrained energy

  Double_t ge_bkg_err = 0; // |cal(precal Ge) - postcal Ge|, a diagnostic only
  BkgGeSimResult postcal;
  TString pair_tag;
  std::vector<Float_t> pb_mus_precal;
  std::vector<Float_t> pb_mu_errs_precal;
  CalFit ref_calfit;
  Bool_t has_ref_cal = kFALSE;
  Double_t ge_gain_err = 0;
  Double_t ge_anchor_err = 0;
};

struct LineCalConfig {
  TString date_label;
  TString am_run;
  Float_t am_lo, am_hi;
  TString ba_run;
  std::vector<TString> day_datasets;
  TString postcal_bkg;
  TString postcal_sig;
  TString ge_label;
};

std::vector<Double_t> LoadInteriorEnergies(TTree *tree) {
  std::vector<Double_t> energies;
  Float_t energy = 0, x = 0, y = 0;
  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  Long64_t n = tree->GetEntries();
  for (Long64_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    if (InInteriorPixels(x, y))
      energies.push_back(energy);
  }
  tree->ResetBranchAddresses();
  return energies;
}

std::vector<Double_t> LoadFiltered(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    return events;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    std::vector<Double_t> ch =
        InteriorPixelsRequested()
            ? LoadInteriorEnergies(tree)
            : RooFitUtils::LoadEventsFromTree(tree, "energykeV");
    events.insert(events.end(), ch.begin(), ch.end());
  }
  file->Close();
  delete file;
  return events;
}

void PrintAnchorFitStatus(const TString &peak_name, const TString &input_name,
                          const FitResult &fr) {
  std::cout << "ANCHOR-FIT " << peak_name << " " << input_name << ": status "
            << fr.fit_status << "  covQual " << fr.cov_qual << "  retried "
            << (fr.fit_retried ? "yes" : "no") << "  chi2/ndf " << std::fixed
            << std::setprecision(3) << fr.reduced_chi2;
  if (fr.peaks.size() == 2 && fr.peaks[1].gaus_amplitude > 0)
    std::cout << "  yield ratio " << std::setprecision(4)
              << fr.peaks[0].gaus_amplitude / fr.peaks[1].gaus_amplitude
              << "  tau/sigma " << std::setprecision(2)
              << fr.peaks[1].low_exp_tail_ratio;
  std::cout << std::endl;
  if (fr.fit_status >= RooFitUtils::kFirstInvalidMinimumStatus)
    std::cerr << "WARNING: " << peak_name << " " << input_name
              << " ended without a valid minimum" << std::endl;
}

FitResult FitCalPeak(const std::vector<Double_t> &events,
                     const TString &input_name, const TString &peak_name,
                     Float_t fit_low, Float_t fit_high, Bool_t double_peak,
                     Double_t mu1, Double_t mu2, Bool_t link_sigma = kFALSE,
                     Bool_t link_tails = kFALSE,
                     Double_t gaus_yield_ratio = 0.0) {
  if (events.empty())
    return {};
  RooFitUtils fitter(events, fit_low, fit_high, Constants::BIN_WIDTH_KEV, kTRUE,
                     kFALSE, kTRUE, kFALSE, kFALSE);
  fitter.SetTailRatioMax(TAIL_RATIO_MAX);
  fitter.SetInteractive();
  FitResult fr;
  if (double_peak)
    fr = fitter.FitDoublePeak(input_name, peak_name, mu1, mu2, link_sigma,
                              link_tails, gaus_yield_ratio);
  else
    fr = fitter.FitSinglePeak(input_name, peak_name);
  PrintAnchorFitStatus(peak_name, input_name, fr);
  return fr;
}

void AddCalPoint(CalibrationData &cal_data, const TString &name, Float_t mu,
                 Float_t mu_err, Float_t true_e, Float_t chi2,
                 Float_t energy_err = 0) {
  cal_data.run_names.push_back(name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_err);
  cal_data.calibration_values_keV.push_back(true_e);
  cal_data.reduced_chi2.push_back(chi2);
  cal_data.energy_errors.push_back(energy_err);
}

Float_t GuardedMuError(const TString &name, const PeakFitResult &pk,
                       const FitResult &fr) {
  if (!(pk.sigma > 0) || !(pk.gaus_amplitude > 0))
    return pk.mu_error;
  const Float_t floor_err = pk.sigma / TMath::Sqrt(pk.gaus_amplitude);
  if (pk.mu_error > floor_err)
    return pk.mu_error;
  std::cerr << "WARNING: " << name << " centroid error " << pk.mu_error * 1000
            << " eV is below the Poisson floor " << floor_err * 1000
            << " eV (sigma " << pk.sigma << ", N " << pk.gaus_amplitude
            << ", covQual " << fr.cov_qual << ") -- raising to the floor. "
            << "Re-tune this fit's saved state under TAIL_RATIO_MAX."
            << std::endl;
  return floor_err;
}

void PrintCalSummary(const CalibrationData &cal_data,
                     const TString &date_label) {
  std::cout << "Calibration points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); i++) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0)
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
    std::cout << std::endl;
  }
}

Double_t CalEnergyError(const CalFit &c, Double_t mu) {
  Double_t var = c.var_p0 + mu * mu * c.var_p1 + 2.0 * mu * c.cov_p0p1;
  return (var > 0) ? std::sqrt(var) : 0.0;
}

// The calibration uncertainty at its best-constrained energy, the minimum of
// CalEnergyError over the anchor region. For a straight line this is the
// intercept error at the pivot; sqrt(cal^2 - off^2) is then the
// multiplicative term.
Double_t CalOffsetError(const CalFit &c) {
  Double_t best = 0;
  for (Double_t x = 50.0; x <= 100.0; x += 0.01) {
    Double_t e = CalEnergyError(c, x);
    if (e > 0 && (best <= 0 || e < best))
      best = e;
  }
  return best;
}

// Parameter covariance as (J^T W J)^-1 with the effective variances
// ey^2 + (f'(x) ex)^2 at the fitted parameters, exact for the linearised
// problem.
void FillCalCovariance(CalFit &c, const std::vector<Float_t> &x,
                       const std::vector<Float_t> &ex,
                       const std::vector<Float_t> &ey) {
  Int_t nf = c.npar;
  TMatrixD H(nf, nf);
  for (size_t i = 0; i < x.size(); i++) {
    Double_t xi = x[i];
    Double_t fp = c.func->Derivative(xi);
    Double_t s2 = (Double_t)ey[i] * ey[i] + fp * fp * (Double_t)ex[i] * ex[i];
    for (Int_t a = 0; a < nf; a++)
      for (Int_t b = 0; b < nf; b++)
        H(a, b) += std::pow(xi, a) * std::pow(xi, b) / s2;
  }
  H.Invert();
  c.var_p0 = H(0, 0);
  c.var_p1 = H(1, 1);
  c.cov_p0p1 = H(0, 1);
  for (Int_t k = 0; k < nf; k++)
    c.func->SetParError(k, std::sqrt(std::max(0.0, H(k, k))));
}

CalFit CreateAndSavePol1Cal(const CalibrationData &cal_data,
                            const TString &date_label) {
  Int_t n = (Int_t)cal_data.mu.size();
  std::vector<Float_t> ex(cal_data.mu_errors);
  std::vector<Float_t> ey(n);
  for (Int_t i = 0; i < n; i++) {
    Float_t e = (i < (Int_t)cal_data.energy_errors.size())
                    ? cal_data.energy_errors[i]
                    : 0.0f;
    ey[i] = (e > 0)
                ? e
                : 1e-4f; // tiny floor so the point isn't infinitely weighted
    if (!(ex[i] > 0))
      ex[i] = 1e-4f;
  }
  TGraphErrors *graph = new TGraphErrors(n, cal_data.mu.data(),
                                         cal_data.calibration_values_keV.data(),
                                         ex.data(), ey.data());
  Double_t draw_lo = CAL_RANGE_LOW, draw_hi = CAL_RANGE_HIGH;
  for (Int_t i = 0; i < n; i++) {
    draw_lo = std::min({draw_lo, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
    draw_hi = std::max({draw_hi, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
  }
  draw_lo -= 3.0;
  draw_hi += 3.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      graph, kBlue, "; Precalibrated Energy [keV]; Deposited Energy [keV]");
  graph->GetXaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetYaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetMarkerStyle(5);
  graph->SetMarkerSize(2);
  graph->Draw("AP");

  TF1 *cal = new TF1("cal_low_" + date_label, "pol1", draw_lo, draw_hi);
  cal->SetParameter(0, 0);
  cal->SetParameter(1, 1);
  cal->SetNpx(1000);
  TFitResultPtr fr = graph->Fit(cal, "E S Q");
  cal->Draw("SAME");

  CalFit out;
  out.func = cal;
  out.npar = cal->GetNpar();
  out.p0 = cal->GetParameter(0);
  out.p1 = cal->GetParameter(1);
  FillCalCovariance(out, cal_data.mu, ex, ey);
  if (fr.Get() && fr->Ndf() > 0)
    std::cout << "CAL-FIT " << date_label << " (pol" << (out.npar - 1) << ", "
              << n << " anchors): chi2/ndf = " << std::fixed
              << std::setprecision(2) << fr->Chi2() / fr->Ndf()
              << "   sigma_E at best-constrained point = "
              << std::setprecision(1) << CalOffsetError(out) * 1000.0 << " eV"
              << std::endl;

  PlottingUtils::SaveFigure(canvas, "calibration_low_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return out;
}

struct GainResult {
  Double_t g = 1.0; // error-weighted mean relative gain
  Double_t sigma_g = 0.0;
  Double_t chi2_ndf = 0.0;
  Int_t n = 0;
};

GainResult PbGainVsReferenceWithError(const std::vector<Float_t> &mu_run,
                                      const std::vector<Float_t> &err_run,
                                      const std::vector<Float_t> &mu_ref,
                                      const std::vector<Float_t> &err_ref) {
  GainResult out;
  Double_t sw = 0, swg = 0;
  std::vector<Double_t> gs, ws;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    Double_t g = (Double_t)mu_run[i] / (Double_t)mu_ref[i];
    Double_t fr = (err_run[i] > 0) ? err_run[i] / mu_run[i] : 0.0;
    Double_t ff = (err_ref[i] > 0) ? err_ref[i] / mu_ref[i] : 0.0;
    Double_t rel = std::sqrt(fr * fr + ff * ff);
    Double_t sg = (rel > 0) ? g * rel : 1e-9;
    Double_t w = 1.0 / (sg * sg);
    sw += w;
    swg += w * g;
    gs.push_back(g);
    ws.push_back(w);
  }
  if (sw <= 0)
    return out;
  out.g = swg / sw;
  out.sigma_g = std::sqrt(1.0 / sw);
  out.n = (Int_t)gs.size();
  if (out.n > 1) {
    Double_t chi2 = 0;
    for (size_t i = 0; i < gs.size(); i++)
      chi2 += ws[i] * (gs[i] - out.g) * (gs[i] - out.g);
    out.chi2_ndf = chi2 / (out.n - 1);
    if (out.chi2_ndf > 1.0)
      out.sigma_g *= std::sqrt(out.chi2_ndf);
  }
  return out;
}

CalFit GainScaleCal(const CalFit &ref, Double_t g, const TString &label) {
  TF1 *src = ref.func;
  Int_t npar = src->GetNpar();
  Double_t lo = 0, hi = 0;
  src->GetRange(lo, hi);
  TF1 *scaled = new TF1("cal_low_" + label, src->GetExpFormula(), lo, hi);
  Double_t gk = 1.0;
  for (Int_t k = 0; k < npar; k++) {
    scaled->SetParameter(k, src->GetParameter(k) / gk);
    gk *= g;
  }
  scaled->SetNpx(1000);
  CalFit out;
  out.func = scaled;
  out.p0 = scaled->GetParameter(0);
  out.p1 = scaled->GetParameter(1);
  out.npar = npar;
  out.var_p0 = ref.var_p0;
  out.var_p1 = ref.var_p1 / (g * g);
  out.cov_p0p1 = ref.cov_p0p1 / g;
  if (npar > 0)
    scaled->SetParError(0, std::sqrt(std::max(0.0, out.var_p0)));
  if (npar > 1)
    scaled->SetParError(1, std::sqrt(std::max(0.0, out.var_p1)));
  return out;
}

Double_t GainTransferEnergyError(const CalFit &ref_cal, Double_t mu_raw,
                                 Double_t g, Double_t sigma_g) {
  if (!(sigma_g > 0) || !(g > 0) || !ref_cal.func)
    return 0.0;
  Double_t fprime = ref_cal.func->Derivative(mu_raw / g);
  return std::fabs(-mu_raw / (g * g) * fprime) * sigma_g;
}

struct OffsetResult {
  Double_t delta = 0; // raw run - raw reference, keV
  Double_t sigma = 0;
  Double_t chi2_ndf = 0;
  Int_t n = 0;
};

OffsetResult PbOffsetVsReference(const std::vector<Float_t> &mu_run,
                                 const std::vector<Float_t> &err_run,
                                 const std::vector<Float_t> &mu_ref,
                                 const std::vector<Float_t> &err_ref) {
  OffsetResult out;
  Double_t sw = 0, swd = 0;
  std::vector<Double_t> ds, ws;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    Double_t d = (Double_t)mu_run[i] - (Double_t)mu_ref[i];
    Double_t var =
        (Double_t)err_run[i] * err_run[i] + (Double_t)err_ref[i] * err_ref[i];
    Double_t w = (var > 0) ? 1.0 / var : 1e18;
    sw += w;
    swd += w * d;
    ds.push_back(d);
    ws.push_back(w);
  }
  if (sw <= 0)
    return out;
  out.delta = swd / sw;
  out.sigma = std::sqrt(1.0 / sw);
  out.n = (Int_t)ds.size();
  if (out.n > 1) {
    Double_t chi2 = 0;
    for (size_t i = 0; i < ds.size(); i++)
      chi2 += ws[i] * (ds[i] - out.delta) * (ds[i] - out.delta);
    out.chi2_ndf = chi2 / (out.n - 1);
    if (out.chi2_ndf > 1.0)
      out.sigma *= std::sqrt(out.chi2_ndf);
  }
  return out;
}

// E(x) = ref(x - delta) for the pol1 reference: the slope and its error stay,
// and the intercept moves by -p1 delta, its covariance with it.
CalFit OffsetShiftCal(const CalFit &ref, Double_t delta, const TString &label) {
  Double_t lo = 0, hi = 0;
  ref.func->GetRange(lo, hi);
  TF1 *shifted = new TF1("cal_low_" + label, "pol1", lo, hi);
  CalFit out;
  out.func = shifted;
  out.npar = 2;
  out.p0 = ref.p0 - ref.p1 * delta;
  out.p1 = ref.p1;
  out.var_p1 = ref.var_p1;
  out.var_p0 =
      ref.var_p0 + delta * delta * ref.var_p1 - 2.0 * delta * ref.cov_p0p1;
  out.cov_p0p1 = ref.cov_p0p1 - delta * ref.var_p1;
  shifted->SetParameter(0, out.p0);
  shifted->SetParameter(1, out.p1);
  shifted->SetParError(0, std::sqrt(std::max(0.0, out.var_p0)));
  shifted->SetParError(1, std::sqrt(std::max(0.0, out.var_p1)));
  shifted->SetNpx(1000);
  return out;
}

Float_t LoadRunStartTime(const TString &input_name) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name
              << ".root for run start time" << std::endl;
    if (file) {
      file->Close();
      delete file;
    }
    return 0;
  }
  Float_t t = 0;
  TTree *hdr = static_cast<TTree *>(file->Get("bef_header_tree"));
  if (hdr && hdr->GetEntries() > 0) {
    UInt_t timeStart = 0;
    hdr->SetBranchAddress("timeStart", &timeStart);
    hdr->GetEntry(0);
    t = (Float_t)timeStart;
  }
  file->Close();
  delete file;
  return t;
}

void ApplyPol1Cal(const std::vector<TString> &input_names, TF1 *cal,
                  const TString &date_label) {
  TFile *cal_file = IO::OpenForWriting(
      OUT_SUBDIR + "/calibration_function_low_" + date_label + ".root");
  cal->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  for (size_t i = 0; i < input_names.size(); i++) {
    TString input_name = input_names[i];
    TFile *in_file = IO::OpenForReading("filtered/" + input_name + ".root");
    if (!in_file || in_file->IsZombie()) {
      std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
                << std::endl;
      delete in_file;
      continue;
    }
    TFile *out_file =
        IO::OpenForWriting(OUT_SUBDIR + "/" + input_name + ".root");

    Float_t deposited_energy = 0;
    TTree *cal_tree = new TTree("calibrated_tree", "Calibrated events");
    cal_tree->SetDirectory(out_file);
    cal_tree->Branch("depositedEnergykeV", &deposited_energy,
                     "depositedEnergykeV/F");

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TTree *tree = static_cast<TTree *>(
          in_file->Get(Form("crystal%d_filtered_tree", c)));
      if (!tree)
        continue;
      Float_t energy = 0, x = 0, y = 0;
      tree->SetBranchAddress("energykeV", &energy);
      tree->SetBranchAddress("xmm", &x);
      tree->SetBranchAddress("ymm", &y);
      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        if (InteriorPixelsRequested() && !InInteriorPixels(x, y))
          continue;
        deposited_energy = cal->Eval(energy);
        cal_tree->Fill();
      }
    }

    out_file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    std::cout << "Calibrated " << input_name << " -> " << OUT_SUBDIR << "/"
              << input_name << ".root" << std::endl;
    out_file->Close();
    delete out_file;
    in_file->Close();
    delete in_file;
  }
}

void AddGeResult(DayResult &result, const BkgGeSimResult &r,
                 const TString &label, Float_t cal_err = 0,
                 Float_t gain_err = 0, Float_t bkg_err = 0, Float_t cal_off = 0,
                 Float_t anchor_err = 0) {
  if (!r.valid || r.sig_channel.peaks.empty())
    return;
  const PeakFitResult &ge = r.sig_channel.peaks.back();
  result.ge_labels.push_back(label);
  result.ge_mus.push_back(ge.mu);
  result.ge_errs.push_back(ge.mu_error);
  result.ge_chi2s.push_back(r.sig_channel.reduced_chi2);
  result.ge_cal_errs.push_back(cal_err);
  result.ge_cal_offs.push_back(cal_off);
  result.ge_gain_errs.push_back(gain_err);
  result.ge_anchor_errs.push_back(anchor_err);
  result.ge_bkg_errs.push_back(bkg_err);
}

Float_t LoadCrystalLiveTime(TFile *file, Int_t crystal) {
  TParameter<Float_t> *p = static_cast<TParameter<Float_t> *>(
      file->Get(Form("LiveTime_Filtered_Crystal%d_s", crystal)));
  if (!p) {
    std::cerr << "WARNING: missing LiveTime_Filtered_Crystal" << crystal << "_s"
              << std::endl;
    return 0;
  }
  return p->GetVal();
}

// Rate histograms (each crystal's counts over its live time, summed) with
// bins of width bin starting at n_phases points spaced bin / n_phases apart,
// all filled in one pass over the run. Histogram k covers
// [lo + k bin / n_phases, hi + k bin / n_phases).
std::vector<TH1F *> BuildRateHists(const TString &input_name, TF1 *cal,
                                   const TString &hname, Double_t lo,
                                   Double_t hi, Double_t bin, Int_t n_phases) {
  std::vector<TH1F *> sums;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    delete file;
    return sums;
  }
  Int_t nbins = (Int_t)std::lround((hi - lo) / bin);
  for (Int_t k = 0; k < n_phases; k++) {
    Double_t start = lo + k * bin / n_phases;
    TH1F *h = new TH1F(hname + Form("_phase%d", k),
                       ";Deposited Energy [keV]; Rate [counts/s]", nbins, start,
                       start + nbins * bin);
    h->SetDirectory(nullptr);
    h->Sumw2();
    sums.push_back(h);
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t lt = LoadCrystalLiveTime(file, c);
    if (lt <= 0)
      continue;
    Float_t energy = 0, x = 0, y = 0;
    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    Long64_t n = tree->GetEntries();
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      if (InteriorPixelsRequested() && !InInteriorPixels(x, y))
        continue;
      Double_t e = cal->Eval(energy);
      for (Int_t k = 0; k < n_phases; k++)
        sums[k]->Fill(e, 1.0 / lt);
    }
    tree->ResetBranchAddresses();
  }
  file->Close();
  delete file;
  return sums;
}

// Line-cal day: the Ge tail interpolated linearly in energy between the
// Am-241 and Ba-133 80.98 anchor fits. They are in raw units, but tau/sigma
// and the height ratio do not depend on the energy scale.
RateSubShape RateSubShapeFromLines(const PeakFitResult &am,
                                   const PeakFitResult &ba) {
  RateSubShape s;
  Double_t w = (E_NEAR_GE - E_AM241) / (E_BA133_81 - E_AM241);
  Double_t k_am = TailHeightRatio(am);
  Double_t k_ba = TailHeightRatio(ba);
  Double_t rho_am =
      am.low_exp_tail_ratio > 0 ? am.low_exp_tail_ratio : ba.low_exp_tail_ratio;
  Double_t rho_ba = ba.low_exp_tail_ratio > 0 ? ba.low_exp_tail_ratio : rho_am;
  s.ge_tail_height_ratio = k_am + w * (k_ba - k_am);
  s.ge_tau_over_sigma = rho_am + w * (rho_ba - rho_am);
  if (!(s.ge_tau_over_sigma > 0)) {
    s.ge_tail_height_ratio = 0;
    s.ge_tau_over_sigma = 1;
  }
  return s;
}

void PrintRateSubShape(const TString &tag, const RateSubShape &s) {
  std::cout << "RATESUB-SHAPE " << tag << ": window " << std::fixed
            << std::setprecision(1) << RATESUB_LO << "-" << RATESUB_HI
            << " keV   Ge tail/gaus height " << std::setprecision(4)
            << s.ge_tail_height_ratio << "   tau/sigma " << std::setprecision(3)
            << s.ge_tau_over_sigma << std::endl;
}

// Method 2 for one day: the residual's 95.902 keV line (for its escapes) and
// the phase-averaged rate-sub fit of the live-time residual, signal minus
// background, from each run's rate histograms (FitRateSubPhases in
// GainTransferFit.hpp), then the plot of the first starting point's fit.
RateSubResult FitRateSubtractedGe(const TString &sig_run,
                                  const TString &bkg_run, TF1 *cal,
                                  const TString &tag,
                                  const RateSubShape &tail_shape) {
  const Float_t lo = RATESUB_LO, hi = RATESUB_HI;
  PrintRateSubShape(tag, tail_shape);

  const Double_t bin = COMB_PERIOD_RAW_KEV * cal->GetParameter(1);
  RateSubShape shape = tail_shape;
  if (EscapeLinesRequested()) {
    std::vector<TH1F *> capture_sigs = BuildRateHists(
        sig_run, cal, "capture_sig_" + tag, 88.0, 104.0, bin, N_RATESUB_PHASES);
    std::vector<TH1F *> capture_bkgs = BuildRateHists(
        bkg_run, cal, "capture_bkg_" + tag, 88.0, 104.0, bin, N_RATESUB_PHASES);
    SetResidualCaptureLine(capture_sigs, capture_bkgs, tag, shape);
    for (size_t k = 0; k < capture_sigs.size(); k++)
      delete capture_sigs[k];
    for (size_t k = 0; k < capture_bkgs.size(); k++)
      delete capture_bkgs[k];
  }
  std::vector<TH1F *> sigs =
      BuildRateHists(sig_run, cal, "rate_sig_" + tag, lo - bin, hi + bin, bin,
                     N_RATESUB_PHASES);
  std::vector<TH1F *> bkgs =
      BuildRateHists(bkg_run, cal, "rate_bkg_" + tag, lo - bin, hi + bin, bin,
                     N_RATESUB_PHASES);
  RateSubResult out;
  TF1 *model = nullptr;
  TH1F *res = nullptr;
  if ((Int_t)sigs.size() == N_RATESUB_PHASES &&
      (Int_t)bkgs.size() == N_RATESUB_PHASES)
    out = FitRateSubPhases(sigs, bkgs, shape, tag, &model, &res);
  for (size_t k = 0; k < sigs.size(); k++)
    delete sigs[k];
  for (size_t k = 0; k < bkgs.size(); k++)
    delete bkgs[k];
  if (!out.valid || !model || !res) {
    delete model;
    delete res;
    return out;
  }

  res->GetYaxis()->SetTitle(
      Form("Counts/s / %.0f eV", 1000.0 * res->GetBinWidth(1)));
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(res, kBlack);
  PlottingUtils::SaveFigure(canvas, "rate_residual_" + tag, "fits",
                            PlotSaveOptions::kLINEAR);

  const Int_t kNPts = 400;
  Double_t xstep = (hi - lo) / (kNPts - 1);
  TGraph *total_graph = new TGraph(kNPts);
  TGraph *bkg_graph = new TGraph(kNPts);
  TGraph *escape_graph = new TGraph(kNPts);
  for (Int_t i = 0; i < kNPts; i++) {
    Double_t xv = lo + i * xstep;
    Double_t bkg_v =
        model->GetParameter(5) + model->GetParameter(6) * (xv - 68.75);
    total_graph->SetPoint(i, xv, model->Eval(xv));
    bkg_graph->SetPoint(i, xv, bkg_v);
    escape_graph->SetPoint(i, xv,
                           bkg_v + model->GetParameter(8) *
                                       EscapeDensity(xv, model->GetParameter(9),
                                                     model->GetParameter(2),
                                                     model->GetParameter(3),
                                                     model->GetParameter(4),
                                                     ESCAPES_INSIDE));
  }
  bkg_graph->SetLineStyle(2);
  escape_graph->SetLineColor(kRed);
  escape_graph->SetLineStyle(2);
  std::vector<TGraph *> components = {bkg_graph};
  if (model->GetParameter(8) != 0)
    components.push_back(escape_graph);

  res->GetXaxis()->SetRangeUser(lo, hi);
  Double_t ymin = 1e30, ymax = -1e30;
  for (Int_t b = res->FindBin(lo); b <= res->FindBin(hi - 1e-6); b++) {
    Double_t v = res->GetBinContent(b);
    Double_t e = res->GetBinError(b);
    if (v - e < ymin)
      ymin = v - e;
    if (v + e > ymax)
      ymax = v + e;
  }
  Double_t pad = 0.12 * (ymax - ymin);
  res->SetMinimum(ymin - pad);
  res->SetMaximum(ymax + pad);
  PlottingUtils::PlotFitWithResiduals(res, total_graph, components, lo, hi,
                                      "ratesub_fit_" + tag, "fits",
                                      "Rate-subtracted " + tag, kFALSE);
  delete res;
  return out;
}

// Map seed peaks from RAW units into CALIBRATED units through cal. Under
// E = f(x) with s = f'(mu), a constant for a straight-line calibration:
//   sigma (a length in x)            -> multiply by s
//   LowLinTailSlope (per unit x)     -> divide by s
//   *TailRatio (tau/sigma), *TailAmplitude, step -> unchanged: tau and sigma
//   scale together, and the amplitudes are ratios
void ConvertSeedPeaksToCalibrated(TF1 *cal, const FitResult &src,
                                  FitResult &dst) {
  for (size_t i = 0; i < dst.peaks.size() && i < src.peaks.size(); i++) {
    Double_t mu_raw = src.peaks[i].mu;
    Double_t s = cal->Derivative(mu_raw);
    if (!(s > 0))
      s = 1.0;
    dst.peaks[i].mu = cal->Eval(mu_raw);
    dst.peaks[i].sigma = src.peaks[i].sigma * s;
    dst.peaks[i].low_exp_tail_ratio = src.peaks[i].low_exp_tail_ratio;
    dst.peaks[i].high_exp_tail_ratio = src.peaks[i].high_exp_tail_ratio;
    dst.peaks[i].low_lin_tail_slope = src.peaks[i].low_lin_tail_slope / s;
  }
}

TString SanitizeLabel(TString label) {
  label.ReplaceAll(" ", "_");
  return label;
}

struct SignalPbAnchor {
  CalFit calfit;
  Bool_t valid = kFALSE;
  Double_t offset = 0;     // keV taken off the pair calibration
  Double_t offset_err = 0; // keV, at the signal run's Pb lines
};

// The pair calibration shifted by a constant so that the signal run's own Pb
// Ka lines land on their energies, weighted by their centroid errors. The
// slope stays the pair's, and its error only acts over the few keV from the Pb
// lines down to Ge; applied as a gain instead, the shift (at most ~30 eV)
// would move Ge by at most ~2 eV. The covariance is set so that
// var E(x) = offset_err^2 + var_p1 (x - x_pb)^2.
SignalPbAnchor AnchorOnSignalPb(const CalFit &pair_cal,
                                const FitResult &sig_channel,
                                const std::vector<Double_t> &pb_mus,
                                const TString &cal_label) {
  SignalPbAnchor out;
  if (!pair_cal.func || sig_channel.peaks.size() <= pb_mus.size())
    return out;
  Double_t sum_w = 0, sum_wr = 0, sum_wx = 0;
  for (size_t i = 0; i < pb_mus.size(); i++) {
    const PeakFitResult &pk = sig_channel.peaks[i];
    Double_t err_kev = pair_cal.p1 * GuardedMuError(cal_label + " signal Pb-Ka",
                                                    pk, sig_channel);
    if (!(err_kev > 0))
      return out;
    Double_t calibrated = pair_cal.func->Eval(pk.mu);
    Double_t resid = calibrated - pb_mus[i];
    Double_t w = 1.0 / (err_kev * err_kev);
    sum_w += w;
    sum_wr += w * resid;
    sum_wx += w * pk.mu;
    std::cout << "SIGPB " << cal_label << ": signal-run Pb line " << std::fixed
              << std::setprecision(4) << pb_mus[i] << " keV, raw " << pk.mu
              << " -> " << calibrated << " keV on the pair calibration ("
              << std::setprecision(1) << resid * 1000.0 << " +/- "
              << err_kev * 1000.0 << " eV)" << std::endl;
  }
  out.offset = sum_wr / sum_w;
  out.offset_err = std::sqrt(1.0 / sum_w + dE_PB_KA * dE_PB_KA);
  Double_t x_pb = sum_wx / sum_w;

  Double_t lo = 0, hi = 0;
  pair_cal.func->GetRange(lo, hi);
  TF1 *shifted =
      new TF1(TString(pair_cal.func->GetName()) + "_sigPb", "pol1", lo, hi);
  shifted->SetParameter(0, pair_cal.p0 - out.offset);
  shifted->SetParameter(1, pair_cal.p1);
  shifted->SetNpx(1000);
  out.calfit.func = shifted;
  out.calfit.npar = 2;
  out.calfit.p0 = pair_cal.p0 - out.offset;
  out.calfit.p1 = pair_cal.p1;
  out.calfit.var_p1 = pair_cal.var_p1;
  out.calfit.var_p0 =
      out.offset_err * out.offset_err + pair_cal.var_p1 * x_pb * x_pb;
  out.calfit.cov_p0p1 = -pair_cal.var_p1 * x_pb;
  shifted->SetParError(0, std::sqrt(out.calfit.var_p0));
  shifted->SetParError(1, std::sqrt(out.calfit.var_p1));
  out.valid = kTRUE;
  std::cout << "SIGPB " << cal_label << ": signal run shifted by " << std::fixed
            << std::setprecision(1) << -out.offset * 1000.0 << " +/- "
            << out.offset_err * 1000.0 << " eV relative to the pair calibration"
            << std::endl;
  return out;
}

// escapes: the pair's 95.902 keV escape fractions, null to leave them out.
PairCalResult ProcessPbAnchoredPair(
    const TString &bkg_run, const TString &sig_run, const TString &date_label,
    const TString &pair_tag, const TString &am_run, Bool_t is_reference,
    const EscapeFractions *escapes, const CalFit *ref_cal = nullptr,
    const std::vector<Float_t> *ref_pb_mus = nullptr,
    const std::vector<Float_t> *ref_pb_errs = nullptr) {
  PairCalResult out;
  out.pair_tag = pair_tag;

  std::vector<Double_t> bkg_events = LoadFiltered(bkg_run);
  std::vector<Double_t> sig_events = LoadFiltered(sig_run);
  std::vector<Double_t> pb_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  TString cal_label = date_label + "_" + pair_tag;

  FitResult am;
  Bool_t have_am = kFALSE;
  if (is_reference && !am_run.IsNull()) {
    std::vector<Double_t> am_events = LoadFiltered(am_run);
    am = FitCalPeak(am_events, am_run, "Am_59.5keV", 51, 71, kFALSE, 0, 0);
    if (am.valid && !am.peaks.empty()) {
      const PeakFitResult &p = am.peaks.at(0);
      std::cout << "AM-SHAPE chi2/ndf " << am.reduced_chi2 << "  sigma "
                << p.sigma << "  step " << p.step_amplitude << "  lowExpAmp "
                << p.low_exp_tail_amplitude << "  lowExpTau "
                << p.low_exp_tail_ratio << "  lowLinAmp "
                << p.low_lin_tail_amplitude << "  lowLinSlope "
                << p.low_lin_tail_slope << std::endl;
    }
    have_am = am.valid;
    if (!have_am)
      std::cerr << "WARNING: Am-241 fit failed for " << am_run
                << "; reference line will use Pb peaks only" << std::endl;
  }

  InSituEscapes raw_escapes;
  raw_escapes.fractions = escapes;
  raw_escapes.units_per_kev = 1.0 / RAW_KEV_PER_UNIT_0113;
  raw_escapes.comb_period = COMB_PERIOD_RAW_KEV;
  BkgGeSimResult pre =
      RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_run, sig_run,
                  "PreCal_" + pair_tag, nullptr, nullptr, raw_escapes);
  if (!pre.valid)
    return out;

  if (pre.bkg_channel.peaks.size() >= 2) {
    Double_t mu_lo = pre.bkg_channel.peaks.at(0).mu;
    Double_t mu_hi = pre.bkg_channel.peaks.at(1).mu;
    Double_t g = (mu_lo + mu_hi) / (Double_t)(E_PB_KA2 + E_PB_KA1);
    Double_t delta_raw = (Double_t)(E_PB_KA1 - E_PB_KA2) * g;
    std::cout << "PRECAL-SEP " << pair_tag << ": fitted " << std::fixed
              << std::setprecision(5) << (mu_hi - mu_lo) << "  expected "
              << delta_raw << "  (excess " << std::setprecision(1)
              << ((mu_hi - mu_lo) - delta_raw) * 1000 << " eV)"
              << std::setprecision(4) << std::endl;
  }

  for (size_t i = 0; i < pre.bkg_channel.peaks.size(); i++) {
    out.pb_mus_precal.push_back(pre.bkg_channel.peaks[i].mu);
    out.pb_mu_errs_precal.push_back(GuardedMuError(
        pair_tag + " Pb-Ka", pre.bkg_channel.peaks[i], pre.bkg_channel));
  }

  CalibrationData cal_data;
  AddCalPoint(cal_data, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
              GuardedMuError(pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0),
                             pre.bkg_channel),
              E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
  AddCalPoint(cal_data, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
              GuardedMuError(pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1),
                             pre.bkg_channel),
              E_PB_KA1, -1, dE_PB_KA);
  Bool_t transfer_mode =
      (ref_cal != nullptr && ref_pb_mus != nullptr && ref_pb_mus->size() >= 2 &&
       pre.bkg_channel.peaks.size() >= 2);
  Bool_t use_am_anchor = (!transfer_mode && have_am);

  CalFit calfit;
  GainResult gain;
  OffsetResult offset;

  if (transfer_mode && OffsetTransferRequested()) {
    std::vector<Float_t> mu_run = {(Float_t)pre.bkg_channel.peaks.at(0).mu,
                                   (Float_t)pre.bkg_channel.peaks.at(1).mu};
    std::vector<Float_t> err_run = {
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka2",
                                pre.bkg_channel.peaks.at(0), pre.bkg_channel),
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka1",
                                pre.bkg_channel.peaks.at(1), pre.bkg_channel)};
    std::vector<Float_t> er =
        ref_pb_errs ? *ref_pb_errs
                    : std::vector<Float_t>(ref_pb_mus->size(), 0.f);
    offset = PbOffsetVsReference(mu_run, err_run, *ref_pb_mus, er);
    calfit = OffsetShiftCal(*ref_cal, offset.delta, cal_label);
    std::cout << "OFFSET-TRANSFER " << cal_label << ": delta = " << std::fixed
              << std::setprecision(1) << offset.delta * 1000.0 << " +/- "
              << offset.sigma * 1000.0 << " eV raw  (n = " << offset.n
              << ", chi2/ndf = " << std::setprecision(2) << offset.chi2_ndf
              << ")" << std::endl;
  } else if (transfer_mode) {
    std::vector<Float_t> mu_run = {(Float_t)pre.bkg_channel.peaks.at(0).mu,
                                   (Float_t)pre.bkg_channel.peaks.at(1).mu};
    std::vector<Float_t> err_run = {
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka2",
                                pre.bkg_channel.peaks.at(0), pre.bkg_channel),
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka1",
                                pre.bkg_channel.peaks.at(1), pre.bkg_channel)};
    std::vector<Float_t> er =
        ref_pb_errs ? *ref_pb_errs
                    : std::vector<Float_t>(ref_pb_mus->size(), 0.f);

    gain = PbGainVsReferenceWithError(mu_run, err_run, *ref_pb_mus, er);
    calfit = GainScaleCal(*ref_cal, gain.g, cal_label);
    std::cout << "GAIN-TRANSFER " << cal_label << ": g = " << std::fixed
              << std::setprecision(6) << gain.g << " +/- " << gain.sigma_g
              << "  (n = " << gain.n << ", chi2/ndf = " << std::setprecision(2)
              << gain.chi2_ndf << ")" << std::endl;
  } else if (use_am_anchor) {
    CalibrationData cal_am;
    AddCalPoint(cal_am, pair_tag + " Am-241 59.5", am.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Am-241", am.peaks.at(0), am),
                E_AM241, am.reduced_chi2, dE_AM241);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Ka2",
                               pre.bkg_channel.peaks.at(0), pre.bkg_channel),
                E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
                GuardedMuError(pair_tag + " Pb-Ka1",
                               pre.bkg_channel.peaks.at(1), pre.bkg_channel),
                E_PB_KA1, -1, dE_PB_KA);
    PrintCalSummary(cal_am, cal_label);
    calfit = CreateAndSavePol1Cal(cal_am, cal_label);
  } else {
    PrintCalSummary(cal_data, cal_label);
    calfit = CreateAndSavePol1Cal(cal_data, cal_label);
  }
  TF1 *cal = calfit.func;

  CalFit pb_probe = calfit;
  if (have_am)
    pb_probe = CreateAndSavePol1Cal(cal_data, cal_label + "_PbOnly");

  if (have_am) {
    out.ref_calfit = calfit;
    out.has_ref_cal = kTRUE;
  }
  out.cal_func = cal;

  SignalPbAnchor sig_anchor;
  if (SignalPbAnchorRequested()) {
    sig_anchor = AnchorOnSignalPb(calfit, pre.sig_channel, pb_mus, cal_label);
    if (!sig_anchor.valid)
      std::cerr << "WARNING: no signal-run Pb anchor for " << cal_label
                << "; the signal run keeps the pair calibration" << std::endl;
  }
  const CalFit &sig_calfit = sig_anchor.valid ? sig_anchor.calfit : calfit;
  TF1 *sig_cal = sig_calfit.func;

  if (!pre.sig_channel.peaks.empty() && sig_anchor.valid) {
    // The signal run sets its own offset, so it inherits no reference offset
    // (cal_off 0) and no transfer (gain 0); the anchor column carries the
    // signal-run Pb offset error, and the rest of cal_err is the pair slope
    // over the Pb-to-Ge distance.
    Double_t mu_ge = pre.sig_channel.peaks.back().mu;
    out.ge_cal_err = CalEnergyError(sig_calfit, mu_ge);
    out.ge_gain_err = 0;
    out.ge_anchor_err = sig_anchor.offset_err;
    out.ge_cal_off = 0;
    std::cout << "SIGPB-CAL " << cal_label << ": at Ge, offset term "
              << std::fixed << std::setprecision(1)
              << out.ge_anchor_err * 1000.0 << " eV   total cal "
              << out.ge_cal_err * 1000.0 << " eV" << std::endl;
  } else if (!pre.sig_channel.peaks.empty()) {
    Double_t mu_ge = pre.sig_channel.peaks.back().mu;
    out.ge_cal_err = CalEnergyError(calfit, mu_ge);
    out.ge_cal_off = CalOffsetError(calfit);

    if (transfer_mode) {
      out.ge_gain_err =
          OffsetTransferRequested()
              ? calfit.p1 * offset.sigma
              : GainTransferEnergyError(*ref_cal, mu_ge, gain.g, gain.sigma_g);
      out.ge_cal_err = std::hypot(out.ge_cal_err, out.ge_gain_err);

      std::cout << "GAIN-CAL " << cal_label << ": transfer term " << std::fixed
                << std::setprecision(1) << out.ge_gain_err * 1000.0
                << " eV   total cal " << out.ge_cal_err * 1000.0 << " eV"
                << std::endl;
    }
  }

  // Cross-checks on the pair that also has Am, against a calibration through
  // its Pb lines alone. Both are evaluated at the same raw Ge centroid, so only
  // the calibrations differ; and the Pb-only line predicts Am, which plays no
  // part in it, so the LINEARITY check is not circular.
  if (have_am && !pre.sig_channel.peaks.empty()) {
    Double_t ge_pre = pre.sig_channel.peaks.back().mu;
    Double_t e_nominal = cal->Eval(ge_pre);
    Double_t e_pbonly = pb_probe.func->Eval(ge_pre);
    std::cout << "VALIDATION " << pair_tag << " (Ge precal " << std::fixed
              << std::setprecision(4) << ge_pre << "):  nominal -> "
              << e_nominal << " keV   Pb-only -> " << e_pbonly
              << " keV   (nominal - Pb-only) " << std::setprecision(1)
              << (e_nominal - e_pbonly) * 1000.0 << " eV" << std::endl;

    Double_t am_raw = am.peaks.at(0).mu;
    Double_t am_pred = pb_probe.func->Eval(am_raw);
    Double_t sig_cal = CalEnergyError(pb_probe, am_raw);
    Double_t sig_amfit = pb_probe.p1 * am.peaks.at(0).mu_error;
    Double_t sig_pred = std::sqrt(sig_cal * sig_cal + sig_amfit * sig_amfit);
    Double_t resid = am_pred - E_AM241;
    std::cout << "LINEARITY " << pair_tag
              << ": Pb-only line predicts Am = " << std::fixed
              << std::setprecision(4) << am_pred << " +/- " << sig_pred
              << " keV   (truth " << E_AM241 << ", residual "
              << std::setprecision(1) << resid * 1000.0
              << " eV = " << std::setprecision(2)
              << (sig_pred > 0 ? resid / sig_pred : 0) << " sigma)"
              << std::endl;
  }

  std::vector<TString> pair_runs = {bkg_run};
  if (!sig_anchor.valid)
    pair_runs.push_back(sig_run);
  if (have_am)
    pair_runs.push_back(am_run);
  ApplyPol1Cal(pair_runs, cal, cal_label);
  if (sig_anchor.valid) {
    std::vector<TString> sig_runs = {sig_run};
    ApplyPol1Cal(sig_runs, sig_cal, cal_label + "_sigPb");
  }

  std::vector<Double_t> bc = LoadCalibrated(bkg_run);
  std::vector<Double_t> sc = LoadCalibrated(sig_run);

  Double_t cal_p0 = cal->GetParameter(0);
  Double_t cal_p1 = cal->GetParameter(1);
  Double_t sig_cal_p0 = sig_cal->GetParameter(0);
  FitResult bkg_seed_post = pre.bkg_channel;
  FitResult sig_seed_post = pre.sig_channel;
  if (cal_p1 > 0) {
    Double_t s_b = pre.bkg_channel.lin_bkg_slope;
    Double_t s_s = pre.sig_channel.lin_bkg_slope;
    bkg_seed_post.lin_bkg_slope = s_b / (cal_p1 - s_b * cal_p0);
    sig_seed_post.lin_bkg_slope = s_s / (cal_p1 - s_s * sig_cal_p0);
    bkg_seed_post.bkg_constant = pre.bkg_channel.bkg_constant / cal_p1;
    sig_seed_post.bkg_constant = pre.sig_channel.bkg_constant / cal_p1;
  }
  // Map the seed peak mus through the calibration so the Pb and Ge peaks start
  // at their calibrated positions, and convert the SHAPE into calibrated units
  // as well: sigma and the tail lengths were fitted in RAW units and carry the
  // local derivative f'(mu) ~ 0.989-0.995, so handing them over unconverted
  // makes every width ~1% too large (sigma wide by 7-10 eV).
  ConvertSeedPeaksToCalibrated(cal, pre.bkg_channel, bkg_seed_post);
  ConvertSeedPeaksToCalibrated(sig_cal, pre.sig_channel, sig_seed_post);

  // Pb peaks left free post-cal: where the calibrated Pb doublet actually
  // lands is a blind check at the energy of interest. The whole chain (precal
  // -> cal -> postcal) uses one shared tail model, so this run is internally
  // consistent.
  InSituEscapes calibrated_escapes;
  calibrated_escapes.fractions = escapes;
  calibrated_escapes.units_per_kev = 1.0;
  calibrated_escapes.comb_period = COMB_PERIOD_RAW_KEV * cal_p1;
  out.postcal =
      RunBkgGeSim(bc, sc, pb_mus, bkg_run + "_postcal_" + pair_tag,
                  sig_run + "_postcal_" + pair_tag, "PostCal_" + pair_tag,
                  &bkg_seed_post, &sig_seed_post, calibrated_escapes);

  if (out.postcal.valid) {
    std::cout << "Post-cal FREE Pb mu for " << cal_label << ":" << std::endl;
    for (size_t i = 0;
         i < out.postcal.bkg_channel.peaks.size() && i < pb_mus.size(); i++)
      std::cout << "  Pb peak " << i << ": " << std::fixed
                << std::setprecision(4) << out.postcal.bkg_channel.peaks[i].mu
                << " keV (truth " << pb_mus[i] << ", shift "
                << (out.postcal.bkg_channel.peaks[i].mu - pb_mus[i]) << ")"
                << std::endl;

    if (!pre.sig_channel.peaks.empty()) {
      Double_t ge_precal_mapped =
          sig_cal->Eval(pre.sig_channel.peaks.back().mu);
      Double_t ge_postcal = out.postcal.sig_channel.peaks.back().mu;
      out.ge_bkg_err = std::fabs(ge_precal_mapped - ge_postcal);
      std::cout << "BACKGROUND " << cal_label
                << ": cal(precal Ge)=" << std::fixed << std::setprecision(4)
                << ge_precal_mapped << "  postcal Ge=" << ge_postcal
                << "  |diff|=" << std::setprecision(1)
                << out.ge_bkg_err * 1000.0 << " eV" << std::endl;
    }
  }
  return out;
}

void AppendPairToDay(DayResult &result, const PairCalResult &p,
                     const TString &date_label, const TString &ge_label) {
  if (p.cal_func) {
    result.cal_labels.push_back(date_label + " " + p.pair_tag);
    result.cal_funcs.push_back(p.cal_func);
  }
  AddGeResult(result, p.postcal, ge_label, p.ge_cal_err, p.ge_gain_err,
              p.ge_bkg_err, p.ge_cal_off, p.ge_anchor_err);
}

DayResult ProcessLineCalDay(const LineCalConfig &cfg) {
  CalibrationData cal_data;

  std::vector<Double_t> am_events = LoadFiltered(cfg.am_run);
  FitResult am = FitCalPeak(am_events, cfg.am_run, "Am_59.5keV", cfg.am_lo,
                            cfg.am_hi, kFALSE, 0, 0);
  if (am.valid)
    AddCalPoint(cal_data, "Am-241 59.5", am.peaks.at(0).mu,
                am.peaks.at(0).mu_error, E_AM241, am.reduced_chi2, dE_AM241);

  // 79.61 and 81.00 keV are 1.4 keV apart at sigma 0.82: two free tails are
  // degenerate (the weak line's railed, the strong line's came out at
  // tau/sigma 13.7 +/- 15.3), so the doublet carries one width and one tail.
  // With the yields free as well, the weak line takes 8% (01/15) and 33%
  // (01/16) more than the intensity ratio allows, trading against the shared
  // tail; the fixed ratio leaves chi2/ndf where it was.
  std::vector<Double_t> ba_events = LoadFiltered(cfg.ba_run);
  const Float_t ba_lo = 75, ba_hi = 90;
  FitResult ba81 =
      FitCalPeak(ba_events, cfg.ba_run, "Ba_80.98keV", ba_lo, ba_hi, kTRUE,
                 E_BA133_79, E_BA133_81, kTRUE, kTRUE, BA133_79_OVER_81);
  if (ba81.valid)
    AddCalPoint(cal_data, "Ba-133 80.98", ba81.peaks.at(1).mu,
                ba81.peaks.at(1).mu_error, E_BA133_81, ba81.reduced_chi2,
                dE_LINE);

  PrintCalSummary(cal_data, cfg.date_label);
  CalFit calfit = CreateAndSavePol1Cal(cal_data, cfg.date_label);
  TF1 *cal = calfit.func;
  Double_t mu_at_ge =
      (calfit.p1 != 0) ? (E_NEAR_GE - calfit.p0) / calfit.p1 : E_NEAR_GE;
  Double_t ge_cal_err = CalEnergyError(calfit, mu_at_ge);
  Double_t ge_cal_off = CalOffsetError(calfit);
  ApplyPol1Cal(cfg.day_datasets, cal, cfg.date_label);

  DayResult result;
  result.cal_labels.push_back(cfg.date_label);
  result.cal_funcs.push_back(cal);

  // Without the Am and Ba tails the default shape, a Gaussian, is used.
  RateSubShape shape;
  if (am.valid && ba81.valid && ba81.peaks.size() > 1)
    shape = RateSubShapeFromLines(am.peaks.at(0), ba81.peaks.at(1));
  else
    std::cerr << "WARNING: no Am/Ba tail for " << cfg.date_label
              << "; rate-sub falls back to a Gaussian" << std::endl;
  RateSubResult rs = FitRateSubtractedGe(cfg.postcal_sig, cfg.postcal_bkg, cal,
                                         cfg.date_label, shape);
  if (rs.valid) {
    result.rs_labels.push_back(cfg.ge_label + " [rate-sub]");
    result.rs_mus.push_back(rs.mu);
    result.rs_errs.push_back(rs.mu_err);
    result.rs_chi2s.push_back(rs.chi2);
    result.rs_cal_errs.push_back(ge_cal_err);
    result.rs_cal_offs.push_back(ge_cal_off);
    result.rs_gain_errs.push_back(0.0f);
    result.rs_sys_errs.push_back(rs.sys_err);
    std::cout << "Rate-sub Ge mu for " << cfg.date_label << ": " << std::fixed
              << std::setprecision(4) << rs.mu << " +/- " << rs.mu_err
              << " (fit) +/- " << rs.sys_err
              << " (sys) keV (chi2/ndf = " << std::setprecision(3) << rs.chi2
              << ")" << std::endl;
  }
  return result;
}

struct ShieldPair {
  TString bkg_run;
  TString sig_run;
  TString tag;
  TString ge_label;
  // Where the 95.902 keV line comes from: captures inside the crystal, or in
  // the Cd shield.
  const EscapeFractions *escapes = nullptr;
};

const EscapeFractions *PairEscapes(const ShieldPair &pair) {
  return EscapeLinesRequested() ? pair.escapes : nullptr;
}

DayResult ProcessDay_20260113(PairCalResult &ref_pair_out) {
  DayResult result;

  std::vector<ShieldPair> pairs;
  ShieldPair cu;
  cu.bkg_run = Constants::CUSHIELDBACKGROUND_10PERCENT_20260113;
  cu.sig_run = Constants::CUSHIELDSIGNAL_10PERCENT_20260113;
  cu.tag = "CuShield10";
  cu.ge_label = "Cu Shield Signal 10% (01/13)";
  cu.escapes = &ESCAPES_INSIDE;
  pairs.push_back(cu);
  ShieldPair cd10;
  cd10.bkg_run = Constants::CDSHIELDBACKGROUND_10PERCENT_20260113;
  cd10.sig_run = Constants::CDSHIELDSIGNAL_10PERCENT_20260113;
  cd10.tag = "CdShield10";
  cd10.ge_label = "Cd Shield Signal 10% (01/13)";
  cd10.escapes = &ESCAPES_OUTSIDE;
  pairs.push_back(cd10);
  ShieldPair cd25;
  cd25.bkg_run = Constants::CDSHIELDBACKGROUND_25PERCENT_20260113;
  cd25.sig_run = Constants::CDSHIELDSIGNAL_25PERCENT_20260113;
  cd25.tag = "CdShield25";
  cd25.ge_label = "Cd Shield Signal 25% (01/13)";
  cd25.escapes = &ESCAPES_OUTSIDE;
  pairs.push_back(cd25);

  TString am_run = Constants::POSTREACTOR_AM241_20260113;
  Float_t am_time = LoadRunStartTime(am_run);
  Int_t ref_idx = 0;
  Float_t best_dt = 1e30f;
  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    Float_t dt = TMath::Abs(LoadRunStartTime(pairs[i].bkg_run) - am_time);
    if (dt < best_dt) {
      best_dt = dt;
      ref_idx = i;
    }
  }
  std::cout << "Reference pair: " << pairs[ref_idx].tag << " (Am from "
            << am_run << ")" << std::endl;

  std::vector<PairCalResult> results(pairs.size());

  results[ref_idx] = ProcessPbAnchoredPair(
      pairs[ref_idx].bkg_run, pairs[ref_idx].sig_run, "20260113",
      pairs[ref_idx].tag, am_run, kTRUE, PairEscapes(pairs[ref_idx]));

  const CalFit *ref_cal =
      results[ref_idx].has_ref_cal ? &results[ref_idx].ref_calfit : nullptr;
  const std::vector<Float_t> *ref_mus = &results[ref_idx].pb_mus_precal;
  const std::vector<Float_t> *ref_errs = &results[ref_idx].pb_mu_errs_precal;

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    if (i == ref_idx)
      continue;
    results[i] = ProcessPbAnchoredPair(
        pairs[i].bkg_run, pairs[i].sig_run, "20260113", pairs[i].tag, "",
        kFALSE, PairEscapes(pairs[i]), ref_cal, ref_mus, ref_errs);
  }

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++)
    AppendPairToDay(result, results[i], "20260113", pairs[i].ge_label);

  ref_pair_out = results[ref_idx];
  return result;
}

DayResult ProcessDay_20260115() {
  LineCalConfig cfg;
  cfg.date_label = "20260115";
  cfg.am_run = Constants::POSTREACTOR_AM241_20260115;
  cfg.am_lo = 51;
  cfg.am_hi = 71;
  cfg.ba_run = Constants::POSTREACTOR_BA133_20260115;
  cfg.day_datasets = {Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::POSTREACTOR_AM241_20260115,
                      Constants::POSTREACTOR_BA133_20260115,
                      Constants::SHUTTERCLOSED_20260115};
  cfg.postcal_bkg = Constants::NOSHIELDBACKGROUND_5PERCENT_20260115;
  cfg.postcal_sig = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;
  cfg.ge_label = "No Shield Signal 5% (01/15)";
  return ProcessLineCalDay(cfg);
}

DayResult ProcessDay_20260116() {
  LineCalConfig cfg;
  cfg.date_label = "20260116";
  cfg.am_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.am_lo = 55;
  cfg.am_hi = 70;
  cfg.ba_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.day_datasets = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  cfg.postcal_bkg =
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116;
  cfg.postcal_sig = Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116;
  cfg.ge_label = "Graphite Castle Signal 10% (01/16)";
  return ProcessLineCalDay(cfg);
}

void PrintDayCals(const DayResult &d) {
  for (size_t i = 0; i < d.cal_funcs.size(); i++) {
    TF1 *c = d.cal_funcs[i];
    if (!c)
      continue;
    std::cout << "Cal pol1 " << d.cal_labels[i] << ": p0 = " << std::fixed
              << std::setprecision(6) << c->GetParameter(0) << " +/- "
              << c->GetParError(0) << ", p1 = " << c->GetParameter(1) << " +/- "
              << c->GetParError(1) << std::endl;
  }
}

void CalibrateAndFit() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  const TString out_tag = RuntimeResultTag();
  std::cout << "OUTPUT-TAG " << out_tag
            << (out_tag == RESULT_TAG ? " (default)" : " (RESULT_TAG_OVERRIDE)")
            << std::endl;
  std::cout << "CONFIG tail tau/sigma <= " << TAIL_RATIO_MAX << "   hard fit "
            << (HardFitRequested() ? "on" : "off") << "   Pb anchor "
            << (SignalPbAnchorRequested() ? "signal run" : "background run")
            << "   transfer " << (OffsetTransferRequested() ? "offset" : "gain")
            << "   pixels "
            << (InteriorPixelsRequested() ? "interior 9 x 9" : "all")
            << "   in-situ background " << (USE_FLAT_BKG ? "flat" : "linear")
            << "   95.902 keV escapes "
            << (EscapeLinesRequested() ? "on" : "off") << "   ROOT "
            << gROOT->GetVersion() << std::endl;
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG, project_root + "/plots/calibration_low_" + out_tag,
      project_root + "/root_files");
  PairCalResult ref_pair;
  DayResult d13 = ProcessDay_20260113(ref_pair);
  DayResult d15 = ProcessDay_20260115();
  DayResult d16 = ProcessDay_20260116();

  PrintDayCals(d13);
  PrintDayCals(d15);
  PrintDayCals(d16);

  // 15th/16th: line-cal (Am + Ba-133) days measured by rate subtraction.
  // Independent of the 13th in both calibration source and background
  // treatment, which is what makes the method systematic meaningful.
  // One .result row per run; no budget math here.
  // The 01/14 Cu 10% pair is not used. The detector was probably unbiased
  // overnight and re-biased that morning, and there was no Am run on the 14th,
  // so neither a gain nor an offset carries the 01/13 calibration across.
  std::vector<DayResult *> day_vec = {&d13, &d15, &d16};
  DayResult **days = day_vec.data();
  const Int_t n_days = (Int_t)day_vec.size();

  const TString results_dir = project_root + "/results";
  gSystem->mkdir(results_dir, kTRUE);
  const TString result_path = results_dir + "/ge_" + out_tag + ".result";
  std::ofstream out(result_path.Data());
  out << "# Ge-73m gamma energy results   tag=" << out_tag << std::endl;
  out << "# method  ge_mu  fit_err  cal_err  bkg_err  chi2  gain_err  cal_off  "
         "sys_err  anchor_err  label"
      << std::endl;
  out << std::fixed << std::setprecision(6);

  std::cout << std::endl;
  std::cout << "Method 1 -- in-situ simultaneous fit (Ge Peak mu, post-cal):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->ge_labels.size(); i++) {
      Float_t ce =
          i < days[d]->ge_cal_errs.size() ? days[d]->ge_cal_errs[i] : 0;
      Float_t be =
          i < days[d]->ge_bkg_errs.size() ? days[d]->ge_bkg_errs[i] : 0;
      Float_t gz =
          i < days[d]->ge_gain_errs.size() ? days[d]->ge_gain_errs[i] : 0;
      Float_t co =
          i < days[d]->ge_cal_offs.size() ? days[d]->ge_cal_offs[i] : 0;
      Float_t ae =
          i < days[d]->ge_anchor_errs.size() ? days[d]->ge_anchor_errs[i] : 0;
      out << "insitu  " << days[d]->ge_mus[i] << "  " << days[d]->ge_errs[i]
          << "  " << ce << "  " << be << "  " << days[d]->ge_chi2s[i] << "  "
          << gz << "  " << co << "  " << 0.0 << "  " << ae << "  "
          << SanitizeLabel(days[d]->ge_labels[i]) << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->ge_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->ge_mus[i]
                << " +/- " << days[d]->ge_errs[i] << " (fit) +/- " << ce
                << " (cal) keV (chi2/ndf = " << std::setprecision(3)
                << days[d]->ge_chi2s[i] << ")" << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "Method 2 -- live-time rate subtraction (Ge Peak mu):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->rs_labels.size(); i++) {
      Float_t ce =
          i < days[d]->rs_cal_errs.size() ? days[d]->rs_cal_errs[i] : 0;
      Float_t se =
          i < days[d]->rs_sys_errs.size() ? days[d]->rs_sys_errs[i] : 0;
      out << "ratesub  " << days[d]->rs_mus[i] << "  " << days[d]->rs_errs[i]
          << "  " << ce << "  " << 0.0 << "  " << days[d]->rs_chi2s[i] << "  "
          << (i < days[d]->rs_gain_errs.size() ? days[d]->rs_gain_errs[i]
                                               : 0.0f)
          << "  "
          << (i < days[d]->rs_cal_offs.size() ? days[d]->rs_cal_offs[i] : 0.0f)
          << "  " << se << "  " << 0.0 << "  "
          << SanitizeLabel(days[d]->rs_labels[i]) << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->rs_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->rs_mus[i]
                << " +/- " << days[d]->rs_errs[i] << " (fit) +/- " << ce
                << " (cal) +/- " << se << " (sys) keV" << std::endl;
    }
  }
  out.close();
  std::cout << std::endl;
  std::cout << "Wrote " << result_path << std::endl;
  std::cout << "Run CombineGeResult.cpp for the budget." << std::endl;
}
