// Rate-subtracted Ge-73m line: subtraction-scale scan and peak-shape test.
//
//   D(alpha) = S - alpha * B          (both in counts/s; alpha = 1 is nominal)
//
// Two fit models are run on every residual spectrum:
//
//   gaus : Gaussian + linear background            (what the pipeline does now)
//   tail : Gaussian + low-energy exponential tail  (FittingFunctions::LowTail,
//          the same peak shape the non-rate-sub measurements use) +
//          linear background -- the Am-241 case uses a flat background, but
//          the residual here keeps the linear term because the subtracted
//          background is not the calibration-peak background.
//
// Questions it answers:
//   1. does mu move when the S/B normalisation is wrong?  (alpha sloape)
//   2. does the low-energy tail fix the poor chi2 without moving mu?
//
// At alpha = 1 the `gaus` row reproduces the pipeline's RATESUB-FIT value.
//
// Usage:
//   root -b -q 'RateSubAlpha.cpp++'
//   root -b -q 'RateSubAlpha.cpp++("20260116")'
//   root -b -q 'RateSubAlpha.cpp++("all", 300)'   // bin width in eV (default
//   255)

#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "PlottingUtils.hpp"
#include <RtypesCore.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

const Float_t kLo = 62.0;
const Float_t kHi = 71.5;
const Double_t kPivot = 68.75;

struct DayCfg {
  const char *label;
  TString sig;
  TString bkg;
  Double_t p0;
  Double_t p1;
  const char *tag;
};

struct Point {
  Double_t alpha = 0;
  Double_t mu = 0;
  Double_t mu_err = 0;
  Double_t sigma = 0;
  Double_t bkg0 = 0;
  Double_t slope = 0;
  Double_t chi2 = -1;
  Double_t area = 0;
  Double_t tail_amp = -1;
  Double_t tail_ratio = -1;
  std::vector<Double_t> pars; // fitted parameters, for rebuilding the model
  Bool_t ok = kFALSE;
};

Float_t LiveTime(TFile *file, Int_t crystal) {
  TParameter<Float_t> *p = static_cast<TParameter<Float_t> *>(
      file->Get(Form("LiveTime_Filtered_Crystal%d_s", crystal)));
  return p ? p->GetVal() : 0.f;
}

// Same construction as CalibrateAndFit.cpp:BuildRateHist, except the input path
// is resolved from the macro's own location so this runs from macros/ or from
// the project root alike.
TH1F *RateHist(const TString &run, TF1 *cal, const TString &name,
               Double_t bin_w) {
  const TString path =
      Paths::ProjectRootOf(__FILE__) + "/root_files/filtered/" + run + ".root";
  TFile *file = TFile::Open(path, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    delete file;
    return nullptr;
  }
  Int_t nbins = (Int_t)std::lround((kHi - kLo) / bin_w);
  TH1F *sum = new TH1F(name, ";Deposited Energy [keV];Rate [counts/s]", nbins,
                       kLo, kHi);
  sum->SetDirectory(nullptr);
  sum->Sumw2();

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t lt = LiveTime(file, c);
    if (lt <= 0)
      continue;
    TH1F *ch = new TH1F(name + Form("_c%d", c), "", nbins, kLo, kHi);
    ch->SetDirectory(nullptr);
    ch->Sumw2();
    Float_t energy = 0;
    tree->SetBranchAddress("energykeV", &energy);
    Int_t n = (Int_t)tree->GetEntries();
    for (Int_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      ch->Fill(cal->Eval(energy));
    }
    ch->Scale(1.0 / lt);
    sum->Add(ch);
    delete ch;
  }
  file->Close();
  delete file;
  return sum;
}

// Gaussian + linear background. par: mu, sigma, amp, bkg0, slope.
Double_t ModelGaus(Double_t *x, Double_t *p) {
  Double_t g[3] = {p[0], p[1], p[2]};
  return FittingFunctions::Gaussian(x, g) + p[3] + p[4] * (x[0] - kPivot);
}

// Gaussian + low-energy exponential tail + linear background.
// par: mu, sigma, amp, tail_amp, tail_ratio, lin_tail_amp, lin_tail_slope,
//      bkg0, slope.
Double_t ModelTail(Double_t *x, Double_t *p) {
  Double_t g[3] = {p[0], p[1], p[2]};
  Double_t t[6] = {p[0], p[1], p[3], p[4], p[5], p[6]};
  return FittingFunctions::Gaussian(x, g) + FittingFunctions::LowTail(x, t) +
         p[7] + p[8] * (x[0] - kPivot);
}

Point FitResidual(TH1F *res, Double_t alpha, Bool_t use_tail) {
  static Int_t counter = 0;
  Point out;
  out.alpha = alpha;

  Double_t edge_lo = res->Integral(res->FindBin(kLo + 0.2), res->FindBin(66.0));
  Double_t nedge = res->FindBin(66.0) - res->FindBin(kLo + 0.2) + 1;
  Double_t bkg0 = (nedge > 0) ? edge_lo / nedge : 0.0;
  Int_t pk_bin = res->GetMaximumBin();
  Double_t pk_val = res->GetBinContent(pk_bin) - bkg0;
  Double_t sigma0 = 0.85;

  TF1 *model = nullptr;
  if (use_tail) {
    model = new TF1(Form("ratesub_tail_%d", counter++), ModelTail, kLo, kHi, 9);
    model->SetParNames("Mu", "Sigma", "GausAmp", "LowExpAmp", "LowExpRatio",
                       "LowLinAmp", "LowLinSlope", "Bkg0", "BkgSlope");
    model->SetParameters(68.75, sigma0, pk_val * sigma0 * 2.5066,
                         pk_val * sigma0 * 0.1, 5.0, 0.0, 0.0, bkg0, 0.0);
    model->SetParLimits(0, kLo + 1.0, kHi - 1.0);
    model->SetParLimits(1, 0.2, 3.0);
    model->SetParLimits(2, 0.0, pk_val * sigma0 * 2.5066 * 5.0);
    model->SetParLimits(3, 0.0, pk_val * sigma0 * 2.5066);
    model->SetParLimits(4, 0.5, 8.0);
  } else {
    model = new TF1(Form("ratesub_gaus_%d", counter++), ModelGaus, kLo, kHi, 5);
    model->SetParNames("Mu", "Sigma", "Area", "Bkg0", "BkgSlope");
    model->SetParameters(68.75, sigma0, pk_val * sigma0 * 2.5066, bkg0, 0.0);
    model->SetParLimits(0, kLo + 1.0, kHi - 1.0);
    model->SetParLimits(1, 0.2, 3.0);
  }

  TFitResultPtr fr = res->Fit(model, "S Q R N");
  Bool_t ok = fr.Get() && fr->IsValid() && model->GetNDF() > 0;
  if (ok) {
    out.mu = model->GetParameter(0);
    out.mu_err = model->GetParError(0);
    out.sigma = model->GetParameter(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.area = model->GetParameter(2);
    if (use_tail) {
      out.tail_amp = model->GetParameter(3);
      out.tail_ratio = model->GetParameter(4);
      out.bkg0 = model->GetParameter(7);
      out.slope = model->GetParameter(8);
    } else {
      out.bkg0 = model->GetParameter(3);
      out.slope = model->GetParameter(4);
    }
    out.ok = kTRUE;
    out.pars.resize(model->GetNpar());
    for (Int_t k = 0; k < model->GetNpar(); k++)
      out.pars[k] = model->GetParameter(k);
  }
  delete model;
  return out;
}

// Per-bin (data - model)/sigma at the fitted parameters, plus the chi2 split
// into low / core / high regions. A smooth excess across the window points at
// the background shape; a few large pulls next to the peak point at a
// contaminant line or a mis-modelled shoulder.
void PrintPullMap(TH1F *res, const Point &p, const char *label,
                  const TString &out_path) {
  if (!p.ok || p.pars.empty())
    return;
  std::ofstream out(out_path.Data());
  if (!out) {
    std::cerr << "WARNING: cannot write pull map to " << out_path << std::endl;
    return;
  }
  TF1 model("pull_map_model", p.pars.size() == 9 ? ModelTail : ModelGaus, kLo,
            kHi, (Int_t)p.pars.size());
  for (size_t k = 0; k < p.pars.size(); k++)
    model.SetParameter((Int_t)k, p.pars[k]);

  out << "pull map, " << label << " model, alpha = 1" << std::endl;
  out << std::setw(10) << "E[keV]" << std::setw(12) << "data" << std::setw(12)
      << "model" << std::setw(9) << "pull" << std::setw(9) << "err"
      << "  bar" << std::endl;

  Int_t b_lo = res->FindBin(kLo + 1e-3);
  Int_t b_hi = res->FindBin(kHi - 1e-3);
  Double_t chi2_low = 0, chi2_core = 0, chi2_high = 0;
  Int_t n_low = 0, n_core = 0, n_high = 0;
  Double_t worst_pull = 0, worst_x = 0;
  std::vector<Double_t> pull_x, pull_y;
  static Int_t ks_counter = 0;
  for (Int_t b = b_lo; b <= b_hi; b++) {
    Double_t x = res->GetBinCenter(b);
    Double_t y = res->GetBinContent(b);
    Double_t e = res->GetBinError(b);
    if (e <= 0)
      continue;
    Double_t m = model.Eval(x);
    Double_t pull = (y - m) / e;

    Int_t nbar = (Int_t)std::lround(std::fabs(pull) * 4);
    if (nbar > 60)
      nbar = 60;
    TString bar;
    for (Int_t k = 0; k < nbar; k++)
      bar += (pull > 0 ? "+" : "-");

    out << std::setw(10) << std::fixed << std::setprecision(3) << x
        << std::setw(12) << std::setprecision(4) << y << std::setw(12)
        << std::setprecision(4) << m << std::setw(9) << std::setprecision(2)
        << pull << std::setw(9) << std::setprecision(4) << e << "  " << bar
        << std::endl;

    Double_t c2 = pull * pull;
    pull_x.push_back(x);
    pull_y.push_back(pull);
    if (x < 66.5) {
      chi2_low += c2;
      n_low++;
    } else if (x < 70.5) {
      chi2_core += c2;
      n_core++;
    } else {
      chi2_high += c2;
      n_high++;
    }
    if (std::fabs(pull) > std::fabs(worst_pull)) {
      worst_pull = pull;
      worst_x = x;
    }
  }

  // Pull histogram + KS test, mirroring PlottingUtils::PlotPullHistogram
  // (private): 82 bins over [-5.5, 5.5], unit-Gaussian reference scaled to the
  // number of points. The figure goes next to the rate_residual_* plots.
  Double_t ks_p = -1;
  if (!pull_y.empty()) {
    TH1D *ph = new TH1D(Form("ks_pull_%d", ++ks_counter), ";#delta/#sigma;bins",
                        82, -5.5, 5.5);
    ph->SetDirectory(nullptr);
    for (size_t i = 0; i < pull_y.size(); i++)
      ph->Fill(pull_y[i]);
    TH1D *ref = new TH1D(Form("ks_ref_%d", ks_counter), "", 82, -5.5, 5.5);
    ref->SetDirectory(nullptr);
    for (Int_t i = 1; i <= ref->GetNbinsX(); i++) {
      Double_t lo = ref->GetBinLowEdge(i);
      Double_t hi = lo + ref->GetBinWidth(i);
      ref->SetBinContent(i, (Double_t)pull_y.size() *
                                (TMath::Freq(hi) - TMath::Freq(lo)));
    }
    ks_p = ph->KolmogorovTest(ref);

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(ph, kAzure);
    PlottingUtils::AddText(TString::Format("KS p = %.3f", ks_p), 0.85, 0.85);
    PlottingUtils::SaveFigure(canvas, TString("pull_hist_") + label,
                              "ratesub_alpha", PlotSaveOptions::kLINEAR);
    gROOT->GetListOfCanvases()->Remove(canvas);
    canvas->SetBatch(kTRUE);
    delete canvas;
    delete ref;
    delete ph;
  }

  out << std::endl
      << "KS p-value (pulls vs unit Gaussian) = " << std::fixed
      << std::setprecision(3) << ks_p << std::endl;
  // Pull sequence statistics: rms tells you whether the errors are scaled
  // right, lag-1 and Durbin-Watson whether adjacent bins are correlated
  // (bin-boundary / quantisation structure) rather than random.
  Double_t pmean = 0, prms = 0, plag1 = 0, pdw = 0;
  if (pull_y.size() > 1) {
    for (size_t i = 0; i < pull_y.size(); i++)
      pmean += pull_y[i];
    pmean /= (Double_t)pull_y.size();
    Double_t den = 0, num = 0, dsum = 0;
    for (size_t i = 0; i < pull_y.size(); i++) {
      Double_t d = pull_y[i] - pmean;
      den += d * d;
      if (i + 1 < pull_y.size()) {
        num += d * (pull_y[i + 1] - pmean);
        Double_t dd = pull_y[i + 1] - pull_y[i];
        dsum += dd * dd;
      }
    }
    if (den > 0) {
      prms = std::sqrt(den / (Double_t)pull_y.size());
      plag1 = num / den;
      pdw = dsum / den;
    }
  }
  out << std::endl
      << "pull sequence: mean " << std::fixed << std::setprecision(3) << pmean
      << "  rms " << prms << "  lag-1 " << plag1 << "  DW " << pdw << std::endl;

  out << std::endl << "chi2 by region:" << std::endl;
  out << "  62.0 - 66.5 : " << std::fixed << std::setprecision(1) << chi2_low
      << " / " << n_low << " bins" << std::endl;
  out << "  66.5 - 70.5 : " << chi2_core << " / " << n_core << " bins  (peak)"
      << std::endl;
  out << "  70.5 - 71.5 : " << chi2_high << " / " << n_high << " bins"
      << std::endl;
  out << "worst pull " << worst_pull << " at " << std::setprecision(3)
      << worst_x << " keV" << std::endl;
  out.close();

  std::cout << "  pull map -> " << out_path << std::endl;
  std::cout << "    chi2/ndf by region  low " << std::fixed
            << std::setprecision(1) << chi2_low << "/" << n_low << "  core "
            << chi2_core << "/" << n_core << "  high " << chi2_high << "/"
            << n_high << "   KS p " << std::setprecision(3) << ks_p
            << "   worst pull " << worst_pull << " at " << std::setprecision(3)
            << worst_x << " keV" << std::endl;
  std::cout << "    pull rms " << std::setprecision(3) << prms << "  lag-1 "
            << plag1 << "  DW " << pdw << std::endl;
}

// Least-squares slope of mu vs alpha, returned in eV per unit alpha.
Double_t SlopePerAlpha(const std::vector<Point> &v) {
  Double_t sx = 0, sy = 0, sxx = 0, sxy = 0;
  Int_t n = 0;
  for (size_t i = 0; i < v.size(); i++) {
    if (!v[i].ok)
      continue;
    sx += v[i].alpha;
    sy += v[i].mu;
    sxx += v[i].alpha * v[i].alpha;
    sxy += v[i].alpha * v[i].mu;
    n++;
  }
  if (n < 2)
    return 0;
  Double_t det = n * sxx - sx * sx;
  if (std::fabs(det) < 1e-12)
    return 0;
  return (n * sxy - sx * sy) / det * 1000.0; // keV/alpha -> eV/alpha
}

Bool_t FitFailed(const std::vector<Point> &v) {
  for (size_t i = 0; i < v.size(); i++)
    if (!v[i].ok)
      return kTRUE;
  return kFALSE;
}

void ScanDay(const DayCfg &cfg, Double_t bin_w) {
  TF1 cal("cal_alpha", "pol1", 0, 200);
  cal.SetParameters(cfg.p0, cfg.p1);

  TH1F *sig = RateHist(cfg.sig, &cal, "alpha_sig", bin_w);
  TH1F *bkg = RateHist(cfg.bkg, &cal, "alpha_bkg", bin_w);
  if (!sig || !bkg) {
    delete sig;
    delete bkg;
    return;
  }

  std::cout << std::endl << "===== " << cfg.label << " =====" << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "  S = " << cfg.sig << std::endl
            << "  B = " << cfg.bkg << std::endl
            << "  cal = " << cfg.p0 << " + " << cfg.p1 << "*raw" << std::endl;

  std::cout << std::endl
            << "  " << std::setw(6) << "alpha" << std::setw(11) << "mu_gaus"
            << std::setw(8) << "chi2" << std::setw(11) << "mu_tail"
            << std::setw(8) << "chi2" << std::setw(10) << "tail/sig"
            << std::setw(9) << "mu_diff" << std::endl;

  std::vector<Point> gaus_pts, tail_pts;
  Point gaus1;
  static Int_t pull_counter = 0;
  for (Int_t k = -4; k <= 4; k++) {
    Double_t alpha = 1.0 + 0.01 * k;
    TH1F *res = static_cast<TH1F *>(sig->Clone("alpha_res"));
    res->SetDirectory(nullptr);
    res->Add(bkg, -alpha);

    Point gp = FitResidual(res, alpha, kFALSE);
    Point tp = FitResidual(res, alpha, kTRUE);
    gaus_pts.push_back(gp);
    tail_pts.push_back(tp);
    if (std::fabs(alpha - 1.0) < 1e-9)
      gaus1 = gp;

    std::cout << "  " << std::setw(6) << std::setprecision(2) << alpha
              << std::setw(11) << (gp.ok ? Form("%.5f", gp.mu) : "fail")
              << std::setw(8) << (gp.ok ? Form("%.2f", gp.chi2) : "-")
              << std::setw(11) << (tp.ok ? Form("%.5f", tp.mu) : "fail")
              << std::setw(8) << (tp.ok ? Form("%.2f", tp.chi2) : "-")
              << std::setw(10) << (tp.ok ? Form("%.2f", tp.tail_ratio) : "-")
              << std::setw(9)
              << ((gp.ok && tp.ok) ? Form("%.1f", (tp.mu - gp.mu) * 1000) : "-")
              << std::endl;
    delete res;
  }

  // Values at alpha = 1.
  Bool_t have_g = kFALSE, have_t = kFALSE;
  Double_t gmu = 0, gerr = 0, gchi = -1, tmu = 0, terr = 0, tchi = -1;
  Double_t tramp = 0, tratio = 0;
  for (size_t i = 0; i < gaus_pts.size(); i++) {
    if (std::fabs(gaus_pts[i].alpha - 1.0) > 1e-9)
      continue;
    if (gaus_pts[i].ok) {
      gmu = gaus_pts[i].mu;
      gerr = gaus_pts[i].mu_err * 1000;
      gchi = gaus_pts[i].chi2;
      have_g = kTRUE;
    }
    if (tail_pts[i].ok) {
      tmu = tail_pts[i].mu;
      terr = tail_pts[i].mu_err * 1000;
      tchi = tail_pts[i].chi2;
      tramp = tail_pts[i].tail_amp;
      tratio = tail_pts[i].tail_ratio;
      have_t = kTRUE;
    }
  }

  std::cout << std::endl << "  at alpha = 1" << std::endl;
  if (have_g)
    std::cout << "    gaus : mu = " << std::setprecision(5) << gmu << " +/- "
              << std::setprecision(1) << gerr << " eV   chi2/ndf "
              << std::setprecision(2) << gchi << std::endl;
  else
    std::cout << "    gaus : fit failed" << std::endl;
  if (have_t) {
    std::cout << "    tail : mu = " << std::setprecision(5) << tmu << " +/- "
              << std::setprecision(1) << terr << " eV   chi2/ndf "
              << std::setprecision(2) << tchi << "   tail amp "
              << std::setprecision(3) << tramp << "   tail/sigma "
              << std::setprecision(2) << tratio << std::endl;
    if (have_g)
      std::cout << "    mu(tail) - mu(gaus) = " << std::setprecision(1)
                << (tmu - gmu) * 1000 << " eV" << std::endl;
  } else {
    std::cout << "    tail : fit failed" << std::endl;
  }

  Double_t gslope = SlopePerAlpha(gaus_pts);
  std::cout << std::endl
            << "  dmu/dalpha  gaus = " << std::setprecision(1) << gslope
            << " eV/unit (" << std::setprecision(2) << gslope / 100.0
            << " eV per 1%)" << std::endl;
  if (FitFailed(tail_pts)) {
    std::cout << "  dmu/dalpha  tail = not computed ("
              << "tail model did not fit at every alpha)" << std::endl;
  } else {
    Double_t tslope = SlopePerAlpha(tail_pts);
    std::cout << "  dmu/dalpha  tail = " << std::setprecision(1) << tslope
              << " eV/unit (" << std::setprecision(2) << tslope / 100.0
              << " eV per 1%)" << std::endl;
    if (std::fabs(tslope) > 1e-9)
      std::cout << "  normalisation error needed for an 80 eV move (tail) = "
                << std::setprecision(0) << 80.0 / std::fabs(tslope) * 100.0
                << " %" << std::endl;
  }
  if (FitFailed(gaus_pts))
    std::cout << "  WARNING: the Gaussian model failed at least one alpha"
              << std::endl;

  if (have_g) {
    TH1F *res1 =
        static_cast<TH1F *>(sig->Clone(Form("pullres_%d", ++pull_counter)));
    res1->SetDirectory(nullptr);
    res1->Add(bkg, -1.0);
    TString pulldir = Paths::ProjectRootOf(__FILE__) + "/plots/ratesub_alpha";
    gSystem->mkdir(pulldir, kTRUE);
    TString lbl = Form("%s_bw%d", cfg.tag, (Int_t)std::lround(bin_w * 1000));
    TString pullpath = pulldir + "/pull_map_" + lbl + ".log";
    PrintPullMap(res1, gaus1, lbl, pullpath);
    delete res1;
  }

  delete sig;
  delete bkg;
}

} // namespace

void RateSubAlpha(const char *which = "all", Double_t bin_w_ev = 255.0) {
  const TString want(which);
  const Double_t bin_w = bin_w_ev / 1000.0;
  gROOT->SetBatch(kTRUE);
  PlottingUtils::SetStylePreferences(PlotSaveFormat::kPNG);
  PlottingUtils::SetPlotsBaseDir(Paths::ProjectRootOf(__FILE__) + "/plots");

  std::vector<DayCfg> days;
  // Calibration constants are the ones the pipeline prints as
  // "Cal pol1 <day>: p0 = ..., p1 = ...".
  days.push_back({"01/14 CuShield 10%",
                  Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                  Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, 0.501012,
                  0.987896, "20260114"});
  days.push_back({"01/15 NoShield 5%",
                  Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                  Constants::NOSHIELDBACKGROUND_5PERCENT_20260115, 0.471026,
                  0.990649, "20260115"});
  days.push_back(
      {"01/16 Graphite 10%",
       Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
       Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
       0.947697, 0.984380, "20260116"});

  Bool_t any = kFALSE;
  for (size_t i = 0; i < days.size(); i++) {
    // Accept "all", a tag ("20260115") or a label fragment ("01/15",
    // "NoShield", "Graphite"). Contains() is a literal substring test, so check
    // both ways.
    const TString tag(days[i].tag);
    if (want != "all" && !want.Contains(days[i].label) && !tag.Contains(want) &&
        !want.Contains(tag))
      continue;
    any = kTRUE;
    ScanDay(days[i], bin_w);
  }
  if (!any)
    std::cout << "no day matched \"" << which
              << "\" -- use all, a tag (e.g. 20260115) or a label fragment "
                 "(01/15, NoShield, Graphite)"
              << std::endl;
}
