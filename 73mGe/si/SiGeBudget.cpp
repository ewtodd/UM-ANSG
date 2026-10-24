// The 68.75 keV 72Ge(n,g) line in the PSI Amptek Si, from the processed
// waveforms (ProcessWaveforms.cpp), in either energy of the same records:
// the trapezoid ("trapezoid") or the digitiser's own on-board energy
// ("digitizer"). Done the collaborator's way: a Gaussian on a linear
// background over 66.5-71 keV, calibrated with Am-241 59.5409 keV (source
// runs) and the beam runs' own Pb Ka2/Ka1 (cal A, the result). Alternatives:
// Am-241 and Ba-133 53.1622, 79.6142, 80.9979 keV (cal B, source lines
// only) and all six lines (cal C).
// Uncertainties: the Ge centroid fit, the calibration offset and slope (the
// pivot at the calibration points' weighted centroid, where the two are
// uncorrelated), the background (|linear - flat|) and the fit window (RMS
// over a window scan). Checks: per run, per day, the width fixed, the W Kb
// group from the collimator, the 53.44 keV line, calibration variants,
// isolated events, and the run-to-run stability of the line positions.
//
// Writes results/si_ge_<energy>.result, the spectra to
// /labdata/ANSG/PSI/processed/Si/ge_spectra_<energy>.root, and the figures
// to plots/studies/psi_si/ge_<energy>. "compare" (after both) draws the two
// energies' Ge spectra together; "bootstrap" measures how far two estimates
// of the same records scatter against each other (Bootstrap()).
//
// Usage: root -b -q 'SiGeBudget.cpp+("trapezoid")'
//        root -b -q 'SiGeBudget.cpp+("digitizer")'
//        root -b -q 'SiGeBudget.cpp+("compare")'
//        root -b -q 'SiGeBudget.cpp+("bootstrap", <replicas>)'
#include "SiAnalysis.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TVectorD.h>

#include <fstream>
#include <map>

namespace SiGe {

using namespace Si;

const char *const SPECTRA_DIR = "/labdata/ANSG/PSI/processed/Si";
const Double_t WIN_LO = 66.5, WIN_HI = 71.0;
const Double_t E_GE_GUESS = 68.75;

// Reference energies' uncertainties (keV): Am-241 and Ba-133 from
// ENSDF/DDEP, Pb K X-rays from Deslattes et al. (2003).
const Double_t BA133_53 = 53.1622;
const Double_t AM241_ERR = 0.0001, BA133_53_ERR = 0.0006, BA133_79_ERR = 0.0019,
               BA133_81_ERR = 0.0011, PB_KA_ERR = 0.0009;
// W K X-rays from the Amptek's collimator.
const Double_t W_KA2 = 57.9817, W_KA1 = 59.3182;
const Double_t W_KB3 = 66.9514, W_KB1 = 67.2443, W_KB2 = 69.067;
// The 53.44 keV line of 73mGe, as a check of the scale below Am.
const Double_t GE73M_53 = 53.440;

// ------------------------------------------------------------------ runs

struct GeRun {
  TString run, label;
  Int_t day;
  // One of the named Ge beam runs; kFALSE for BeamTurnOn_4.
  Bool_t main;
};

std::vector<GeRun> GeRuns() {
  std::vector<GeRun> runs;
  runs.push_back(
      {"PSI_Beam_GeTarget_20251010_1525_1", "GeTarget_1", 10, kTRUE});
  runs.push_back(
      {"PSI_Beam_GeTarget_20251010_1525_2", "GeTarget_2", 10, kTRUE});
  runs.push_back({"PSI_Beam_CdTe_NoGe_Si_YesGe_20251011_1056_1", "SiYesGe_1056",
                  11, kTRUE});
  runs.push_back({"PSI_Beam_Ge_20251011_1715_1", "Ge_1715_1", 11, kTRUE});
  runs.push_back({"PSI_Beam_Ge_20251011_1715_2", "Ge_1715_2", 11, kTRUE});
  runs.push_back({"PSI_Beam_LongWf_CdTe_500V_NoGe_Si_YesGe_20251012_1001_1_1",
                  "LongWf500V_1001", 12, kTRUE});
  runs.push_back({"PSI_Beam_LongWf_CdTe_NoGe_Si_YesGe_20251012_1517_1",
                  "LongWf_1517", 12, kTRUE});
  runs.push_back(
      {"PSI_Running_BeamTurnOn_20251009_1914_4", "BeamTurnOn_4", 9, kFALSE});
  return runs;
}

Int_t DayOf(const TString &run) {
  Int_t at = run.Index("_202510");
  if (at < 0)
    return 0;
  return TString(run(at + 7, 2)).Atoi();
}

std::vector<TString> FilterDay(const char *const *names, Int_t n, Int_t day) {
  std::vector<TString> runs;
  for (Int_t i = 0; i < n; i++)
    if (day == 0 || DayOf(names[i]) == day)
      runs.push_back(TString("PSI_") + names[i]);
  return runs;
}

// Every source run with the Am-241 line on the Si at the final gain setting.
// The "Si_Sm" run shows the Am line of the CdTe's Am source, at the same
// position, so it is kept.
std::vector<TString> AmRuns(Int_t day) {
  const Int_t n = 13;
  const char *const names[n] = {
      "Calibration_CdTe_Am_Si_Am_NoGe_20251010_1014_2",
      "Calibration_CdTe_Am_Si_Am_20251010_1422_1",
      "Calibration_CdTe_Am_Si_Am_20251010_1422_2",
      "Calib_CdTe_Am_Si_Am_20251011_0920_1",
      "Calib_CdTe_ba_Si_am_20251011_1043_1",
      "Calib_CdTe_ba_Si_am_20251011_1043_2",
      "Calibration_CdTe_Am_Si_Am_20251011_1617_1",
      "Calibration_CdTe_Co_Si_Am_20251011_1617_1",
      "Calib_Long_CdTe_500V_Am_Si_Am_20251012_0917_1",
      "Calib_Long_CdTe_500V_Ba_Si_Am_20251012_0917_1",
      "Calib_Long_CdTe_500V_Co_Si_Am_20251012_0917_1",
      "Calib_Long_CdTe_Am_Si_Sm_20251012_0917_1",
      "Calibration_CdTe_Am_Si_Am_20251013_0652_1"};
  return FilterDay(names, n, day);
}

std::vector<TString> BaRuns(Int_t day) {
  const Int_t n = 13;
  const char *const names[n] = {
      "Calibration_CdTe_Co_Si_Ba_NoGe_20251010_1032_1",
      "Calibration_CdTe_Co_Si_Ba_NoGe_20251010_1056_1",
      "Calibration_CdTe_Co_Si_Ba_20251010_1449_1",
      "Calib_CdTe_co_Si_ba_20251011_0938_1",
      "Calibration_CdTe_Co_Si_Ba_20251011_1617_1",
      "Calibration_CdTe_none_Si_Ba_20251011_1617_1",
      "Calib_Long_CdTe_Cs_Si_Ba_20251012_0917_1",
      "Calib_Long_CdTe_500V_Cs_Si_Ba_20251012_0917_1",
      "Calib_Long_CdTe_500V_Cs_Si_Ba_20251012_0917_2",
      "Calib_Long_CdTe_500V_Cs_Si_Ba_20251012_0917_3",
      "Calib_Long_CdTe_500V_Cs_Si_Ba_20251012_0917_4",
      "Calibration_CdTe_Co_Si_Ba_20251013_0652_1",
      "Calibration_CdTe_Cs_Si_Ba_20251013_0652_1"};
  return FilterDay(names, n, day);
}

// ------------------------------------------------------- calibration lines

struct Point {
  TString name;
  Double_t ch = 0, ch_err = 0, e = 0, e_err = 0, sigma = 0, area = 0;
};

// Variations of the calibration-line fits: channels added to each side of
// the default window, and the background order.
struct LineSetting {
  Double_t widen = 0;
  Int_t background = 1;
};

Point MakePoint(const TString &name, const FitOutput &f, Int_t k, Double_t e,
                Double_t e_err) {
  Point p;
  p.name = name;
  p.ch = f.mean[k];
  p.ch_err = f.mean_err[k];
  p.e = e;
  p.e_err = e_err;
  p.sigma = f.sigma[0];
  p.area = f.area[k];
  return p;
}

Point FitAm(TH1D *h, const LineSetting &s) {
  std::vector<LineSpec> one(1);
  one[0].mean_guess = 1428;
  FitOutput f =
      FitLines(h, 1395 - s.widen, 1465 + s.widen, one, 1, s.background);
  return MakePoint("Am 59.54", f, 0, ReferenceLines::AM241, AM241_ERR);
}

Point FitBa53(TH1D *h, const LineSetting &s) {
  std::vector<LineSpec> one(1);
  one[0].mean_guess = 1275;
  FitOutput f =
      FitLines(h, 1245 - s.widen, 1310 + s.widen, one, 1, s.background);
  return MakePoint("Ba 53.16", f, 0, BA133_53, BA133_53_ERR);
}

void FitBa8081(TH1D *h, const LineSetting &s, Point &p80, Point &p81) {
  std::vector<LineSpec> two(2);
  two[0].mean_guess = 1908;
  two[1].mean_guess = 1941;
  FitOutput f =
      FitLines(h, 1870 - s.widen, 1985 + s.widen, two, 1, s.background);
  p80 = MakePoint("Ba 79.61", f, 0, ReferenceLines::BA133_79, BA133_79_ERR);
  p81 = MakePoint("Ba 81.00", f, 1, ReferenceLines::BA133_81, BA133_81_ERR);
}

void FitPb(TH1D *h, const LineSetting &s, Point &ka2, Point &ka1) {
  std::vector<LineSpec> two(2);
  two[0].mean_guess = 1745;
  two[1].mean_guess = 1797;
  FitOutput f =
      FitLines(h, 1710 - s.widen, 1840 + s.widen, two, 1, s.background, 7.0);
  ka2 = MakePoint("Pb Ka2", f, 0, ReferenceLines::PB_KA2, PB_KA_ERR);
  ka1 = MakePoint("Pb Ka1", f, 1, ReferenceLines::PB_KA1, PB_KA_ERR);
}

Calibration Calibrate(const std::vector<Point> &points) {
  std::vector<Double_t> ch, ch_err, e, e_err;
  for (size_t i = 0; i < points.size(); i++) {
    ch.push_back(points[i].ch);
    ch_err.push_back(points[i].ch_err);
    e.push_back(points[i].e);
    e_err.push_back(points[i].e_err);
  }
  return FitCalibration(ch, ch_err, e, e_err, -1);
}

// ----------------------------------------------------------------- Ge fits

struct GeFit {
  Bool_t ok = kFALSE;
  Double_t mu = 0, mu_err = 0, sigma = 0, sigma_err = 0, area = 0, area_err = 0;
  Double_t chi2 = 0;
  Int_t ndf = 0;
  Double_t lo = 0, hi = 0;
  Double_t c0 = 0, c1 = 0, c2 = 0;
  // The W Kb group, when fitted.
  Double_t w_area = 0, w_area_err = 0;
  Int_t status = -1;
};

GeFit FromOutput(const FitOutput &f) {
  GeFit g;
  g.status = f.status;
  g.ok =
      f.sigma[0] > 1.05 && f.sigma[0] < 39.0 && f.area[0] > 2.0 * f.area_err[0];
  g.mu = f.mean[0];
  g.mu_err = f.mean_err[0];
  g.sigma = f.sigma[0];
  g.sigma_err = f.sigma_err[0];
  g.area = f.area[0];
  g.area_err = f.area_err[0];
  g.chi2 = f.chi2;
  g.ndf = f.ndf;
  g.lo = f.lo;
  g.hi = f.hi;
  g.c0 = f.c0;
  g.c1 = f.c1;
  g.c2 = f.c2;
  return g;
}

// One Gaussian on a background of order background over [lo_kev, hi_kev],
// the window turned into channels with the given calibration.
GeFit FitGe(TH1D *h, const Calibration &cal, Double_t lo_kev, Double_t hi_kev,
            Int_t background, Double_t fixed_sigma = -1) {
  std::vector<LineSpec> one(1);
  one[0].mean_guess = cal.Channel(E_GE_GUESS);
  one[0].mean_range = 12;
  Bool_t fix = fixed_sigma > 0;
  FitOutput f = FitLines(h, cal.Channel(lo_kev), cal.Channel(hi_kev), one, 1,
                         background, fix ? fixed_sigma : 6.2, fix);
  return FromOutput(f);
}

// Ge plus the W Kb3/Kb1/Kb2 group from the Amptek collimator, the W lines at
// their calibrated positions with one free (signed) amplitude and fixed
// relative intensities, all sharing the Ge width.
const Double_t W_RATIO_KB3 = 0.52, W_RATIO_KB2 = 0.34;
struct GeWModel {
  Double_t x_ref;
  Double_t position[3];
  Double_t ratio[3];
  Double_t operator()(Double_t *x, Double_t *p) {
    Double_t v = p[0] + p[1] * (x[0] - x_ref);
    Double_t norm = 1.0 / (TMath::Sqrt(2.0 * TMath::Pi()) * p[2]);
    Double_t z = (x[0] - p[4]) / p[2];
    v += p[3] * norm * TMath::Exp(-0.5 * z * z);
    for (Int_t k = 0; k < 3; k++) {
      Double_t zw = (x[0] - position[k]) / p[2];
      v += p[5] * ratio[k] * norm * TMath::Exp(-0.5 * zw * zw);
    }
    return v;
  }
};

GeFit FitGeWithW(TH1D *h, const Calibration &cal, Double_t lo_kev,
                 Double_t hi_kev) {
  GeFit seed = FitGe(h, cal, lo_kev, hi_kev, 1);
  GeWModel model;
  Double_t lo = cal.Channel(lo_kev), hi = cal.Channel(hi_kev);
  model.x_ref = 0.5 * (lo + hi);
  model.position[0] = cal.Channel(W_KB3);
  model.position[1] = cal.Channel(W_KB1);
  model.position[2] = cal.Channel(W_KB2);
  model.ratio[0] = W_RATIO_KB3;
  model.ratio[1] = 1.0;
  model.ratio[2] = W_RATIO_KB2;
  TF1 f(Form("gew_%d", FitCounter()++), model, lo, hi, 6);
  f.SetParameters(seed.c0, seed.c1, seed.sigma, seed.area, seed.mu, 0.0);
  f.SetParLimits(2, 1.0, 40.0);
  f.SetParLimits(4, seed.mu - 12, seed.mu + 12);
  h->Fit(&f, "QLRN0");
  TFitResultPtr r = h->Fit(&f, "QLRN0SM");
  GeFit g;
  g.ok = r.Get() && r->Status() == 0;
  g.mu = f.GetParameter(4);
  g.mu_err = f.GetParError(4);
  g.sigma = f.GetParameter(2);
  g.area = f.GetParameter(3);
  g.area_err = f.GetParError(3);
  g.w_area = f.GetParameter(5);
  g.w_area_err = f.GetParError(5);
  return g;
}

// ------------------------------------------------------- energy breakdown

struct GeEnergy {
  Double_t e = 0, fit = 0, off = 0, slope = 0, cal = 0, background = 0,
           window = 0, total = 0;
};

GeEnergy Energy(const GeFit &lin, const GeFit &flat, const Calibration &cal,
                Double_t window_err) {
  GeEnergy r;
  r.e = cal.Energy(lin.mu);
  r.fit = TMath::Abs(cal.c1) * lin.mu_err;
  r.off = cal.OffsetError();
  r.slope = cal.SlopeError(lin.mu);
  r.cal = cal.Error(lin.mu);
  r.background = TMath::Abs(cal.Energy(lin.mu) - cal.Energy(flat.mu));
  r.window = window_err;
  r.total = TMath::Sqrt(r.fit * r.fit + r.cal * r.cal +
                        r.background * r.background + r.window * r.window);
  return r;
}

void PrintEnergy(const TString &label, const GeEnergy &r, const GeFit &lin,
                 const Calibration &cal) {
  std::cout << "SIGE " << std::left << std::setw(46) << label << std::right
            << std::fixed << std::setprecision(4) << r.e << " keV  +- "
            << std::setprecision(1) << 1000 * r.fit << " fit +- "
            << 1000 * r.off << " off +- " << 1000 * r.slope << " slope +- "
            << 1000 * r.background << " bkg +- " << 1000 * r.window
            << " win = " << 1000 * r.total << " eV  chi2/ndf "
            << std::setprecision(1) << lin.chi2 << "/" << lin.ndf << "  sigma "
            << std::setprecision(0) << 1000 * cal.c1 * lin.sigma << " eV  area "
            << lin.area << " +- " << lin.area_err << "  status " << lin.status
            << std::endl;
}

// Energy of channel mu under a quadratic calibration through the points,
// with its error from the fit covariance.
void QuadraticEnergy(const std::vector<Point> &points, Double_t slope_guess,
                     Double_t mu, Double_t &e, Double_t &e_err, Double_t &chi2,
                     Int_t &ndf) {
  Double_t pivot = 0;
  for (size_t i = 0; i < points.size(); i++)
    pivot += points[i].ch / points.size();
  TGraphErrors g((Int_t)points.size());
  for (size_t i = 0; i < points.size(); i++) {
    g.SetPoint(i, points[i].ch - pivot, points[i].e);
    g.SetPointError(
        i, 0,
        TMath::Sqrt(TMath::Power(slope_guess * points[i].ch_err, 2) +
                    TMath::Power(points[i].e_err, 2)));
  }
  TF1 f("quad_cal", "pol2", -1000, 1000);
  f.SetParameters(points[0].e, slope_guess, 0);
  TFitResultPtr r = g.Fit(&f, "QSN0");
  Double_t x = mu - pivot;
  e = f.Eval(x);
  Double_t grad[3] = {1.0, x, x * x};
  Double_t var = 0;
  for (Int_t a = 0; a < 3; a++)
    for (Int_t b = 0; b < 3; b++)
      var += grad[a] * grad[b] * r->CovMatrix(a, b);
  e_err = TMath::Sqrt(var);
  chi2 = f.GetChisquare();
  ndf = f.GetNDF();
}

void PrintCalibration(const TString &label, const Calibration &cal,
                      const std::vector<Point> &all) {
  std::cout << "SIGE CAL " << std::left << std::setw(30) << label << std::right
            << std::fixed << std::setprecision(6) << " E = " << cal.c0 << " + "
            << cal.c1 << " (ch - " << std::setprecision(2) << cal.pivot
            << ")  +- " << 1000 * cal.OffsetError() << " eV, slope +- "
            << std::setprecision(1) << 1e6 * TMath::Sqrt(cal.v11) / cal.c1
            << " ppm  chi2/ndf " << std::setprecision(2) << cal.chi2 << "/"
            << cal.ndf << std::endl;
  for (size_t i = 0; i < all.size(); i++) {
    Double_t resid = 1000 * (cal.Energy(all[i].ch) - all[i].e);
    Double_t err = 1000 * TMath::Abs(cal.c1) * all[i].ch_err;
    std::cout << "SIGE CAL   " << std::left << std::setw(30) << label << " "
              << std::setw(9) << all[i].name << std::right << std::fixed
              << std::setprecision(2) << " ch " << std::setw(8) << all[i].ch
              << " +- " << std::setw(5) << all[i].ch_err << "  E "
              << std::setprecision(4) << cal.Energy(all[i].ch) << "  resid "
              << std::setprecision(1) << std::setw(6) << resid << " +- " << err
              << " eV  sigma " << std::setprecision(0)
              << 1000 * cal.c1 * all[i].sigma << " eV" << std::endl;
  }
}

// ----------------------------------------------------- window scan helper

struct ScanResult {
  Double_t rms = 0, max_dev = 0;
  Int_t n = 0;
};

ScanResult ScanWindows(TH1D *h, const Calibration &cal, Double_t e_default) {
  ScanResult s;
  std::vector<Double_t> los, his;
  for (Double_t half = 1.0; half < 3.01; half += 0.25) {
    los.push_back(E_GE_GUESS - half);
    his.push_back(E_GE_GUESS + half);
  }
  for (Double_t shift = -0.5; shift < 0.51; shift += 0.25) {
    if (TMath::Abs(shift) < 1e-6)
      continue;
    los.push_back(WIN_LO + shift);
    his.push_back(WIN_HI + shift);
  }
  los.push_back(WIN_LO);
  his.push_back(72.0);
  los.push_back(65.0);
  his.push_back(WIN_HI);
  Double_t sum2 = 0;
  for (size_t i = 0; i < los.size(); i++) {
    GeFit g = FitGe(h, cal, los[i], his[i], 1);
    Double_t e = cal.Energy(g.mu);
    Double_t d = e - e_default;
    sum2 += d * d;
    s.max_dev = TMath::Max(s.max_dev, TMath::Abs(d));
    std::cout << "SIGE WIN sum " << std::fixed << std::setprecision(2) << los[i]
              << "-" << his[i] << " keV: " << std::setprecision(4) << e
              << " +- " << std::setprecision(1) << 1000 * cal.c1 * g.mu_err
              << " eV  (" << std::showpos << 1000 * d << std::noshowpos
              << " eV)  chi2/ndf " << g.chi2 << "/" << g.ndf << "  sigma "
              << std::setprecision(0) << 1000 * cal.c1 * g.sigma << " eV"
              << std::endl;
  }
  s.n = (Int_t)los.size();
  s.rms = TMath::Sqrt(sum2 / los.size());
  return s;
}

// ------------------------------------------------------------- plotting

// The channel histogram rebinned and turned into keV with a calibration.
TH1D *ToKeV(TH1D *h, const Calibration &cal, Int_t rebin, Double_t lo_kev,
            Double_t hi_kev, const TString &name) {
  Int_t first = h->FindBin(cal.Channel(lo_kev));
  Int_t n = (h->FindBin(cal.Channel(hi_kev)) - first) / rebin;
  Double_t e_lo = cal.Energy(h->GetBinLowEdge(first));
  Double_t e_hi = cal.Energy(h->GetBinLowEdge(first + n * rebin));
  TH1D *out = new TH1D(name, "", n, e_lo, e_hi);
  out->SetDirectory(nullptr);
  for (Int_t i = 0; i < n; i++) {
    Double_t sum = 0;
    for (Int_t k = 0; k < rebin; k++)
      sum += h->GetBinContent(first + i * rebin + k);
    out->SetBinContent(i + 1, sum);
    out->SetBinError(i + 1, TMath::Sqrt(sum));
  }
  return out;
}

Double_t Gauss(Double_t x, Double_t area, Double_t mu, Double_t sigma) {
  Double_t z = (x - mu) / sigma;
  return area * TMath::Exp(-0.5 * z * z) /
         (TMath::Sqrt(2.0 * TMath::Pi()) * sigma);
}

Double_t BackgroundPerChannel(const GeFit &g, Double_t ch) {
  Double_t d = ch - 0.5 * (g.lo + g.hi);
  return g.c0 + g.c1 * d + g.c2 * d * d;
}

Double_t ModelPerChannel(const GeFit &g, Double_t ch) {
  return BackgroundPerChannel(g, ch) + Gauss(ch, g.area, g.mu, g.sigma);
}

void DrawSumFit(TH1D *h, const Calibration &cal, const GeFit &lin,
                const GeFit &flat, const GeEnergy &result,
                const TString &label) {
  const Int_t rebin = 2;
  const Double_t show_lo = 65.5, show_hi = 72.0;
  TH1D *data = ToKeV(h, cal, rebin, show_lo, show_hi, "sum_kev");
  Double_t bin_ev = 1000.0 * cal.c1 * rebin;
  const Int_t n_points = 150;
  Double_t e_lo = cal.Energy(lin.lo), e_hi = cal.Energy(lin.hi);
  TGraph *g_lin = new TGraph(n_points), *g_flat = new TGraph(n_points),
         *g_background = new TGraph(n_points);
  for (Int_t i = 0; i < n_points; i++) {
    Double_t e = e_lo + (e_hi - e_lo) * i / (n_points - 1.0);
    Double_t ch = cal.Channel(e);
    g_lin->SetPoint(i, e, rebin * ModelPerChannel(lin, ch));
    g_flat->SetPoint(i, e, rebin * ModelPerChannel(flat, ch));
    g_background->SetPoint(i, e, rebin * BackgroundPerChannel(lin, ch));
  }
  TH1D *pulls = static_cast<TH1D *>(data->Clone("pulls"));
  pulls->Reset();
  for (Int_t b = 1; b <= data->GetNbinsX(); b++) {
    Double_t ch_lo = cal.Channel(data->GetBinLowEdge(b));
    if (ch_lo < lin.lo - 0.01 || ch_lo + rebin > lin.hi + 1.01)
      continue;
    Double_t expected = 0;
    for (Int_t k = 0; k < rebin; k++)
      expected += ModelPerChannel(lin, ch_lo + 0.5 + k);
    pulls->SetBinContent(b, (data->GetBinContent(b) - expected) /
                                TMath::Sqrt(expected));
  }
  std::vector<TString> labels;
  labels.push_back(label);
  labels.push_back("Gaussian + linear background");
  labels.push_back("Gaussian + flat background");
  labels.push_back("Linear background");
  Double_t legend_w = PlottingUtils::LegendWidthPx(labels);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)legend_w, 1100);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(plot, 300, upper, lower);
  upper->cd();
  PlottingUtils::ConfigureHistogram(
      data, kBlack, Form(";Energy [keV];Counts / %.0f eV", bin_ev));
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  Double_t data_max = data->GetMaximum();
  data->SetMinimum(0.9 * data->GetMinimum());
  data->SetMaximum(1.2 * data_max);
  data->GetXaxis()->SetLabelSize(0);
  data->GetXaxis()->SetTitleSize(0);
  data->Draw("E");
  PlottingUtils::ConfigureGraph(g_lin, kRed + 1);
  PlottingUtils::ConfigureGraph(g_flat, kAzure + 2);
  PlottingUtils::ConfigureGraph(g_background, kGray + 2);
  g_flat->SetLineStyle(2);
  g_background->SetLineStyle(3);
  g_background->Draw("L SAME");
  g_flat->Draw("L SAME");
  g_lin->Draw("L SAME");
  Double_t w_lines[2] = {W_KB1, W_KB2};
  const char *w_names[2] = {"W K#beta_{1}", "W K#beta_{2}"};
  for (Int_t k = 0; k < 2; k++) {
    TLine *l =
        new TLine(w_lines[k], data->GetMinimum(), w_lines[k], 1.04 * data_max);
    l->SetLineStyle(3);
    l->SetLineColor(kGreen + 3);
    l->Draw();
    TLatex *t = new TLatex(w_lines[k] + 0.05, 1.02 * data_max, w_names[k]);
    t->SetTextFont(43);
    t->SetTextSize(24);
    t->SetTextColor(kGreen + 3);
    t->Draw();
  }
  TLatex *res = new TLatex(
      show_lo + 0.15, 1.12 * data_max,
      Form("E_{#gamma} = %.4f keV (fit #pm %.0f eV), #chi^{2}/ndf = %.0f/%d",
           result.e, 1000 * result.fit, lin.chi2, lin.ndf));
  res->SetTextFont(43);
  res->SetTextSize(28);
  res->Draw();
  lower->cd();
  PlottingUtils::ConfigureHistogram(pulls, kBlack, ";Energy [keV];Pull");
  pulls->SetMinimum(-3.9);
  pulls->SetMaximum(3.9);
  pulls->GetYaxis()->SetNdivisions(504);
  pulls->SetFillColor(kGray + 1);
  pulls->SetFillStyle(1001);
  pulls->Draw("HIST");
  TLine *zero = new TLine(show_lo, 0, data->GetXaxis()->GetXmax(), 0);
  zero->Draw();
  legend_pad->cd();
  TLegend *leg = new TLegend(0.05, 0.92, 0.98, 0.62);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetTextSize(30);
  leg->SetTextFont(43);
  leg->AddEntry(data, labels[0], "lep");
  leg->AddEntry(g_lin, labels[1], "l");
  leg->AddEntry(g_flat, labels[2], "l");
  leg->AddEntry(g_background, labels[3], "l");
  leg->Draw();
  PlottingUtils::ScaleFigure(plot);
  PlottingUtils::SaveFigure(canvas, "psi_si_ge_sum_fit", "",
                            PlotSaveOptions::kLINEAR);
}

void DrawCalibration(const Calibration &cal_a, const Calibration &cal_b,
                     const std::vector<Point> &lines) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TGraphErrors *ga = new TGraphErrors((Int_t)lines.size());
  TGraphErrors *gb = new TGraphErrors((Int_t)lines.size());
  for (size_t i = 0; i < lines.size(); i++) {
    ga->SetPoint(i, lines[i].e - 0.25,
                 1000 * (cal_a.Energy(lines[i].ch) - lines[i].e));
    ga->SetPointError(i, 0, 1000 * cal_a.c1 * lines[i].ch_err);
    gb->SetPoint(i, lines[i].e + 0.25,
                 1000 * (cal_b.Energy(lines[i].ch) - lines[i].e));
    gb->SetPointError(i, 0, 1000 * cal_b.c1 * lines[i].ch_err);
  }
  PlottingUtils::ConfigureGraph(
      ga, kBlack, ";Energy [keV];Fitted #minus reference energy [eV]");
  PlottingUtils::ConfigureGraph(gb, kRed + 1);
  gb->SetMarkerStyle(21);
  Double_t y_lo = -60, y_hi = 80;
  for (Int_t i = 0; i < ga->GetN(); i++) {
    y_lo = TMath::Min(y_lo, TMath::Min(ga->GetY()[i] - ga->GetEY()[i],
                                       gb->GetY()[i] - gb->GetEY()[i]) -
                                20);
    y_hi = TMath::Max(y_hi, TMath::Max(ga->GetY()[i] + ga->GetEY()[i],
                                       gb->GetY()[i] + gb->GetEY()[i]) +
                                60);
  }
  ga->GetXaxis()->SetLimits(50, 84);
  ga->SetMinimum(y_lo);
  ga->SetMaximum(y_hi);
  ga->Draw("AP");
  gb->Draw("P SAME");
  TLine *zero = new TLine(50, 0, 84, 0);
  zero->SetLineStyle(2);
  zero->Draw();
  TLine *ge = new TLine(ReferenceLines::GE73M_ENSDF, y_lo,
                        ReferenceLines::GE73M_ENSDF, 0.5 * y_hi);
  ge->SetLineColor(kGreen + 3);
  ge->SetLineStyle(7);
  ge->Draw();
  TLatex *t = new TLatex(ReferenceLines::GE73M_ENSDF + 0.4, y_lo + 10, "Ge");
  t->SetTextFont(43);
  t->SetTextSize(28);
  t->SetTextColor(kGreen + 3);
  t->Draw();
  const char *names[6] = {"Ba 53.2",        "Am 59.5", "Pb K#alpha_{2}",
                          "Pb K#alpha_{1}", "Ba 79.6", "Ba 81.0"};
  for (size_t i = 0; i < lines.size() && i < 6; i++) {
    TLatex *l =
        new TLatex(lines[i].e, y_lo + (i == 3 || i == 5 ? 8 : 20), names[i]);
    l->SetTextFont(43);
    l->SetTextSize(22);
    l->SetTextAlign(21);
    l->Draw();
  }
  TLegend *leg = PlottingUtils::AddLegend(0.16, 0.73, 0.73, 0.89);
  leg->AddEntry(ga, "Am + Pb K#alpha (beam runs) calibration", "lep");
  leg->AddEntry(gb, "Am + ^{133}Ba calibration", "lep");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_si_calibration", "",
                            PlotSaveOptions::kLINEAR);
}

struct RunValue {
  TString label;
  Double_t e_free, err_free, e_fixed, err_fixed;
};

// Per-run energies (Am + the run's own Pb Ka calibration) with the width
// free and fixed, against the summed-spectrum result and the per-run
// weighted means.
void DrawPerRun(const std::vector<RunValue> &values, Double_t sigma_fixed_ev,
                Double_t sum_e, Double_t sum_err, Double_t mean_free,
                Double_t mean_fixed) {
  Int_t n = (Int_t)values.size();
  std::vector<TString> labels;
  labels.push_back("Per run, width free");
  labels.push_back(Form("Per run, width fixed (%.0f eV)", sigma_fixed_ev));
  labels.push_back("Weighted mean, width free");
  labels.push_back("Weighted mean, width fixed");
  labels.push_back("Summed runs, total error");
  Double_t legend_w = PlottingUtils::LegendWidthPx(labels);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)legend_w);
  plot->cd();
  plot->SetBottomMargin(0.27);
  const Double_t y_lo = 68.45, y_hi = 69.12;
  TH1D *frame = new TH1D("per_run_frame", "", n, 0.5, n + 0.5);
  frame->SetDirectory(nullptr);
  for (Int_t i = 0; i < n; i++)
    frame->GetXaxis()->SetBinLabel(i + 1, values[i].label);
  PlottingUtils::ConfigureHistogram(frame, kBlack, ";;E_{#gamma} [keV]");
  frame->GetXaxis()->LabelsOption("v");
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->SetMinimum(y_lo);
  frame->SetMaximum(y_hi);
  frame->SetStats(0);
  frame->Draw("AXIS");
  TLine *divider = new TLine(n - 0.5, y_lo, n - 0.5, y_hi);
  divider->SetLineStyle(3);
  divider->SetLineColor(kGray + 2);
  divider->Draw();
  TBox *band = new TBox(0.5, sum_e - sum_err, n + 0.5, sum_e + sum_err);
  band->SetFillColorAlpha(kGray, 0.6);
  band->SetLineWidth(0);
  band->Draw();
  TLine *sum_line = new TLine(0.5, sum_e, n + 0.5, sum_e);
  sum_line->SetLineColor(kGray + 2);
  sum_line->SetLineWidth(2);
  sum_line->Draw();
  TLine *free_line = new TLine(0.5, mean_free, n + 0.5, mean_free);
  free_line->SetLineColor(kBlack);
  free_line->SetLineStyle(2);
  free_line->SetLineWidth(2);
  free_line->Draw();
  TLine *fixed_line = new TLine(0.5, mean_fixed, n + 0.5, mean_fixed);
  fixed_line->SetLineColor(kAzure + 2);
  fixed_line->SetLineStyle(7);
  fixed_line->SetLineWidth(2);
  fixed_line->Draw();
  TGraphErrors *g_free = new TGraphErrors(n), *g_fixed = new TGraphErrors(n);
  for (Int_t i = 0; i < n; i++) {
    g_free->SetPoint(i, i + 1 - 0.12, values[i].e_free);
    g_free->SetPointError(i, 0, values[i].err_free);
    g_fixed->SetPoint(i, i + 1 + 0.12, values[i].e_fixed);
    g_fixed->SetPointError(i, 0, values[i].err_fixed);
  }
  PlottingUtils::ConfigureGraph(g_free, kBlack);
  PlottingUtils::ConfigureGraph(g_fixed, kAzure + 2);
  g_fixed->SetMarkerStyle(24);
  g_free->Draw("P SAME");
  g_fixed->Draw("P SAME");
  frame->Draw("AXIS SAME");
  legend_pad->cd();
  TLegend *leg = new TLegend(0.05, 0.92, 0.98, 0.55);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetTextSize(30);
  leg->SetTextFont(43);
  leg->AddEntry(g_free, labels[0], "lep");
  leg->AddEntry(g_fixed, labels[1], "lep");
  leg->AddEntry(free_line, labels[2], "l");
  leg->AddEntry(fixed_line, labels[3], "l");
  leg->AddEntry(band, labels[4], "f");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_si_ge_per_run", "",
                            PlotSaveOptions::kLINEAR);
}

// ------------------------------------------------------------ combination

struct Combined {
  Double_t e = 0, fit = 0, chi2 = 0;
  Int_t ndf = 0;
};

Combined WeightedMean(const std::vector<Double_t> &e,
                      const std::vector<Double_t> &err) {
  Combined c;
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < e.size(); i++) {
    Double_t w = 1.0 / (err[i] * err[i]);
    sw += w;
    swx += w * e[i];
  }
  c.e = swx / sw;
  c.fit = 1.0 / TMath::Sqrt(sw);
  for (size_t i = 0; i < e.size(); i++)
    c.chi2 += TMath::Power((e[i] - c.e) / err[i], 2);
  c.ndf = (Int_t)e.size() - 1;
  return c;
}

// The spread of a set of line positions (one per run) against their
// weighted mean: chi2/ndf, and the extra scatter that brings it to 1.
void PrintStability(const TString &label, const std::vector<Double_t> &ch,
                    const std::vector<Double_t> &ch_err, Double_t ev_per_ch) {
  if (ch.size() < 2)
    return;
  Combined c = WeightedMean(ch, ch_err);
  Double_t extra = 0;
  if (c.chi2 > c.ndf) {
    Double_t lo = 0, hi = 20;
    for (Int_t it = 0; it < 60; it++) {
      Double_t mid = 0.5 * (lo + hi);
      std::vector<Double_t> inflated(ch_err.size());
      for (size_t i = 0; i < ch_err.size(); i++)
        inflated[i] = TMath::Sqrt(ch_err[i] * ch_err[i] + mid * mid);
      if (WeightedMean(ch, inflated).chi2 > c.ndf)
        lo = mid;
      else
        hi = mid;
    }
    extra = 0.5 * (lo + hi);
  }
  std::cout << "SIGE STABILITY " << std::left << std::setw(8) << label
            << std::right << " n " << ch.size() << " mean " << std::fixed
            << std::setprecision(3) << c.e << " +- " << c.fit
            << " ch  chi2/ndf " << std::setprecision(1) << c.chi2 << "/"
            << c.ndf << "  p " << std::setprecision(3)
            << TMath::Prob(c.chi2, c.ndf) << "  extra scatter "
            << std::setprecision(2) << extra << " ch (" << std::setprecision(1)
            << extra * ev_per_ch << " eV)" << std::endl;
}

TString SpectraFile(Estimator estimator) {
  return TString(SPECTRA_DIR) + "/ge_spectra_" + EstimatorName(estimator) +
         ".root";
}

void WriteResult(const TString &path, Estimator estimator,
                 const GeEnergy &result, const GeEnergy &cal_b,
                 const GeEnergy &cal_c, const Combined &days_b) {
  std::ofstream out(path.Data());
  out << "# Si (PSI) Ge energy, SiGeBudget.cpp, " << EstimatorName(estimator)
      << " energy; summed main runs, Am + Pb Ka calibration; keV" << std::endl;
  out << std::fixed << std::setprecision(6);
  out << "ge " << result.e << std::endl;
  out << "total " << result.total << std::endl;
  out << "fit " << result.fit << std::endl;
  out << "calibration " << result.cal << std::endl;
  out << "background " << result.background << std::endl;
  out << "window " << result.window << std::endl;
  out << "# alternatives" << std::endl;
  out << "ge_cal_b " << cal_b.e << std::endl;
  out << "ge_cal_c " << cal_c.e << std::endl;
  out << "ge_days_cal_b " << days_b.e << std::endl;
  out << "ge_days_cal_b_err " << days_b.fit << std::endl;
}

// ----------------------------------------------------------------- budget

void Budget(Estimator estimator, const TString &project_root) {
  // ------------------------------------------------------------ spectra
  std::vector<GeRun> runs = GeRuns();
  std::vector<TH1D *> h_run;
  TH1D *h_sum = NewSpectrum("ge_sum");
  TH1D *h_sum_all = NewSpectrum("ge_sum_all");
  std::map<Int_t, TH1D *> h_day;
  Double_t seconds_main = 0;
  for (size_t i = 0; i < runs.size(); i++) {
    TH1D *h = NewSpectrum("ge_" + runs[i].label);
    FillSpectrum(h, runs[i].run, estimator, Selection());
    h_run.push_back(h);
    Double_t seconds = 0;
    Long64_t entries = 0;
    RunSpan(runs[i].run, seconds, entries);
    std::cout << "SIGE RUN " << std::left << std::setw(18) << runs[i].label
              << std::right << " day " << runs[i].day << std::fixed
              << std::setprecision(0) << "  span " << std::setw(6) << seconds
              << " s  records " << entries << "  rate " << std::setprecision(1)
              << entries / seconds << " Hz  48-90 keV " << std::setprecision(0)
              << h->Integral(1150, 2150) << std::endl;
    if (runs[i].main) {
      h_sum->Add(h);
      seconds_main += seconds;
    }
    h_sum_all->Add(h);
    if (!h_day.count(runs[i].day))
      h_day[runs[i].day] = NewSpectrum(Form("ge_day%d", runs[i].day));
    h_day[runs[i].day]->Add(h);
  }
  std::cout << "SIGE RUN main total span " << std::setprecision(0)
            << seconds_main << " s" << std::endl;
  TH1D *h_am = SumSpectrum(AmRuns(0), "am_all", estimator);
  TH1D *h_ba = SumSpectrum(BaRuns(0), "ba_all", estimator);
  std::map<Int_t, TH1D *> h_am_day, h_ba_day;
  for (Int_t d = 10; d <= 13; d++) {
    h_am_day[d] = SumSpectrum(AmRuns(d), Form("am_d%d", d), estimator);
    h_ba_day[d] = SumSpectrum(BaRuns(d), Form("ba_d%d", d), estimator);
  }

  // ------------------------------------------------------ calibrations
  LineSetting nominal;
  Point am = FitAm(h_am, nominal);
  Point ba53 = FitBa53(h_ba, nominal);
  Point ba80, ba81, ka2, ka1, source_ka2, source_ka1;
  FitBa8081(h_ba, nominal, ba80, ba81);
  FitPb(h_sum, nominal, ka2, ka1);
  FitPb(h_ba, nominal, source_ka2, source_ka1);
  std::vector<Point> points_a, points_b, points_c, all_lines;
  points_a.push_back(am);
  points_a.push_back(ka2);
  points_a.push_back(ka1);
  points_b.push_back(am);
  points_b.push_back(ba53);
  points_b.push_back(ba80);
  points_b.push_back(ba81);
  points_c = points_b;
  points_c.push_back(ka2);
  points_c.push_back(ka1);
  all_lines.push_back(ba53);
  all_lines.push_back(am);
  all_lines.push_back(ka2);
  all_lines.push_back(ka1);
  all_lines.push_back(ba80);
  all_lines.push_back(ba81);
  Calibration cal_a = Calibrate(points_a);
  Calibration cal_b = Calibrate(points_b);
  Calibration cal_c = Calibrate(points_c);
  PrintCalibration("A: Am + Pb Ka (beam)", cal_a, all_lines);
  PrintCalibration("B: Am + Ba (source)", cal_b, all_lines);
  PrintCalibration("C: all six lines", cal_c, all_lines);
  std::cout << "SIGE PB beam Ka2 " << std::setprecision(2) << ka2.ch << " +- "
            << ka2.ch_err << "  Ka1 " << ka1.ch << " +- " << ka1.ch_err
            << " | Ba source runs Ka2 " << source_ka2.ch << " +- "
            << source_ka2.ch_err << "  Ka1 " << source_ka1.ch << " +- "
            << source_ka1.ch_err << "  (beam - source: Ka2 "
            << ka2.ch - source_ka2.ch << ", Ka1 " << ka1.ch - source_ka1.ch
            << " ch)" << std::endl;

  // Width at 68.75 keV from the source lines, sigma^2 = a + b E, as an
  // independent value for the fixed-width variants.
  Double_t s_am = cal_b.c1 * am.sigma, s_ba = cal_b.c1 * ba81.sigma;
  Double_t slope_s2 = (s_ba * s_ba - s_am * s_am) /
                      (ReferenceLines::BA133_81 - ReferenceLines::AM241);
  Double_t sigma_ge_kev = TMath::Sqrt(
      s_am * s_am + slope_s2 * (E_GE_GUESS - ReferenceLines::AM241));
  Double_t sigma_ge_ch = sigma_ge_kev / cal_b.c1;
  std::cout << "SIGE WIDTH Am " << std::setprecision(0) << 1000 * s_am
            << " eV, Ba81 " << 1000 * s_ba << " eV, Pb Ka (beam) "
            << 1000 * cal_a.c1 * ka1.sigma
            << " eV -> interpolated at 68.75: " << 1000 * sigma_ge_kev
            << " eV (" << std::setprecision(2) << sigma_ge_ch << " ch)"
            << std::endl;

  // ------------------------------------------------------- summed Ge fit
  GeFit lin = FitGe(h_sum, cal_a, WIN_LO, WIN_HI, 1);
  GeFit flat = FitGe(h_sum, cal_a, WIN_LO, WIN_HI, 0);
  GeFit quad = FitGe(h_sum, cal_a, WIN_LO, WIN_HI, 2);
  GeFit lin_fixed = FitGe(h_sum, cal_a, WIN_LO, WIN_HI, 1, sigma_ge_ch);
  GeFit flat_fixed = FitGe(h_sum, cal_a, WIN_LO, WIN_HI, 0, sigma_ge_ch);
  GeFit with_w = FitGeWithW(h_sum, cal_a, WIN_LO, WIN_HI);
  ScanResult scan = ScanWindows(h_sum, cal_a, cal_a.Energy(lin.mu));
  std::cout << "SIGE SUM Ge mu " << std::setprecision(3) << lin.mu << " +- "
            << lin.mu_err << " ch  sigma " << lin.sigma << " +- "
            << lin.sigma_err << " ch  area " << std::setprecision(0) << lin.area
            << " +- " << lin.area_err << "  window " << std::setprecision(2)
            << WIN_LO << "-" << WIN_HI << " keV = ch " << lin.lo << "-"
            << lin.hi << std::endl;
  std::cout << "SIGE SUM window scan: RMS " << std::setprecision(1)
            << 1000 * scan.rms << " eV, max |dev| " << 1000 * scan.max_dev
            << " eV over " << scan.n << " windows" << std::endl;
  GeEnergy sum_a = Energy(lin, flat, cal_a, scan.rms);
  GeEnergy sum_b = Energy(lin, flat, cal_b, scan.rms * cal_b.c1 / cal_a.c1);
  GeEnergy sum_c = Energy(lin, flat, cal_c, scan.rms * cal_c.c1 / cal_a.c1);
  PrintEnergy("SUM main, cal A (Am + Pb Ka)", sum_a, lin, cal_a);
  PrintEnergy("SUM main, cal B (Am + Ba)", sum_b, lin, cal_b);
  PrintEnergy("SUM main, cal C (all six)", sum_c, lin, cal_c);
  std::cout << "SIGE SUM variants (cal A): linear " << std::setprecision(4)
            << cal_a.Energy(lin.mu) << "  flat " << cal_a.Energy(flat.mu)
            << "  quadratic " << cal_a.Energy(quad.mu) << "  | sigma fixed "
            << std::setprecision(0) << 1000 * sigma_ge_kev << " eV: linear "
            << std::setprecision(4) << cal_a.Energy(lin_fixed.mu) << " +- "
            << std::setprecision(1) << 1000 * cal_a.c1 * lin_fixed.mu_err
            << " eV, flat " << std::setprecision(4)
            << cal_a.Energy(flat_fixed.mu)
            << "  | with W Kb group: " << cal_a.Energy(with_w.mu) << " +- "
            << std::setprecision(1) << 1000 * cal_a.c1 * with_w.mu_err
            << " eV, W Kb1 area " << std::setprecision(0) << with_w.w_area
            << " +- " << with_w.w_area_err << " (Ge area " << with_w.area << ")"
            << std::endl;
  std::cout << "SIGE SUM chi2/ndf linear " << std::setprecision(1) << lin.chi2
            << "/" << lin.ndf << "  flat " << flat.chi2 << "/" << flat.ndf
            << "  quadratic " << quad.chi2 << "/" << quad.ndf << "  bkg slope "
            << std::setprecision(4) << lin.c1 / lin.c0 * 100 << " %/ch"
            << std::endl;

  // W K alpha in the Ge runs: positions fixed by cal A, width by the Am
  // line.
  {
    std::vector<LineSpec> w(2);
    w[0].mean_guess = cal_a.Channel(W_KA2);
    w[0].fix_mean = kTRUE;
    w[0].allow_negative = kTRUE;
    w[1].mean_guess = cal_a.Channel(W_KA1);
    w[1].fix_mean = kTRUE;
    w[1].allow_negative = kTRUE;
    FitOutput fw = FitLines(h_sum, cal_a.Channel(56.5), cal_a.Channel(62.0), w,
                            1, 1, am.sigma, kTRUE);
    std::cout << "SIGE W Ka2 (57.98) area " << std::setprecision(0)
              << fw.area[0] << " +- " << fw.area_err[0]
              << "  W Ka1 (59.32) area " << fw.area[1] << " +- "
              << fw.area_err[1] << "  chi2/ndf " << std::setprecision(1)
              << fw.chi2 << "/" << fw.ndf << std::endl;
  }
  {
    std::vector<LineSpec> one(1);
    one[0].mean_guess = cal_a.Channel(GE73M_53);
    one[0].mean_range = 8;
    FitOutput f53 = FitLines(h_sum, cal_a.Channel(51.5), cal_a.Channel(55.5),
                             one, 1, 1, am.sigma);
    std::cout << "SIGE 53keV line: ch " << std::setprecision(2) << f53.mean[0]
              << " +- " << f53.mean_err[0] << "  area " << std::setprecision(0)
              << f53.area[0] << " +- " << f53.area_err[0] << "  E(A) "
              << std::setprecision(4) << cal_a.Energy(f53.mean[0]) << "  E(B) "
              << cal_b.Energy(f53.mean[0]) << " +- " << std::setprecision(1)
              << 1000 * cal_b.c1 * f53.mean_err[0] << " eV  sigma "
              << std::setprecision(0) << 1000 * cal_b.c1 * f53.sigma[0] << " eV"
              << std::endl;
  }

  // Calibration variants: the source calibration without the weak Ba 79.61
  // line, and quadratic calibrations through the source lines and all six.
  {
    std::vector<Point> p;
    p.push_back(am);
    p.push_back(ba53);
    p.push_back(ba81);
    Calibration cal = Calibrate(p);
    std::cout << "SIGE CALVAR B without Ba 79.61: Ge " << std::setprecision(4)
              << cal.Energy(lin.mu) << " (" << std::showpos
              << std::setprecision(1) << 1000 * (cal.Energy(lin.mu) - sum_b.e)
              << std::noshowpos << " eV vs B)  chi2/ndf " << cal.chi2 << "/"
              << cal.ndf << "  Pb Ka2/Ka1 resid " << std::showpos
              << 1000 * (cal.Energy(ka2.ch) - ReferenceLines::PB_KA2) << " / "
              << 1000 * (cal.Energy(ka1.ch) - ReferenceLines::PB_KA1)
              << std::noshowpos << " eV" << std::endl;
    Double_t e = 0, e_err = 0, chi2 = 0;
    Int_t ndf = 0;
    QuadraticEnergy(points_b, cal_b.c1, lin.mu, e, e_err, chi2, ndf);
    std::cout << "SIGE CALVAR B quadratic: Ge " << std::setprecision(4) << e
              << " +- " << std::setprecision(1) << 1000 * e_err
              << " eV (cal)  (" << std::showpos << 1000 * (e - sum_b.e)
              << std::noshowpos << " eV vs B)  chi2/ndf " << chi2 << "/" << ndf
              << std::endl;
    QuadraticEnergy(points_c, cal_c.c1, lin.mu, e, e_err, chi2, ndf);
    std::cout << "SIGE CALVAR C quadratic: Ge " << std::setprecision(4) << e
              << " +- " << std::setprecision(1) << 1000 * e_err
              << " eV (cal)  (" << std::showpos << 1000 * (e - sum_c.e)
              << std::noshowpos << " eV vs C)  chi2/ndf " << chi2 << "/" << ndf
              << std::endl;
  }

  // Isolated events: no other processed record within 100 us (the
  // digitiser records none within 10 us).
  {
    Selection isolated;
    isolated.min_gap_us = 100;
    std::vector<TString> main_runs;
    for (size_t i = 0; i < runs.size(); i++)
      if (runs[i].main)
        main_runs.push_back(runs[i].run);
    TH1D *h_isolated =
        SumSpectrum(main_runs, "ge_sum_isolated", estimator, isolated);
    Point k2, k1;
    FitPb(h_isolated, nominal, k2, k1);
    std::vector<Point> p;
    p.push_back(am);
    p.push_back(k2);
    p.push_back(k1);
    Calibration cal = Calibrate(p);
    GeFit l = FitGe(h_isolated, cal, WIN_LO, WIN_HI, 1);
    std::cout << "SIGE ISOLATED (gap > 100 us): kept " << std::setprecision(4)
              << h_isolated->Integral() / h_sum->Integral()
              << " of events; Pb Ka1 " << std::setprecision(2) << k1.ch
              << " (was " << ka1.ch << "), Ge mu " << l.mu << " (was " << lin.mu
              << "); Ge(A) " << std::setprecision(4) << cal.Energy(l.mu) << " ("
              << std::showpos << std::setprecision(1)
              << 1000 * (cal.Energy(l.mu) - sum_a.e) << std::noshowpos
              << " eV), sigma " << std::setprecision(0)
              << 1000 * cal.c1 * l.sigma << " eV, Pb sigma "
              << 1000 * cal.c1 * k1.sigma << " eV" << std::endl;
    delete h_isolated;
  }

  // The two rise bands alone, each on its own Am + Pb Ka calibration: the
  // trapezoid's gain to the digitiser's differs by 0.1% between them.
  {
    const Int_t n_bands = 2;
    const Double_t band_lo[n_bands] = {150, 280}, band_hi[n_bands] = {230, 360};
    std::vector<TString> main_runs;
    for (size_t i = 0; i < runs.size(); i++)
      if (runs[i].main)
        main_runs.push_back(runs[i].run);
    for (Int_t b = 0; b < n_bands; b++) {
      Selection band;
      band.rise_min_ns = band_lo[b];
      band.rise_max_ns = band_hi[b];
      TH1D *h_band = SumSpectrum(main_runs, "ge_sum_band", estimator, band);
      TH1D *h_am_band = SumSpectrum(AmRuns(0), "am_band", estimator, band);
      Point am_band = FitAm(h_am_band, nominal), k2, k1;
      FitPb(h_band, nominal, k2, k1);
      std::vector<Point> p;
      p.push_back(am_band);
      p.push_back(k2);
      p.push_back(k1);
      Calibration cal = Calibrate(p);
      GeFit l = FitGe(h_band, cal, WIN_LO, WIN_HI, 1);
      GeFit fl = FitGe(h_band, cal, WIN_LO, WIN_HI, 0);
      PrintEnergy(Form("RISE %.0f-%.0f ns (%.3f of Ge runs), cal A", band_lo[b],
                       band_hi[b], h_band->Integral() / h_sum->Integral()),
                  Energy(l, fl, cal, 0), l, cal);
      delete h_band;
      delete h_am_band;
    }
  }

  // The sum including BeamTurnOn_4.
  {
    Point k2, k1;
    FitPb(h_sum_all, nominal, k2, k1);
    std::vector<Point> p;
    p.push_back(am);
    p.push_back(k2);
    p.push_back(k1);
    Calibration cal = Calibrate(p);
    GeFit l = FitGe(h_sum_all, cal, WIN_LO, WIN_HI, 1);
    GeFit fl = FitGe(h_sum_all, cal, WIN_LO, WIN_HI, 0);
    PrintEnergy("SUM incl. BeamTurnOn_4, cal A", Energy(l, fl, cal, 0), l, cal);
    PrintEnergy("SUM incl. BeamTurnOn_4, cal B", Energy(l, fl, cal_b, 0), l,
                cal_b);
  }

  // Calibration-line fit variations: every calibration line refitted with a
  // narrower or wider window or a flat background.
  {
    LineSetting settings[4];
    settings[0].widen = -10;
    settings[1].widen = 10;
    settings[2].background = 0;
    settings[3].widen = 20;
    const char *names[4] = {"windows -10 ch", "windows +10 ch", "flat bkg",
                            "windows +20 ch"};
    for (Int_t k = 0; k < 4; k++) {
      Point v_am = FitAm(h_am, settings[k]);
      Point v_53 = FitBa53(h_ba, settings[k]);
      Point v_80, v_81, v_k2, v_k1;
      FitBa8081(h_ba, settings[k], v_80, v_81);
      FitPb(h_sum, settings[k], v_k2, v_k1);
      std::vector<Point> pa, pb;
      pa.push_back(v_am);
      pa.push_back(v_k2);
      pa.push_back(v_k1);
      pb.push_back(v_am);
      pb.push_back(v_53);
      pb.push_back(v_80);
      pb.push_back(v_81);
      Calibration ca = Calibrate(pa), cb = Calibrate(pb);
      std::cout << "SIGE CALVAR " << std::left << std::setw(16) << names[k]
                << std::right << std::setprecision(2) << " Am " << v_am.ch
                << " Ka2 " << v_k2.ch << " Ka1 " << v_k1.ch << " Ba53 "
                << v_53.ch << " Ba80 " << v_80.ch << " Ba81 " << v_81.ch
                << " -> Ge(A) " << std::showpos << std::setprecision(1)
                << 1000 * (ca.Energy(lin.mu) - sum_a.e) << " eV, Ge(B) "
                << 1000 * (cb.Energy(lin.mu) - sum_b.e) << std::noshowpos
                << " eV" << std::endl;
    }
  }

  // ----------------------------------------------------------- per run
  std::vector<RunValue> plot_values;
  std::vector<Double_t> e_a, f_a, e_a_flat, f_a_flat, e_b, f_b, e_b_flat,
      f_b_flat, e_a_fix, f_a_fix;
  std::vector<Double_t> e_a8, f_a8, e_b8, f_b8;
  std::vector<Double_t> d_am, d_k2, d_k1, s_k2, s_k1;
  std::vector<Double_t> cal_b_err;
  std::vector<Double_t> pb_ka1_ch, pb_ka1_err;
  for (size_t i = 0; i < runs.size(); i++) {
    TH1D *h = h_run[i];
    Point k2, k1;
    FitPb(h, nominal, k2, k1);
    pb_ka1_ch.push_back(k1.ch);
    pb_ka1_err.push_back(k1.ch_err);
    std::vector<Point> p;
    p.push_back(am);
    p.push_back(k2);
    p.push_back(k1);
    Calibration cal = Calibrate(p);
    GeFit l = FitGe(h, cal, WIN_LO, WIN_HI, 1);
    GeFit fl = FitGe(h, cal, WIN_LO, WIN_HI, 0);
    GeFit lfix = FitGe(h, cal, WIN_LO, WIN_HI, 1, sigma_ge_ch);
    GeEnergy ra = Energy(l, fl, cal, 0);
    GeEnergy rb = Energy(l, fl, cal_b, 0);
    TString tag = runs[i].label + (l.ok ? "" : " [weak fit]");
    PrintEnergy("RUN " + tag + " cal A", ra, l, cal);
    PrintEnergy("RUN " + tag + " cal B", rb, l, cal_b);
    GeFit flfix = FitGe(h, cal, WIN_LO, WIN_HI, 0, sigma_ge_ch);
    std::cout << "SIGE RUN " << std::left << std::setw(18) << runs[i].label
              << std::right << " Pb Ka2 " << std::setprecision(2) << k2.ch
              << " +- " << k2.ch_err << " Ka1 " << k1.ch << " +- " << k1.ch_err
              << " (cal B resid " << std::showpos << std::setprecision(1)
              << 1000 * (cal_b.Energy(k2.ch) - ReferenceLines::PB_KA2) << " / "
              << 1000 * (cal_b.Energy(k1.ch) - ReferenceLines::PB_KA1)
              << std::noshowpos << " eV)  sigma-fixed Ge(A) "
              << std::setprecision(4) << cal.Energy(lfix.mu) << " +- "
              << std::setprecision(1) << 1000 * cal.c1 * lfix.mu_err
              << " eV (flat " << std::showpos
              << 1000 * (cal.Energy(flfix.mu) - cal.Energy(lfix.mu))
              << std::noshowpos << " eV)  Ge(B) " << std::setprecision(4)
              << cal_b.Energy(lfix.mu) << "  area " << std::setprecision(0)
              << lfix.area << " +- " << lfix.area_err << " ("
              << std::setprecision(1) << lfix.area / lfix.area_err
              << " sigma)  chi2/ndf " << lfix.chi2 << "/" << lfix.ndf
              << "  flat-bkg Ge(A) " << std::setprecision(4)
              << cal.Energy(fl.mu) << std::endl;
    // Derivatives of this run's cal-A energy with respect to the line
    // positions.
    const Double_t delta = 0.1;
    std::vector<Point> q = p;
    q[0].ch += delta;
    Double_t dam = (Calibrate(q).Energy(l.mu) - ra.e) / delta;
    q = p;
    q[1].ch += delta;
    Double_t dk2 = (Calibrate(q).Energy(l.mu) - ra.e) / delta;
    q = p;
    q[2].ch += delta;
    Double_t dk1 = (Calibrate(q).Energy(l.mu) - ra.e) / delta;
    Double_t fix_err = TMath::Sqrt(TMath::Power(cal.c1 * lfix.mu_err, 2) +
                                   TMath::Power(cal.Error(lfix.mu), 2));
    TString short_label = runs[i].label;
    short_label.ReplaceAll("GeTarget", "GeTgt")
        .ReplaceAll("SiYesGe_1056", "SiYesGe")
        .ReplaceAll("LongWf500V_1001", "LWf500V");
    short_label.ReplaceAll("LongWf_1517", "LWf")
        .ReplaceAll("BeamTurnOn", "TurnOn")
        .ReplaceAll("Ge_1715", "Ge");
    RunValue rv = {Form("10/%02d %s", runs[i].day, short_label.Data()), ra.e,
                   TMath::Sqrt(ra.fit * ra.fit + ra.cal * ra.cal),
                   cal.Energy(lfix.mu), fix_err};
    plot_values.push_back(rv);
    e_a8.push_back(ra.e);
    f_a8.push_back(ra.fit);
    e_b8.push_back(rb.e);
    f_b8.push_back(rb.fit);
    if (!runs[i].main)
      continue;
    e_a.push_back(ra.e);
    f_a.push_back(ra.fit);
    e_a_flat.push_back(cal.Energy(fl.mu));
    f_a_flat.push_back(cal.c1 * fl.mu_err);
    e_b.push_back(rb.e);
    f_b.push_back(rb.fit);
    e_b_flat.push_back(cal_b.Energy(fl.mu));
    f_b_flat.push_back(cal_b.c1 * fl.mu_err);
    e_a_fix.push_back(cal.Energy(lfix.mu));
    f_a_fix.push_back(cal.c1 * lfix.mu_err);
    d_am.push_back(dam);
    d_k2.push_back(dk2);
    d_k1.push_back(dk1);
    s_k2.push_back(k2.ch_err);
    s_k1.push_back(k1.ch_err);
    cal_b_err.push_back(rb.cal);
  }
  Combined ca = WeightedMean(e_a, f_a),
           ca_flat = WeightedMean(e_a_flat, f_a_flat);
  Combined cb = WeightedMean(e_b, f_b),
           cb_flat = WeightedMean(e_b_flat, f_b_flat);
  Combined ca_fix = WeightedMean(e_a_fix, f_a_fix);
  Combined ca8 = WeightedMean(e_a8, f_a8), cb8 = WeightedMean(e_b8, f_b8);
  // Calibration error of the cal-A combination: the Am position is shared
  // by every run, each run's Pb positions are its own.
  Double_t sw = 0, sum_dam = 0, var_pb = 0, sum_calb = 0;
  for (size_t i = 0; i < e_a.size(); i++) {
    Double_t w = 1.0 / (f_a[i] * f_a[i]);
    sw += w;
    sum_dam += w * d_am[i];
    var_pb += w * w *
              (d_k2[i] * d_k2[i] * s_k2[i] * s_k2[i] +
               d_k1[i] * d_k1[i] * s_k1[i] * s_k1[i]);
  }
  Double_t cal_comb_a =
      TMath::Sqrt(TMath::Power(sum_dam * am.ch_err, 2) + var_pb) / sw;
  Double_t cal_comb_a_am = TMath::Abs(sum_dam * am.ch_err) / sw;
  Double_t sw_b = 0;
  for (size_t i = 0; i < e_b.size(); i++) {
    Double_t w = 1.0 / (f_b[i] * f_b[i]);
    sw_b += w;
    sum_calb += w * cal_b_err[i];
  }
  Double_t cal_comb_b = sum_calb / sw_b;
  Double_t bkg_comb_a = TMath::Abs(ca.e - ca_flat.e),
           bkg_comb_b = TMath::Abs(cb.e - cb_flat.e);
  Double_t tot_a = TMath::Sqrt(ca.fit * ca.fit + cal_comb_a * cal_comb_a +
                               bkg_comb_a * bkg_comb_a + scan.rms * scan.rms);
  Double_t tot_b = TMath::Sqrt(cb.fit * cb.fit + cal_comb_b * cal_comb_b +
                               bkg_comb_b * bkg_comb_b + scan.rms * scan.rms);
  std::cout << "SIGE COMB per-run (" << e_a.size()
            << " runs) cal A: " << std::setprecision(4) << ca.e << " keV +- "
            << std::setprecision(1) << 1000 * ca.fit << " fit +- "
            << 1000 * cal_comb_a << " cal (Am part " << 1000 * cal_comb_a_am
            << ") +- " << 1000 * bkg_comb_a << " bkg +- " << 1000 * scan.rms
            << " win = " << 1000 * tot_a << " eV  chi2/ndf " << ca.chi2 << "/"
            << ca.ndf << "  p " << std::setprecision(3)
            << TMath::Prob(ca.chi2, ca.ndf) << std::endl;
  std::cout << "SIGE COMB per-run (" << e_b.size()
            << " runs) cal B: " << std::setprecision(4) << cb.e << " keV +- "
            << std::setprecision(1) << 1000 * cb.fit << " fit +- "
            << 1000 * cal_comb_b << " cal +- " << 1000 * bkg_comb_b
            << " bkg +- " << 1000 * scan.rms << " win = " << 1000 * tot_b
            << " eV  chi2/ndf " << cb.chi2 << "/" << cb.ndf << "  p "
            << std::setprecision(3) << TMath::Prob(cb.chi2, cb.ndf)
            << std::endl;
  std::cout << "SIGE COMB per-run sigma-fixed cal A: " << std::setprecision(4)
            << ca_fix.e << " +- " << std::setprecision(1) << 1000 * ca_fix.fit
            << " eV  chi2/ndf " << ca_fix.chi2 << "/" << ca_fix.ndf
            << std::endl;
  std::cout << "SIGE COMB per-run incl. BeamTurnOn_4 (" << e_a8.size()
            << " runs): cal A " << std::setprecision(4) << ca8.e << " +- "
            << std::setprecision(1) << 1000 * ca8.fit << " eV chi2/ndf "
            << ca8.chi2 << "/" << ca8.ndf << " | cal B " << std::setprecision(4)
            << cb8.e << " +- " << std::setprecision(1) << 1000 * cb8.fit
            << " eV chi2/ndf " << cb8.chi2 << "/" << cb8.ndf << std::endl;

  // ----------------------------------------------------------- per day
  std::vector<Double_t> day_e_a, day_f_a, day_e_b, day_f_b;
  for (std::map<Int_t, TH1D *>::iterator it = h_day.begin(); it != h_day.end();
       ++it) {
    Int_t d = it->first;
    TH1D *h = it->second;
    Point k2, k1;
    FitPb(h, nominal, k2, k1);
    Bool_t own_sources = h_am_day.count(d) > 0;
    Point am_d = own_sources ? FitAm(h_am_day[d], nominal) : am;
    std::vector<Point> pa;
    pa.push_back(am_d);
    pa.push_back(k2);
    pa.push_back(k1);
    Calibration cal_ad = Calibrate(pa);
    Calibration cal_bd = cal_b;
    if (own_sources) {
      Point b53 = FitBa53(h_ba_day[d], nominal), b80, b81;
      FitBa8081(h_ba_day[d], nominal, b80, b81);
      std::vector<Point> pb;
      pb.push_back(am_d);
      pb.push_back(b53);
      pb.push_back(b80);
      pb.push_back(b81);
      cal_bd = Calibrate(pb);
      std::cout << "SIGE DAY " << d << " Am " << std::setprecision(2) << am_d.ch
                << " +- " << am_d.ch_err << "  Ba53 " << b53.ch << " Ba80 "
                << b80.ch << " Ba81 " << b81.ch << " +- " << b81.ch_err
                << "  Pb Ka2 " << k2.ch << " Ka1 " << k1.ch << " +- "
                << k1.ch_err << std::endl;
    }
    GeFit l = FitGe(h, cal_ad, WIN_LO, WIN_HI, 1);
    GeFit fl = FitGe(h, cal_ad, WIN_LO, WIN_HI, 0);
    GeEnergy ra = Energy(l, fl, cal_ad, 0), rb = Energy(l, fl, cal_bd, 0);
    PrintEnergy(Form("DAY 10/%02d cal A%s", d,
                     own_sources ? " (own Am)" : " (all-days Am)"),
                ra, l, cal_ad);
    PrintEnergy(Form("DAY 10/%02d cal B%s", d,
                     own_sources ? " (own Am+Ba)" : " (global)"),
                rb, l, cal_bd);
    GeFit lfix = FitGe(h, cal_ad, WIN_LO, WIN_HI, 1, sigma_ge_ch);
    std::cout << "SIGE DAY 10/" << d << " sigma-fixed: cal A "
              << std::setprecision(4) << cal_ad.Energy(lfix.mu) << " cal B "
              << cal_bd.Energy(lfix.mu) << " +- " << std::setprecision(1)
              << 1000 * cal_ad.c1 * lfix.mu_err << " eV  area "
              << std::setprecision(0) << lfix.area << " +- " << lfix.area_err
              << std::endl;
    if (d == 9)
      continue;
    day_e_a.push_back(ra.e);
    day_f_a.push_back(TMath::Sqrt(ra.fit * ra.fit + ra.cal * ra.cal));
    day_e_b.push_back(rb.e);
    day_f_b.push_back(TMath::Sqrt(rb.fit * rb.fit + rb.cal * rb.cal));
  }
  Combined da = WeightedMean(day_e_a, day_f_a),
           db = WeightedMean(day_e_b, day_f_b);
  std::cout << "SIGE DAYS (10-12) cal A: mean " << std::setprecision(4) << da.e
            << " +- " << std::setprecision(1) << 1000 * da.fit
            << " eV chi2/ndf " << da.chi2 << "/" << da.ndf
            << " | cal B: " << std::setprecision(4) << db.e << " +- "
            << std::setprecision(1) << 1000 * db.fit << " eV chi2/ndf "
            << db.chi2 << "/" << db.ndf << std::endl;

  // ------------------------------------------------- line stability
  {
    std::vector<TString> am_runs = AmRuns(0), ba_runs = BaRuns(0);
    std::vector<Double_t> am_ch, am_err, ba_ch, ba_err;
    for (size_t i = 0; i < am_runs.size(); i++) {
      TH1D *h = NewSpectrum("stability_am");
      FillSpectrum(h, am_runs[i], estimator, Selection());
      Point p = FitAm(h, nominal);
      am_ch.push_back(p.ch);
      am_err.push_back(p.ch_err);
      delete h;
    }
    for (size_t i = 0; i < ba_runs.size(); i++) {
      TH1D *h = NewSpectrum("stability_ba");
      FillSpectrum(h, ba_runs[i], estimator, Selection());
      Point p80, p81;
      FitBa8081(h, nominal, p80, p81);
      ba_ch.push_back(p81.ch);
      ba_err.push_back(p81.ch_err);
      delete h;
    }
    Double_t ev_per_ch = 1000 * cal_a.c1;
    PrintStability("Am59", am_ch, am_err, ev_per_ch);
    PrintStability("Ba81", ba_ch, ba_err, ev_per_ch);
    PrintStability("PbKa1", pb_ka1_ch, pb_ka1_err, ev_per_ch);
  }

  // ------------------------------------------------------------- output
  DrawCalibration(cal_a, cal_b, all_lines);
  DrawSumFit(h_sum, cal_a, lin, flat, sum_a,
             EstimatorName(estimator) == "trapezoid"
                 ? "Ge runs, trapezoid, Am + Pb K#alpha calibration"
                 : "Ge runs, digitiser energy, Am + Pb K#alpha calibration");
  DrawPerRun(plot_values, 1000 * sigma_ge_kev, sum_a.e, sum_a.total, ca.e,
             ca_fix.e);
  WriteResult(project_root + "/results/si_ge_" + EstimatorName(estimator) +
                  ".result",
              estimator, sum_a, sum_b, sum_c, db);

  TFile out(SpectraFile(estimator), "RECREATE");
  TVectorD calibration(3);
  calibration[0] = cal_a.c0;
  calibration[1] = cal_a.c1;
  calibration[2] = cal_a.pivot;
  calibration.Write("cal_a");
  h_sum->Write();
  h_sum_all->Write();
  h_am->Write();
  h_ba->Write();
  for (size_t i = 0; i < h_run.size(); i++)
    h_run[i]->Write();
  out.Close();
}

// ---------------------------------------------------------------- compare

// The two energies' summed Ge spectra, each on its own Am + Pb Ka
// calibration, on one figure.
void Compare() {
  const Estimator estimators[2] = {kDigitizer, kTrapezoid};
  const Int_t colors[2] = {kGray + 2, kRed + 1};
  const char *labels[2] = {"Digitiser energy", "Trapezoid"};
  const Int_t rebin = 2;
  const Double_t show_lo = 56.0, show_hi = 77.0;
  TH1D *spectra[2] = {nullptr, nullptr};
  Double_t bin_ev = 0;
  for (Int_t k = 0; k < 2; k++) {
    TFile in(SpectraFile(estimators[k]));
    TH1D *h = static_cast<TH1D *>(in.Get("ge_sum"));
    TVectorD *calibration = static_cast<TVectorD *>(in.Get("cal_a"));
    if (!h || !calibration) {
      std::cout << "SIGE COMPARE missing " << SpectraFile(estimators[k])
                << std::endl;
      return;
    }
    h->SetDirectory(nullptr);
    Calibration cal;
    cal.c0 = (*calibration)[0];
    cal.c1 = (*calibration)[1];
    cal.pivot = (*calibration)[2];
    spectra[k] = ToKeV(h, cal, rebin, show_lo, show_hi, Form("compare_%d", k));
    bin_ev = 1000.0 * cal.c1 * rebin;
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  for (Int_t k = 0; k < 2; k++) {
    PlottingUtils::ConfigureHistogram(
        spectra[k], colors[k], Form(";Energy [keV];Counts / %.0f eV", bin_ev));
    spectra[k]->SetFillStyle(0);
  }
  spectra[0]->SetMinimum(
      0.8 * TMath::Min(spectra[0]->GetMinimum(), spectra[1]->GetMinimum()));
  spectra[0]->SetMaximum(
      1.25 * TMath::Max(spectra[0]->GetMaximum(), spectra[1]->GetMaximum()));
  spectra[0]->GetYaxis()->SetMaxDigits(3);
  spectra[0]->Draw("HIST");
  spectra[1]->Draw("HIST SAME");
  TLegend *legend = PlottingUtils::AddLegend(0.2, 0.5, 0.6, 0.78);
  for (Int_t k = 0; k < 2; k++)
    legend->AddEntry(spectra[k], labels[k], "l");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_si_ge_energies", "",
                            PlotSaveOptions::kLINEAR);
}

// ------------------------------------------------------------- bootstrap

// One record of the summed main Ge runs in both energies (digitiser
// channels) and its rise time.
struct RecordPair {
  Float_t trapezoid, digitizer, rise_ns;
};

std::vector<RecordPair> LoadPairs(const std::vector<TString> &runs) {
  std::vector<RecordPair> pairs;
  for (size_t r = 0; r < runs.size(); r++) {
    TFile file(FeaturesFile(runs[r]));
    TTree *tree = static_cast<TTree *>(file.Get("features"));
    Double_t digitizer_energy = 0;
    Float_t trapezoid_height = 0, baseline_rms = 0, rise_time = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("digitizer_energy", 1);
    tree->SetBranchStatus("trapezoid_height", 1);
    tree->SetBranchStatus("baseline_rms", 1);
    tree->SetBranchStatus("rise_time", 1);
    tree->SetBranchAddress("digitizer_energy", &digitizer_energy);
    tree->SetBranchAddress("trapezoid_height", &trapezoid_height);
    tree->SetBranchAddress("baseline_rms", &baseline_rms);
    tree->SetBranchAddress("rise_time", &rise_time);
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (baseline_rms < FLOOR_BASELINE_RMS)
        continue;
      Float_t trapezoid = trapezoid_height / TRAPEZOID_PER_CHANNEL;
      Bool_t in_range = (trapezoid > 1550 && trapezoid < 1900) ||
                        (digitizer_energy > 1550 && digitizer_energy < 1900);
      if (!in_range)
        continue;
      RecordPair pair = {trapezoid, Float_t(digitizer_energy),
                         Float_t(SAMPLE_NS * rise_time)};
      pairs.push_back(pair);
    }
  }
  return pairs;
}

// Ge on the Am + Pb Ka calibration of a summed spectrum, Am fixed.
Double_t GeOnOwnPb(TH1D *h, const Point &am) {
  LineSetting nominal;
  Point ka2, ka1;
  FitPb(h, nominal, ka2, ka1);
  std::vector<Point> points;
  points.push_back(am);
  points.push_back(ka2);
  points.push_back(ka1);
  Calibration cal = Calibrate(points);
  return cal.Energy(FitGe(h, cal, WIN_LO, WIN_HI, 1).mu);
}

// The two energies of the same records are two estimates that scatter
// against each other: every replica gives each record of the summed main
// runs a Poisson(1) weight, the same in both energies, and refits Pb Ka and
// Ge in each (Am fixed at its nominal fit). The rms of trapezoid - digitiser
// over the replicas is the scale against which their observed difference is
// judged, for all records and for the slow rise band.
void Bootstrap(Int_t n_replicas, const TString &project_root) {
  std::vector<GeRun> runs = GeRuns();
  std::vector<TString> main_runs;
  for (size_t i = 0; i < runs.size(); i++)
    if (runs[i].main)
      main_runs.push_back(runs[i].run);
  std::vector<RecordPair> pairs = LoadPairs(main_runs);
  std::cout << "SIGE BOOT " << pairs.size() << " records near the lines"
            << std::endl;
  // Per record, in the Ge window: the two energies' difference, its median
  // and its robust width (half the 16-84% range).
  {
    std::vector<Double_t> d;
    for (size_t i = 0; i < pairs.size(); i++)
      if (pairs[i].digitizer > 1594 && pairs[i].digitizer < 1702)
        d.push_back(pairs[i].trapezoid - pairs[i].digitizer);
    std::sort(d.begin(), d.end());
    Double_t ev_per_ch = 1000.0 / DPP_PER_KEV;
    Double_t median = d[d.size() / 2];
    Double_t width =
        0.5 * (d[Long64_t(0.84 * d.size())] - d[Long64_t(0.16 * d.size())]);
    std::cout << "SIGE BOOT per record in the Ge window (" << d.size()
              << "): trapezoid - digitiser median " << std::fixed
              << std::setprecision(1) << median * ev_per_ch
              << " eV, robust width " << width * ev_per_ch << " eV"
              << std::endl;
  }
  LineSetting nominal;
  Point am[2];
  for (Int_t e = 0; e < 2; e++) {
    TH1D *h_am =
        SumSpectrum(AmRuns(0), "boot_am", e == 0 ? kTrapezoid : kDigitizer);
    am[e] = FitAm(h_am, nominal);
    delete h_am;
  }
  const Int_t n_sets = 2;
  const char *set_names[n_sets] = {"all", "rise 280-360 ns"};
  TRandom3 random(20261024);
  std::vector<Double_t> difference[n_sets];
  Double_t observed[n_sets] = {0, 0};
  for (Int_t replica = -1; replica < n_replicas; replica++) {
    TH1D *h[n_sets][2];
    for (Int_t s = 0; s < n_sets; s++)
      for (Int_t e = 0; e < 2; e++)
        h[s][e] = NewSpectrum(Form("boot_%d_%d", s, e));
    for (size_t i = 0; i < pairs.size(); i++) {
      Int_t weight = replica < 0 ? 1 : random.Poisson(1.0);
      if (weight == 0)
        continue;
      Bool_t slow = pairs[i].rise_ns >= 280 && pairs[i].rise_ns < 360;
      for (Int_t s = 0; s < n_sets; s++) {
        if (s == 1 && !slow)
          continue;
        h[s][0]->Fill(pairs[i].trapezoid, weight);
        h[s][1]->Fill(pairs[i].digitizer, weight);
      }
    }
    for (Int_t s = 0; s < n_sets; s++) {
      Double_t d = GeOnOwnPb(h[s][0], am[0]) - GeOnOwnPb(h[s][1], am[1]);
      if (replica < 0)
        observed[s] = d;
      else
        difference[s].push_back(d);
      for (Int_t e = 0; e < 2; e++)
        delete h[s][e];
    }
  }
  for (Int_t s = 0; s < n_sets; s++) {
    Double_t mean = TMath::Mean(difference[s].size(), difference[s].data());
    Double_t rms = TMath::RMS(difference[s].size(), difference[s].data());
    std::cout << "SIGE BOOT " << std::left << std::setw(16) << set_names[s]
              << std::right << " trapezoid - digitiser " << std::showpos
              << std::fixed << std::setprecision(1) << 1000 * observed[s]
              << std::noshowpos << " +- " << 1000 * rms << " eV ("
              << std::showpos << observed[s] / rms << std::noshowpos
              << " sigma from agreement; " << n_replicas << " replicas, mean "
              << std::showpos << 1000 * mean << std::noshowpos << " eV)"
              << std::endl;
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TH1D *h_all = new TH1D("boot_difference", "", 30, -45, 45);
  TH1D *h_slow = new TH1D("boot_difference_slow", "", 30, -45, 45);
  h_all->SetDirectory(nullptr);
  h_slow->SetDirectory(nullptr);
  for (Int_t i = 0; i < n_replicas; i++) {
    h_all->Fill(1000 * difference[0][i]);
    h_slow->Fill(1000 * difference[1][i]);
  }
  PlottingUtils::ConfigureHistogram(
      h_all, kGray + 2, ";Trapezoid #minus digitiser Ge [eV];Replicas / 3 eV");
  PlottingUtils::ConfigureHistogram(h_slow, kAzure + 2, "");
  h_all->SetFillStyle(0);
  h_slow->SetFillStyle(0);
  Double_t y_top = TMath::Max(h_all->GetMaximum(), h_slow->GetMaximum());
  h_all->SetMaximum(1.9 * y_top);
  h_all->Draw("HIST");
  h_slow->Draw("HIST SAME");
  TLine *line_all = new TLine(1000 * observed[0], 0, 1000 * observed[0], y_top);
  TLine *line_slow =
      new TLine(1000 * observed[1], 0, 1000 * observed[1], y_top);
  line_all->SetLineColor(kGray + 2);
  line_slow->SetLineColor(kAzure + 2);
  line_all->SetLineWidth(3);
  line_slow->SetLineWidth(3);
  line_slow->SetLineStyle(2);
  line_all->Draw();
  line_slow->Draw();
  TLine *agreement = new TLine(0, 0, 0, y_top);
  agreement->SetLineColor(kBlack);
  agreement->SetLineStyle(3);
  agreement->SetLineWidth(2);
  agreement->Draw();
  TLegend *legend = PlottingUtils::AddLegend(0.16, 0.62, 0.64, 0.9);
  legend->AddEntry(h_all, "All records, replicas", "l");
  legend->AddEntry(line_all, "All records, observed", "l");
  legend->AddEntry(h_slow, "Rise 280-360 ns, replicas", "l");
  legend->AddEntry(line_slow, "Rise 280-360 ns, observed", "l");
  legend->AddEntry(agreement, "The two energies agree", "l");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_si_energy_bootstrap", "",
                            PlotSaveOptions::kLINEAR);
}

} // namespace SiGe

void SiGeBudget(const char *energy = "trapezoid", Int_t n_replicas = 200) {
  using namespace SiGe;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString mode = energy;
  Bool_t both = mode == "compare" || mode == "bootstrap";
  TString plot_dir = project_root + "/plots/studies/psi_si/" +
                     (both ? TString("ge_compare") : "ge_" + mode);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, plot_dir, SPECTRA_DIR);
  if (mode == "compare")
    Compare();
  else if (mode == "bootstrap")
    Bootstrap(n_replicas, project_root);
  else
    Budget(EstimatorFromName(mode), project_root);
}
