// The 68.75 keV 72Ge(n,g) line in the PSI Amptek Si, from the trapezoid of
// the processed waveforms (ProcessWaveforms.cpp; the digitiser's own on-board
// energy is not used). Done the collaborator's way: every Ge beam run of the
// first PSI run (10/09-10/12) in its six day/night groups, each a Gaussian
// on a linear background over 66.5-71 keV, calibrated with Am-241 59.5409
// keV (source runs) and the group's own Pb Ka2/Ka1 (cal A); the result is
// the groups' error-weighted mean. The summed spectrum is a check.
// Alternatives:
// Am-241 and Ba-133 53.1622, 79.6142, 80.9979 keV (cal B, source lines
// only) and all six lines (cal C).
// Uncertainties: the Ge centroid fit, the calibration offset and slope (the
// pivot at the calibration points' weighted centroid, where the two are
// uncorrelated), the background (|linear - flat|) and the fit window (RMS
// over a window scan). Checks: per run, per day, the width fixed, the 53.44
// keV line, calibration variants, isolated events, and the run-to-run
// stability of the line positions.
//
// Writes results/si_ge.result, the spectra to
// /labdata/ANSG/PSI/processed/Si/ge_spectra.root, and the figures to
// plots/si/ge_budget.
//
// Usage: root -b -q 'SiGeBudget.cpp+'
#include "FigureStyle.hpp"
#include "SiAnalysis.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
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
  Bool_t converged = kFALSE;
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
  p.converged = f.converged;
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
  Bool_t converged = kFALSE;
  // A physical width and an area above twice its error.
  Bool_t plausible_line = kFALSE;
  Double_t mu = 0, mu_err = 0, sigma = 0, sigma_err = 0, area = 0, area_err = 0;
  Double_t chi2 = 0;
  Int_t ndf = 0;
  Double_t lo = 0, hi = 0;
  Double_t c0 = 0, c1 = 0, c2 = 0;
  Int_t status = -1;
};

GeFit FromOutput(const FitOutput &f) {
  GeFit g;
  g.status = f.status;
  g.converged = f.converged;
  g.plausible_line =
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
            << lin.area << " +- " << lin.area_err << "  "
            << (lin.converged
                    ? "converged"
                    : TString::Format("NOT CONVERGED (status %d)", lin.status))
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
                const GeFit &flat, const TString &label) {
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
  TGraph *pulls = new TGraph();
  for (Int_t b = 1; b <= data->GetNbinsX(); b++) {
    Double_t ch_lo = cal.Channel(data->GetBinLowEdge(b));
    if (ch_lo < lin.lo - 0.01 || ch_lo + rebin > lin.hi + 1.01)
      continue;
    Double_t expected = 0;
    for (Int_t k = 0; k < rebin; k++)
      expected += ModelPerChannel(lin, ch_lo + 0.5 + k);
    pulls->SetPoint(pulls->GetN(), data->GetBinCenter(b),
                    (data->GetBinContent(b) - expected) /
                        TMath::Sqrt(expected));
  }
  std::vector<TString> labels;
  labels.push_back(label);
  labels.push_back(Form("Gaussian + linear background, #chi^{2}/ndf = %.3f",
                        lin.chi2 / lin.ndf));
  labels.push_back("Gaussian + flat background");
  labels.push_back("Linear background");
  const Int_t canvas_w = 1400, plot_h = 800, panel_h = 300;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(canvas_w, plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(canvas, panel_h, upper, lower);
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
  TLegend *leg = FigureStyle::AddFrameLegend(upper, labels, kFALSE);
  leg->AddEntry(data, labels[0], "lep");
  leg->AddEntry(g_lin, labels[1], "l");
  leg->AddEntry(g_flat, labels[2], "l");
  leg->AddEntry(g_background, labels[3], "l");
  lower->cd();
  FigureStyle::DrawPullPanel(pulls, data->GetXaxis()->GetXmin(),
                             data->GetXaxis()->GetXmax(), "Energy [keV]");
  PlottingUtils::ScaleFigure(canvas, canvas_w, plot_h);
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
  TLine *ge = new TLine(ReferenceLines::GE73_ENSDF, y_lo,
                        ReferenceLines::GE73_ENSDF, 0.5 * y_hi);
  ge->SetLineColor(kGreen + 3);
  ge->SetLineStyle(7);
  ge->Draw();
  TLatex *t = new TLatex(ReferenceLines::GE73_ENSDF + 0.4, y_lo + 10, "Ge");
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

struct GroupValue {
  TString label;
  Double_t e, err, draft_e;
};

// The six day/night groups the result averages (Am + each group's own Pb
// Ka calibration, total errors) beside the collaborator's draft values, with
// the result and the draft's Si-1 mean.
void DrawGroups(const std::vector<GroupValue> &values, Double_t result_e,
                Double_t result_err, Double_t draft_mean) {
  Int_t n = (Int_t)values.size();
  std::vector<TString> labels;
  labels.push_back("This analysis, total error");
  labels.push_back("Collaborator's draft");
  labels.push_back("Weighted mean, total error");
  labels.push_back("Draft's Si-1 mean");
  Double_t legend_w = PlottingUtils::LegendWidthPx(labels);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)legend_w);
  plot->cd();
  plot->SetBottomMargin(0.27);
  const Double_t y_lo = 68.60, y_hi = 68.95;
  TH1D *frame = new TH1D("groups_frame", "", n, 0.5, n + 0.5);
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
  TBox *band =
      new TBox(0.5, result_e - result_err, n + 0.5, result_e + result_err);
  band->SetFillColorAlpha(kGray, 0.6);
  band->SetLineWidth(0);
  band->Draw();
  TLine *result_line = new TLine(0.5, result_e, n + 0.5, result_e);
  result_line->SetLineColor(kGray + 2);
  result_line->SetLineWidth(2);
  result_line->Draw();
  TLine *draft_line = new TLine(0.5, draft_mean, n + 0.5, draft_mean);
  draft_line->SetLineColor(kRed + 1);
  draft_line->SetLineStyle(2);
  draft_line->SetLineWidth(2);
  draft_line->Draw();
  TGraphErrors *g_ours = new TGraphErrors(n);
  TGraph *g_draft = new TGraph(n);
  for (Int_t i = 0; i < n; i++) {
    g_ours->SetPoint(i, i + 1 - 0.1, values[i].e);
    g_ours->SetPointError(i, 0, values[i].err);
    g_draft->SetPoint(i, i + 1 + 0.1, values[i].draft_e);
  }
  PlottingUtils::ConfigureGraph(g_ours, kBlack);
  PlottingUtils::ConfigureGraph(g_draft, kRed + 1);
  g_draft->SetMarkerStyle(24);
  g_draft->SetMarkerColor(kRed + 1);
  g_draft->SetMarkerSize(1.6);
  g_ours->Draw("P SAME");
  g_draft->Draw("P SAME");
  frame->Draw("AXIS SAME");
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  leg->AddEntry(g_ours, labels[0], "lep");
  leg->AddEntry(g_draft, labels[1], "p");
  leg->AddEntry(band, labels[2], "f");
  leg->AddEntry(draft_line, labels[3], "l");
  PlottingUtils::SaveFigure(canvas, "psi_si_ge_groups", "",
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

TString SpectraFile() { return TString(SPECTRA_DIR) + "/ge_spectra.root"; }

void WriteResult(const TString &path, const GeEnergy &result,
                 const GeEnergy &summed, const GeEnergy &cal_b,
                 const GeEnergy &cal_c, const Combined &days_b) {
  std::ofstream out(path.Data());
  out << "# Si (PSI) Ge energy, SiGeBudget.cpp, trapezoid; the six day/night "
         "groups of 10/09-10/12, Am + Pb Ka calibration; keV"
      << std::endl;
  out << std::fixed << std::setprecision(6);
  out << "ge " << result.e << std::endl;
  out << "total " << result.total << std::endl;
  out << "fit " << result.fit << std::endl;
  out << "calibration " << result.cal << std::endl;
  out << "background " << result.background << std::endl;
  out << "window " << result.window << std::endl;
  out << "# alternatives" << std::endl;
  out << "ge_summed_main_runs " << summed.e << std::endl;
  out << "ge_cal_b " << cal_b.e << std::endl;
  out << "ge_cal_c " << cal_c.e << std::endl;
  out << "ge_days_cal_b " << days_b.e << std::endl;
  out << "ge_days_cal_b_err " << days_b.fit << std::endl;
}

// ----------------------------------------------------------------- budget

// The result, the collaborator's Si-1 the way the draft describes it
// (main.tex): six day/night groups of the first PSI run, each calibrated
// with Am-241 and its own Pb K-alpha, its Ge line a Gaussian on a linear
// background with uncertainties from the centroid, the calibration's offset
// and slope, and linear against flat background; then their error-weighted
// mean, the Am position shared by every group's calibration. Printed beside
// the draft's own per-group values.
const Int_t N_GROUPS = 6;
const char *const GROUP_NAMES[N_GROUPS] = {"Thu 10/09",     "Fri night 10/10",
                                           "Sat day 10/11", "Sat night 10/11",
                                           "Sun day 10/12", "Sun night 10/12"};
// Each group's runs, by GeRuns() label.
const char *const GROUP_RUNS[N_GROUPS][2] = {
    {"BeamTurnOn_4", ""},    {"GeTarget_1", "GeTarget_2"},
    {"SiYesGe_1056", ""},    {"Ge_1715_1", "Ge_1715_2"},
    {"LongWf500V_1001", ""}, {"LongWf_1517", ""}};
// The draft's values (its comments to the Si-1 table), keV.
const Double_t DRAFT_GROUP_E[N_GROUPS] = {68.78295, 68.78285, 68.73901,
                                          68.74103, 68.72764, 68.75110};
const Double_t DRAFT_MEAN = 68.7543, DRAFT_MEAN_ERR = 0.0115;

GeEnergy CollaboratorGroups(const std::vector<GeRun> &runs,
                            const std::vector<TH1D *> &h_run, const Point &am,
                            const LineSetting &nominal, Double_t window_err,
                            std::vector<GroupValue> &group_values) {
  std::vector<Double_t> e, e_flat, total_err, fit_err;
  std::vector<Double_t> d_am, pb_var;
  for (Int_t g = 0; g < N_GROUPS; g++) {
    TH1D *h = NewSpectrum(Form("collab_group_%d", g));
    for (size_t i = 0; i < runs.size(); i++)
      for (Int_t k = 0; k < 2; k++)
        if (TString(GROUP_RUNS[g][k]) != "" &&
            runs[i].label == GROUP_RUNS[g][k])
          h->Add(h_run[i]);
    Point k2, k1;
    FitPb(h, nominal, k2, k1);
    std::vector<Point> p;
    p.push_back(am);
    p.push_back(k2);
    p.push_back(k1);
    Calibration cal = Calibrate(p);
    GeFit lin = FitGe(h, cal, WIN_LO, WIN_HI, 1);
    GeFit flat = FitGe(h, cal, WIN_LO, WIN_HI, 0);
    GeEnergy r = Energy(lin, flat, cal, 0);
    PrintEnergy(Form("GROUP %s", GROUP_NAMES[g]), r, lin, cal);
    // The energy's derivatives with respect to the calibration lines'
    // positions: Am's is shared by the groups, the Pb lines' are each
    // group's own.
    const Double_t delta = 0.1;
    Double_t derivative[3];
    for (Int_t k = 0; k < 3; k++) {
      std::vector<Point> q = p;
      q[k].ch += delta;
      derivative[k] = (Calibrate(q).Energy(lin.mu) - r.e) / delta;
    }
    d_am.push_back(derivative[0]);
    pb_var.push_back(TMath::Power(derivative[1] * k2.ch_err, 2) +
                     TMath::Power(derivative[2] * k1.ch_err, 2));
    e_flat.push_back(cal.Energy(flat.mu));
    std::cout << "SIGE GROUP " << std::left << std::setw(16) << GROUP_NAMES[g]
              << std::right << " ours " << std::fixed << std::setprecision(4)
              << r.e << " +- " << std::setprecision(1) << 1000 * r.total
              << " eV, draft " << std::setprecision(4) << DRAFT_GROUP_E[g]
              << ": ours - draft " << std::showpos << std::setprecision(1)
              << 1000 * (r.e - DRAFT_GROUP_E[g]) << std::noshowpos << " eV"
              << std::endl;
    e.push_back(r.e);
    total_err.push_back(r.total);
    fit_err.push_back(r.fit);
    GroupValue value = {GROUP_NAMES[g], r.e, r.total, DRAFT_GROUP_E[g]};
    group_values.push_back(value);
    delete h;
  }
  Combined by_total = WeightedMean(e, total_err);
  Combined by_fit = WeightedMean(e, fit_err);
  GeEnergy result;
  Double_t sw = 0, sum_flat = 0, var_fit = 0, sum_dam = 0, var_pb = 0;
  for (size_t g = 0; g < e.size(); g++) {
    Double_t w = 1.0 / (total_err[g] * total_err[g]);
    sw += w;
    sum_flat += w * e_flat[g];
    var_fit += w * w * fit_err[g] * fit_err[g];
    sum_dam += w * d_am[g];
    var_pb += w * w * pb_var[g];
  }
  result.e = by_total.e;
  result.fit = TMath::Sqrt(var_fit) / sw;
  result.cal = TMath::Sqrt(TMath::Power(sum_dam * am.ch_err, 2) + var_pb) / sw;
  result.background = TMath::Abs(result.e - sum_flat / sw);
  result.window = window_err;
  result.total = TMath::Sqrt(result.fit * result.fit + result.cal * result.cal +
                             result.background * result.background +
                             result.window * result.window);
  std::cout << "SIGE GROUPS error-weighted mean: " << std::setprecision(4)
            << by_total.e << " +- " << std::setprecision(1)
            << 1000 * by_total.fit << " eV (total errors), "
            << std::setprecision(4) << by_fit.e << " +- "
            << std::setprecision(1) << 1000 * by_fit.fit
            << " eV (fit errors); draft " << std::setprecision(4) << DRAFT_MEAN
            << " +- " << std::setprecision(1) << 1000 * DRAFT_MEAN_ERR << " eV"
            << std::endl;
  std::cout << "SIGE RESULT six groups: " << std::setprecision(4) << result.e
            << " keV +- " << std::setprecision(1) << 1000 * result.fit
            << " fit +- " << 1000 * result.cal << " cal +- "
            << 1000 * result.background << " bkg +- " << 1000 * result.window
            << " win = " << 1000 * result.total << " eV, chi2/ndf "
            << by_total.chi2 << "/" << by_total.ndf << std::endl;
  return result;
}

void Budget(const TString &project_root) {
  // ------------------------------------------------------------ spectra
  std::vector<GeRun> runs = GeRuns();
  std::vector<TH1D *> h_run;
  TH1D *h_sum = NewSpectrum("ge_sum");
  TH1D *h_sum_all = NewSpectrum("ge_sum_all");
  std::map<Int_t, TH1D *> h_day;
  Double_t seconds_main = 0;
  for (size_t i = 0; i < runs.size(); i++) {
    TH1D *h = NewSpectrum("ge_" + runs[i].label);
    FillSpectrum(h, runs[i].run, Selection());
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
  TH1D *h_am = SumSpectrum(AmRuns(0), "am_all");
  TH1D *h_ba = SumSpectrum(BaRuns(0), "ba_all");
  std::map<Int_t, TH1D *> h_am_day, h_ba_day;
  for (Int_t d = 10; d <= 13; d++) {
    h_am_day[d] = SumSpectrum(AmRuns(d), Form("am_d%d", d));
    h_ba_day[d] = SumSpectrum(BaRuns(d), Form("ba_d%d", d));
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
            << cal_a.Energy(flat_fixed.mu) << std::endl;
  std::cout << "SIGE SUM chi2/ndf linear " << std::setprecision(1) << lin.chi2
            << "/" << lin.ndf << "  flat " << flat.chi2 << "/" << flat.ndf
            << "  quadratic " << quad.chi2 << "/" << quad.ndf << "  bkg slope "
            << std::setprecision(4) << lin.c1 / lin.c0 * 100 << " %/ch"
            << std::endl;

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
    TH1D *h_isolated = SumSpectrum(main_runs, "ge_sum_isolated", isolated);
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

  // The two rise bands alone, each on its own Am + Pb Ka calibration (the
  // trapezoid's gain differs by 0.1% between them).
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
      TH1D *h_band = SumSpectrum(main_runs, "ge_sum_band", band);
      TH1D *h_am_band = SumSpectrum(AmRuns(0), "am_band", band);
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
    TString tag = runs[i].label + (l.plausible_line ? "" : " [weak fit]");
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

  std::vector<GroupValue> group_values;
  GeEnergy groups =
      CollaboratorGroups(runs, h_run, am, nominal, scan.rms, group_values);

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
      FillSpectrum(h, am_runs[i], Selection());
      Point p = FitAm(h, nominal);
      delete h;
      if (!p.converged) {
        std::cout << "SIGE STABILITY skips " << am_runs[i]
                  << ": its Am fit did not converge" << std::endl;
        continue;
      }
      am_ch.push_back(p.ch);
      am_err.push_back(p.ch_err);
    }
    for (size_t i = 0; i < ba_runs.size(); i++) {
      TH1D *h = NewSpectrum("stability_ba");
      FillSpectrum(h, ba_runs[i], Selection());
      Point p80, p81;
      FitBa8081(h, nominal, p80, p81);
      delete h;
      if (!p81.converged) {
        std::cout << "SIGE STABILITY skips " << ba_runs[i]
                  << ": its Ba fit did not converge" << std::endl;
        continue;
      }
      ba_ch.push_back(p81.ch);
      ba_err.push_back(p81.ch_err);
    }
    Double_t ev_per_ch = 1000 * cal_a.c1;
    PrintStability("Am59", am_ch, am_err, ev_per_ch);
    PrintStability("Ba81", ba_ch, ba_err, ev_per_ch);
    PrintStability("PbKa1", pb_ka1_ch, pb_ka1_err, ev_per_ch);
  }

  // ------------------------------------------------------------- output
  DrawCalibration(cal_a, cal_b, all_lines);
  DrawSumFit(h_sum, cal_a, lin, flat, "Ge runs 10/10#minus10/12, summed");
  DrawGroups(group_values, groups.e, groups.total, DRAFT_MEAN);
  WriteResult(project_root + "/results/si_ge.result", groups, sum_a, sum_b,
              sum_c, db);

  TFile out(SpectraFile(), "RECREATE");
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

  std::cout << "SIGE FITS not converged: " << UnconvergedFitCounter()
            << std::endl;
}

} // namespace SiGe

void SiGeBudget() {
  using namespace SiGe;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/si/ge_budget",
                                SPECTRA_DIR);
  Budget(project_root);
}
