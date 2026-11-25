// Whether the CZT's pixel gains are a fixed property of the crystals or
// change from day to day (the electronics might respond to the irradiations,
// shielded and thermal as they are). Three measures per day and pixel, in eV
// at Ge's energy:
//   Pb K-alpha  each beam run's pixel K-alpha1 position against the run's
//   95.902 keV  summed spectrum (PixelAlignment.hpp), on top of the production
//               pixel gains; the day's map is the weighted mean of its runs',
//               and each run is compared with it. Pb only behind the Cu or
//               Cd shield (01/13-14); the capture line every day. Both lines
//               from the filtered events (z >= 8 mm); their depth profiles
//               near the cathode, background-subtracted, are printed;
//   558 keV     the 558.456 keV capture line per pixel in the day's runs
//               together, selected as PixelCalibration.cpp does (every
//               interaction on a pixel centre, any depth, energies before the
//               pixel gains), as a gain against the production table (01/16).
// For pairs of days, the difference of the two maps against its errors (chi2
// per pixel, after removing each map's mean, which the calibration absorbs)
// and the correlation of the maps: a fixed pattern correlates strongly and
// differs only by its errors.
//
// Usage: root -b -q 'PixelGainStability.cpp+'. Plots go to
// plots/czt/studies/pixel_gain.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PixelAlignment.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLine.h>
#include <TStyle.h>

namespace CztPixelGain {

using namespace PixelAlignment;

const Double_t E_CD558 = 558.456;
const Double_t E_GE = 68.75;
// The 558 keV fit window, before the pixel gains (about 1 keV per unit).
const Double_t LO_558 = 535.0, HI_558 = 580.0;
const Int_t MIN_EVENTS_558 = 300;

struct Day {
  TString label, tag;
  std::vector<TString> runs;
};

std::vector<Day> Days() {
  std::vector<Day> days(4);
  days[0].label = "01/13";
  days[0].tag = "0113";
  days[0].runs.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);
  days[0].runs.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
  days[0].runs.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  days[0].runs.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  days[0].runs.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
  days[0].runs.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
  days[0].runs.push_back(Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113);
  days[0].runs.push_back(Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113);
  days[1].label = "01/14";
  days[1].tag = "0114";
  days[1].runs.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114);
  days[1].runs.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260114);
  days[1].runs.push_back(Constants::CUSHIELDSIGNAL_90PERCENT_20260114);
  days[2].label = "01/15";
  days[2].tag = "0115";
  days[2].runs.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115);
  days[2].runs.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_20260115);
  days[3].label = "01/16";
  days[3].tag = "0116";
  days[3].runs.push_back(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116);
  days[3].runs.push_back(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116);
  days[3].runs.push_back(Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116);
  days[3].runs.push_back(
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116);
  return days;
}

// A value per pixel; err <= 0 where the pixel is not measured.
struct PixelMap {
  std::vector<Double_t> value, err;
  PixelMap() : value(N_PIXELS, 0.0), err(N_PIXELS, 0.0) {}
};

Double_t WeightedMeanOf(const PixelMap &m) {
  Double_t sw = 0, swx = 0;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (m.err[p] <= 0)
      continue;
    sw += 1.0 / (m.err[p] * m.err[p]);
    swx += m.value[p] / (m.err[p] * m.err[p]);
  }
  return sw > 0 ? swx / sw : 0;
}

// ------------------------------------------------------------ line maps

PixelMap LineDayMap(const Day &day, Int_t line) {
  std::vector<PixelMap> run_maps;
  std::vector<TString> run_names;
  for (size_t r = 0; r < day.runs.size(); r++) {
    const RunAlignment &a = AlignmentOf(day.runs[r], line);
    if (!a.valid)
      continue;
    PixelMap m;
    for (Int_t p = 0; p < N_PIXELS; p++) {
      if (!a.measured[p])
        continue;
      m.value[p] = a.OffsetEv(p);
      m.err[p] = a.OffsetErrEv(p);
    }
    run_maps.push_back(m);
    run_names.push_back(day.runs[r]);
  }
  PixelMap day_map;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    Double_t sw = 0, swx = 0;
    for (size_t r = 0; r < run_maps.size(); r++) {
      if (run_maps[r].err[p] <= 0)
        continue;
      Double_t w = 1.0 / (run_maps[r].err[p] * run_maps[r].err[p]);
      sw += w;
      swx += w * run_maps[r].value[p];
    }
    if (sw > 0) {
      day_map.value[p] = swx / sw;
      day_map.err[p] = 1.0 / TMath::Sqrt(sw);
    }
  }
  // Each run against the day's map (which it is part of, hence the
  // difference's variance err_run^2 - err_day^2).
  for (size_t r = 0; r < run_maps.size(); r++) {
    Double_t chi2 = 0;
    Int_t n = 0;
    for (Int_t p = 0; p < N_PIXELS; p++) {
      if (run_maps[r].err[p] <= 0)
        continue;
      Double_t variance = run_maps[r].err[p] * run_maps[r].err[p] -
                          day_map.err[p] * day_map.err[p];
      if (variance <= 0)
        continue;
      chi2 +=
          TMath::Power(run_maps[r].value[p] - day_map.value[p], 2) / variance;
      n++;
    }
    std::cout << "LINE-DAY " << LINE_NAMES[line] << " " << day.label << " "
              << run_names[r] << ": against the day's map chi2/pixel "
              << std::fixed << std::setprecision(2) << (n ? chi2 / n : 0)
              << " (" << n << " pixels)" << std::endl;
  }
  return day_map;
}

// The depth profile of a line's events near the cathode (z = 10 mm): events
// within one sigma of the line, less those 3-5 sigma above it (clear of its
// low tail), over the day's runs.
void PrintLineDepth(const Day &day, Int_t line) {
  const Int_t n_bins = 8;
  TH1D *band = new TH1D("depth_band_" + day.tag + Form("_%d", line), "", n_bins,
                        8.0, 10.0);
  TH1D *side = new TH1D("depth_side_" + day.tag + Form("_%d", line), "", n_bins,
                        8.0, 10.0);
  band->SetDirectory(nullptr);
  side->SetDirectory(nullptr);
  for (size_t r = 0; r < day.runs.size(); r++) {
    const RunAlignment &a = AlignmentOf(day.runs[r], line);
    if (!a.valid)
      continue;
    TFile *file = IO::OpenForReading("filtered/" + day.runs[r] + ".root");
    if (!file || file->IsZombie())
      continue;
    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TTree *tree =
          static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
      if (!tree)
        continue;
      Float_t energy = 0, z = 0;
      tree->SetBranchAddress("energykeV", &energy);
      tree->SetBranchAddress("zmm", &z);
      Long64_t n = tree->GetEntries();
      for (Long64_t i = 0; i < n; i++) {
        tree->GetEntry(i);
        Double_t u = (energy - a.mu) / a.sigma;
        if (TMath::Abs(u) < 1.0)
          band->Fill(z);
        if (u > 3.0 && u < 5.0)
          side->Fill(z);
      }
      tree->ResetBranchAddresses();
    }
    file->Close();
    delete file;
  }
  band->Add(side, -1.0);
  Double_t total = band->Integral();
  Double_t within_half_mm = band->Integral(band->FindBin(9.51), n_bins);
  Double_t within_one_mm = band->Integral(band->FindBin(9.01), n_bins);
  std::cout << "DEPTH " << LINE_NAMES[line] << " " << day.label
            << ": background-subtracted events " << std::fixed
            << std::setprecision(0) << total << ", mean z "
            << std::setprecision(3) << band->GetMean() << " mm, within 0.5 mm "
            << "of the cathode " << within_half_mm / total << ", within 1 mm "
            << within_one_mm / total << std::endl;
  delete band;
  delete side;
}

// ---------------------------------------------------------------- 558 maps

// [0] position [1] sigma [2] core counts [3] tail amplitude [4] tail decay
// [5] continuum per bin [6] bin width
Double_t SingleLineModel(Double_t *x, Double_t *p) {
  return p[5] + Line(x[0], p[0], p[1], p[2], p[3], p[4], p[6]);
}

TH1D *New558Hist(const TString &name) {
  const Double_t bin = Constants::BIN_WIDTH_KEV;
  const Int_t nbins = (Int_t)std::lround((HI_558 - LO_558) / bin);
  TH1D *h = new TH1D(name, "", nbins, LO_558, LO_558 + nbins * bin);
  h->SetDirectory(nullptr);
  return h;
}

// The 558 keV line in h, its tail held unless tail_free; returns the fit.
TF1 *Fit558(TH1D *h, const TString &name, Bool_t tail_free, Double_t tail_amp,
            Double_t tail_decay, Bool_t &ok) {
  TF1 *f = new TF1(name, SingleLineModel, LO_558, HI_558, 7);
  f->FixParameter(6, Constants::BIN_WIDTH_KEV);
  h->GetXaxis()->SetRangeUser(550.0, 570.0);
  Double_t top = h->GetBinCenter(h->GetMaximumBin());
  Double_t height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  const Double_t sigma0 = 2.0;
  f->SetParameter(0, top);
  f->SetParLimits(0, top - 4.0, top + 4.0);
  f->SetParameter(1, sigma0);
  f->SetParLimits(1, 0.5, 6.0);
  f->SetParameter(2, height * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 /
                         Constants::BIN_WIDTH_KEV);
  f->SetParLimits(2, 0, 1e9);
  if (tail_free) {
    f->SetParameter(3, 0.1);
    f->SetParLimits(3, 0, 2.0);
    f->SetParameter(4, 3.0);
    f->SetParLimits(4, 0.5, 20.0);
  } else {
    f->FixParameter(3, tail_amp);
    f->FixParameter(4, tail_decay);
  }
  f->SetParameter(5, TMath::Max(h->GetBinContent(h->GetNbinsX()), 0.5));
  f->SetParLimits(5, 0, 1e9);
  h->Fit(f, "QLRN0");
  TFitResultPtr r = h->Fit(f, "QLRN0MS");
  ok = (Converged(r) && f->GetParError(0) > 0);
  return f;
}

std::vector<Double_t> ProductionGains(const TString &project_root) {
  std::vector<Double_t> gains(N_PIXELS, 0.0);
  TFile *file =
      TFile::Open(project_root + "/root_files/pixel_calibration.root");
  TH2D *table = file ? static_cast<TH2D *>(file->Get("pixel_gain")) : nullptr;
  if (!table) {
    std::cerr << "ERROR: no pixel_gain in root_files/pixel_calibration.root"
              << std::endl;
    return gains;
  }
  for (Int_t ix = 0; ix < N_SIDE; ix++)
    for (Int_t iy = 0; iy < N_SIDE; iy++)
      gains[ix * N_SIDE + iy] = table->GetBinContent(ix + 1, iy + 1);
  file->Close();
  delete file;
  return gains;
}

PixelMap Gain558DayMap(const Day &day, const std::vector<Double_t> &table) {
  TH1D *sum = New558Hist("h558_sum_" + day.tag);
  std::vector<TH1D *> pixels;
  for (Int_t p = 0; p < N_PIXELS; p++)
    pixels.push_back(New558Hist(Form("h558_%s_%d", day.tag.Data(), p)));
  for (size_t r = 0; r < day.runs.size(); r++) {
    TFile *file = IO::OpenForReading("filtered/" + day.runs[r] + ".root");
    TTree *tree = file ? static_cast<TTree *>(file->Get("bef_tree")) : nullptr;
    if (!tree) {
      std::cerr << "WARNING: no bef_tree in " << day.runs[r] << std::endl;
      continue;
    }
    Float_t energy = 0, x = 0, y = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("energykeV", 1);
    tree->SetBranchStatus("xmm", 1);
    tree->SetBranchStatus("ymm", 1);
    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
      tree->GetEntry(i);
      if (energy < LO_558 - 5.0 || energy > HI_558 + 5.0)
        continue;
      Int_t p = PixelIndexOf(x, y);
      if (p < 0 || table[p] <= 0)
        continue;
      pixels[p]->Fill(energy);
      sum->Fill(energy * table[p]);
    }
    file->Close();
    delete file;
  }
  // The shape from the day's spectrum on the production gains; each pixel's
  // line before them, with that tail held.
  Bool_t ok = kFALSE;
  TF1 *shape = Fit558(sum, "f558_sum_" + day.tag, kTRUE, 0, 1, ok);
  const Double_t tail_amp = shape->GetParameter(3);
  const Double_t tail_decay = shape->GetParameter(4);
  std::cout << "GAIN558 " << day.label << ": summed line at " << std::fixed
            << std::setprecision(3) << shape->GetParameter(0) << " keV, sigma "
            << shape->GetParameter(1) << ", tail " << tail_amp << " x "
            << tail_decay << " sigma" << (ok ? "" : " FIT FAILED") << std::endl;
  delete shape;
  PixelMap m;
  Int_t n_measured = 0;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (table[p] <= 0 || pixels[p]->GetEntries() < MIN_EVENTS_558)
      continue;
    Bool_t pixel_ok = kFALSE;
    TF1 *f = Fit558(pixels[p], Form("f558_%s_%d", day.tag.Data(), p), kFALSE,
                    tail_amp, tail_decay, pixel_ok);
    if (pixel_ok) {
      Double_t mu = f->GetParameter(0), mu_err = f->GetParError(0);
      Double_t gain = E_CD558 / mu;
      m.value[p] = 1000.0 * E_GE * (gain / table[p] - 1.0);
      m.err[p] = 1000.0 * E_GE * (gain / table[p]) * mu_err / mu;
      n_measured++;
    }
    delete f;
  }
  std::cout << "GAIN558 " << day.label << ": " << n_measured
            << " pixels measured" << std::endl;
  delete sum;
  for (Int_t p = 0; p < N_PIXELS; p++)
    delete pixels[p];
  return m;
}

// Each crystal's weighted mean of a map (crystal 0: x < 0, y < 0; 1: x > 0,
// y < 0; 2: x < 0, y > 0; 3: x > 0, y > 0), and the pixels' rms about their
// crystal's mean.
void PrintCrystalMeans(const PixelMap &m, const TString &label) {
  Double_t sw[4] = {0, 0, 0, 0}, swx[4] = {0, 0, 0, 0};
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (m.err[p] <= 0)
      continue;
    Int_t ix = p / N_SIDE, iy = p % N_SIDE;
    Int_t crystal = (ix >= N_SIDE / 2 ? 1 : 0) + (iy >= N_SIDE / 2 ? 2 : 0);
    Double_t w = 1.0 / (m.err[p] * m.err[p]);
    sw[crystal] += w;
    swx[crystal] += w * m.value[p];
  }
  Double_t mean[4];
  for (Int_t c = 0; c < 4; c++)
    mean[c] = sw[c] > 0 ? swx[c] / sw[c] : 0;
  Double_t s2 = 0, e2 = 0;
  Int_t n = 0;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (m.err[p] <= 0)
      continue;
    Int_t ix = p / N_SIDE, iy = p % N_SIDE;
    Int_t crystal = (ix >= N_SIDE / 2 ? 1 : 0) + (iy >= N_SIDE / 2 ? 2 : 0);
    s2 += TMath::Power(m.value[p] - mean[crystal], 2);
    e2 += m.err[p] * m.err[p];
    n++;
  }
  std::cout << "CRYSTALS " << label << ": means " << std::showpos << std::fixed
            << std::setprecision(1) << mean[0] << " / " << mean[1] << " / "
            << mean[2] << " / " << mean[3] << std::noshowpos << " eV (+-"
            << (sw[0] > 0 ? 1.0 / TMath::Sqrt(sw[0]) : 0)
            << "); pixels about their crystal's mean rms "
            << (n ? TMath::Sqrt(s2 / n) : 0) << " eV against errors' "
            << (n ? TMath::Sqrt(e2 / n) : 0) << " eV" << std::endl;
}

// ------------------------------------------------------------------ drawing

std::vector<Double_t> PixelEdges(const std::vector<Float_t> &centres) {
  std::vector<Double_t> edges;
  edges.push_back(centres[0] - 0.95);
  for (size_t i = 1; i < centres.size(); i++)
    edges.push_back(0.5 * (centres[i - 1] + centres[i]));
  edges.push_back(centres.back() + 0.95);
  return edges;
}

void DrawMap(const PixelMap &m, Double_t shift, const TString &name,
             const TString &z_title, Double_t z_range) {
  std::vector<Double_t> x_edges = PixelEdges(Constants::PIXEL_CENTERS_X_MM);
  std::vector<Double_t> y_edges = PixelEdges(Constants::PIXEL_CENTERS_Y_MM);
  TH2D *h = new TH2D("map_" + name, "", N_SIDE, x_edges.data(), N_SIDE,
                     y_edges.data());
  h->SetDirectory(nullptr);
  for (Int_t ix = 0; ix < N_SIDE; ix++)
    for (Int_t iy = 0; iy < N_SIDE; iy++) {
      Int_t p = ix * N_SIDE + iy;
      // Clamped into the colour range: a cell below the minimum would be
      // left unpainted, like an unmeasured pixel.
      if (m.err[p] > 0)
        h->SetBinContent(
            ix + 1, iy + 1,
            TMath::Max(-0.999 * z_range,
                       TMath::Min(0.999 * z_range, m.value[p] - shift)));
    }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, canvas,
                                             ";x [mm];y [mm];" + z_title);
  // Offsets of either sign: a linear z axis and a diverging palette.
  canvas->SetLogz(kFALSE);
  gStyle->SetPalette(kLightTemperature);
  canvas->SetFrameFillColor(kGray);
  h->SetMinimum(-z_range);
  h->SetMaximum(z_range);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
  delete canvas;
  delete h;
}

// Two days' maps compared: mean, rms, the difference against its errors after
// each map's mean is removed, and the correlation.
void Compare(const TString &what, const TString &unit, const PixelMap &a,
             const PixelMap &b, const TString &label_a, const TString &label_b,
             const TString &name, Double_t z_range) {
  const Double_t mean_a = WeightedMeanOf(a), mean_b = WeightedMeanOf(b);
  Double_t chi2 = 0, sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0, sdd = 0,
           sexp = 0;
  Int_t n = 0;
  PixelMap difference;
  TGraphErrors *scatter = new TGraphErrors();
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (a.err[p] <= 0 || b.err[p] <= 0)
      continue;
    Double_t va = a.value[p] - mean_a, vb = b.value[p] - mean_b;
    Double_t variance = a.err[p] * a.err[p] + b.err[p] * b.err[p];
    chi2 += (va - vb) * (va - vb) / variance;
    sa += va;
    sb += vb;
    saa += va * va;
    sbb += vb * vb;
    sab += va * vb;
    sdd += (va - vb) * (va - vb);
    sexp += variance;
    difference.value[p] = va - vb;
    difference.err[p] = TMath::Sqrt(variance);
    scatter->SetPoint(n, va, vb);
    scatter->SetPointError(n, a.err[p], b.err[p]);
    n++;
  }
  if (n < 2)
    return;
  Double_t cov = sab / n - (sa / n) * (sb / n);
  Double_t var_a = saa / n - (sa / n) * (sa / n);
  Double_t var_b = sbb / n - (sb / n) * (sb / n);
  std::cout << "COMPARE " << what << " " << label_a << " vs " << label_b << ": "
            << n << " pixels; means " << std::showpos << std::fixed
            << std::setprecision(1) << mean_a << " / " << mean_b
            << std::noshowpos << " " << unit << "; rms " << TMath::Sqrt(var_a)
            << " / " << TMath::Sqrt(var_b) << " " << unit << "; difference rms "
            << TMath::Sqrt(sdd / n) << " against its errors' "
            << TMath::Sqrt(sexp / n) << " " << unit << ", chi2/pixel "
            << std::setprecision(2) << chi2 / n << "; correlation "
            << cov / TMath::Sqrt(var_a * var_b) << std::endl;
  DrawMap(difference, 0, name + "_difference",
          Form("%s - %s [%s]", label_a.Data(), label_b.Data(), unit.Data()),
          z_range);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(scatter, kBlue + 1,
                                Form(";%s, %s [%s];%s, %s [%s]", what.Data(),
                                     label_a.Data(), unit.Data(), what.Data(),
                                     label_b.Data(), unit.Data()));
  scatter->SetMarkerStyle(20);
  scatter->SetMarkerSize(0.5);
  scatter->GetXaxis()->SetLimits(-2 * z_range, 2 * z_range);
  scatter->GetYaxis()->SetRangeUser(-2 * z_range, 2 * z_range);
  scatter->Draw("AP");
  TLine *diagonal =
      new TLine(-2 * z_range, -2 * z_range, 2 * z_range, 2 * z_range);
  diagonal->SetLineStyle(2);
  diagonal->Draw();
  PlottingUtils::SaveFigure(canvas, name + "_scatter", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

} // namespace CztPixelGain

void PixelGainStability() {
  using namespace CztPixelGain;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/czt/studies/pixel_gain",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  std::vector<Day> days = Days();
  std::vector<Double_t> table = ProductionGains(project_root);

  // Pb K-alpha on the days behind a shield, the capture line on every day.
  const Int_t N_PB_DAYS = 2;
  std::vector<PixelMap> pb, capture, gain;
  for (size_t d = 0; d < days.size(); d++) {
    if ((Int_t)d < N_PB_DAYS) {
      pb.push_back(LineDayMap(days[d], LINE_PB_KA));
      DrawMap(pb[d], WeightedMeanOf(pb[d]), "pb_offset_" + days[d].tag,
              "Pb K#alpha_{1} offset [eV], " + days[d].label, 150);
      PrintLineDepth(days[d], LINE_PB_KA);
      PrintCrystalMeans(pb[d], "Pb Ka1 " + days[d].label);
    }
    capture.push_back(LineDayMap(days[d], LINE_CAPTURE_96));
    DrawMap(capture[d], WeightedMeanOf(capture[d]),
            "cap96_offset_" + days[d].tag,
            "95.9 keV offset [eV], " + days[d].label, 150);
    PrintLineDepth(days[d], LINE_CAPTURE_96);
    PrintCrystalMeans(capture[d], "95.902 keV " + days[d].label);
    gain.push_back(Gain558DayMap(days[d], table));
    PrintCrystalMeans(gain[d], "558 keV gain " + days[d].label);
    DrawMap(gain[d], WeightedMeanOf(gain[d]), "gain558_" + days[d].tag,
            "558 keV gain offset [eV], " + days[d].label, 150);
    // The measures on the same day: do they see one pattern?
    if ((Int_t)d < N_PB_DAYS)
      Compare("Pb Ka1 vs 95.902", "eV", pb[d], capture[d],
              "Pb " + days[d].label, "96 " + days[d].label,
              "pb_vs_cap96_" + days[d].tag, 150);
    Compare("95.902 vs 558", "eV", capture[d], gain[d], "96 " + days[d].label,
            "558 " + days[d].label, "cap96_vs_558_" + days[d].tag, 150);
  }
  Compare("Pb Ka1", "eV", pb[0], pb[1], days[0].label, days[1].label,
          "pb_offset_" + days[0].tag + "_" + days[1].tag, 150);
  const Int_t N_PAIRS = 4;
  const Int_t PAIRS[N_PAIRS][2] = {{0, 1}, {1, 2}, {2, 3}, {0, 2}};
  for (Int_t k = 0; k < N_PAIRS; k++) {
    Int_t i = PAIRS[k][0], j = PAIRS[k][1];
    Compare("95.902 keV", "eV", capture[i], capture[j], days[i].label,
            days[j].label, "cap96_offset_" + days[i].tag + "_" + days[j].tag,
            150);
    Compare("558 keV gain", "eV", gain[i], gain[j], days[i].label,
            days[j].label, "gain558_" + days[i].tag + "_" + days[j].tag, 150);
  }
}
