#ifndef PIXEL_ALIGNMENT_HPP
#define PIXEL_ALIGNMENT_HPP

// The day pixel map: each pixel's energy scale against its run's, measured
// on a line of the day's beam runs. The pixel gains (PixelCalibration.cpp)
// are one table, from the 558.456 keV line of the 01/16 graphite-castle run
// at every depth; near 70 keV the crystals still read up to 300 eV apart and
// the pixels within a crystal about 110 eV rms, a pattern the 558 keV gains do
// not see (PixelGainStability.cpp). The lines, from the filtered events
// (z >= 8 mm), absorbed near the cathode as Ge is:
//   Pb K-alpha   6 keV above Ge; strong behind the Cu or Cd shield (01/13), a
//                few percent above the continuum without it;
//   95.902 keV   113Cd(n,g) in the crystal's cathode layer, in every beam run
//                (01/15).
// A pixel's line position against its run's summed spectrum, averaged over
// the day's beam runs, gives its factor E -> E mu_run / mu_pixel, which leaves
// a run's summed scale where it was. Each pixel's line is fitted with the
// summed spectrum's shape held (with it free, a few hundred events cannot fix
// the tail), its position, yield and flat continuum free.

#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "ReferenceLines.hpp"

#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

namespace PixelAlignment {

const Int_t N_SIDE = 22;
const Int_t N_PIXELS = N_SIDE * N_SIDE;

enum { LINE_PB_KA, LINE_CAPTURE_96, N_LINES };
const char *const LINE_NAMES[N_LINES] = {"Pb Ka1", "95.902 keV"};
// The fit windows, in the filtered trees' energies (0.98-1.02 keV per unit).
const Double_t WINDOW_LO[N_LINES] = {70.0, 91.5};
const Double_t WINDOW_HI[N_LINES] = {77.5, 100.5};
// Offsets are quoted in eV at Ge's energy, the one that matters.
const Double_t E_QUOTE = 68.75;
// Events in the window a pixel needs for its own fit.
const Int_t MIN_EVENTS = 150;
// A pixel's factor is applied only where it is measured to this (eV at Ge's
// energy): a worse measurement adds more spread than the pattern it removes.
const Double_t MAX_OFFSET_ERR_EV = 100.0;

// The pixel of an event on a pixel centre (every filtered event), or -1.
inline Int_t PixelIndexOf(Float_t x_mm, Float_t y_mm) {
  Int_t ix = -1, iy = -1;
  for (Int_t i = 0; i < N_SIDE; i++) {
    if (TMath::Abs(x_mm - Constants::PIXEL_CENTERS_X_MM[i]) <=
        Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      ix = i;
    if (TMath::Abs(y_mm - Constants::PIXEL_CENTERS_Y_MM[i]) <=
        Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      iy = i;
  }
  if (ix < 0 || iy < 0)
    return -1;
  return ix * N_SIDE + iy;
}

// A line whose Gaussian core holds `counts` events, with a low exponential
// tail (amplitude a ratio of the core's height, decay in units of sigma), per
// bin of width bin.
inline Double_t Line(Double_t x, Double_t mu, Double_t sigma, Double_t counts,
                     Double_t tail_amp, Double_t tail_decay, Double_t bin) {
  Double_t height = counts * bin / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma);
  Double_t par[12] = {mu, sigma, height, 0,   tail_amp, tail_decay,
                      0,  0,     0,      1.0, 0,        0};
  return FittingFunctions::PeakFunction(&x, par);
}

// [0] main line's position [1] sigma [2] its core counts [3] the partner
// line's yield over it (Pb K-alpha2; 0 for the 95.902 keV line) [4] main -
// partner spacing [5] tail amplitude [6] tail decay [7] continuum per bin
// [8] bin width
inline Double_t LineModel(Double_t *x, Double_t *p) {
  return p[7] + Line(x[0], p[0], p[1], p[2], p[5], p[6], p[8]) +
         Line(x[0], p[0] - p[4], p[1], p[3] * p[2], p[5], p[6], p[8]);
}

struct RunAlignment {
  TString run;
  Int_t line = LINE_PB_KA;
  Bool_t valid = kFALSE;
  // The summed spectrum's line position, and the shape the pixels' fits hold.
  Double_t mu = 0, mu_err = 0, sigma = 0, ka2_ratio = 0, spacing = 0;
  Double_t tail_amp = 0, tail_decay = 1.0;
  std::vector<Double_t> pixel_mu, pixel_mu_err;
  std::vector<Int_t> pixel_events;
  std::vector<Bool_t> measured;

  // The pixel's line against the run's, in eV at Ge's energy.
  Double_t OffsetEv(Int_t pixel) const {
    return 1000.0 * E_QUOTE * (pixel_mu[pixel] / mu - 1.0);
  }
  Double_t OffsetErrEv(Int_t pixel) const {
    return 1000.0 * E_QUOTE * pixel_mu_err[pixel] / mu;
  }
  Double_t Factor(Int_t pixel) const {
    if (!valid || pixel < 0 || !measured[pixel] ||
        OffsetErrEv(pixel) > MAX_OFFSET_ERR_EV)
      return 1.0;
    return mu / pixel_mu[pixel];
  }
};

inline TH1D *NewWindowHist(const TString &name, Int_t line) {
  const Double_t bin = Constants::BIN_WIDTH_KEV;
  const Int_t nbins =
      (Int_t)std::lround((WINDOW_HI[line] - WINDOW_LO[line]) / bin);
  TH1D *h =
      new TH1D(name, "", nbins, WINDOW_LO[line], WINDOW_LO[line] + nbins * bin);
  h->SetDirectory(nullptr);
  return h;
}

// The run's filtered events in the line's window, filled into the summed and
// the per-pixel histograms.
inline Bool_t FillWindow(const TString &run, Int_t line, TH1D *sum,
                         std::vector<TH1D *> &pixels) {
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << run << ".root" << std::endl;
    delete file;
    return kFALSE;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t energy = 0, x = 0, y = 0;
    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
      tree->GetEntry(i);
      if (energy < WINDOW_LO[line] || energy >= WINDOW_HI[line])
        continue;
      Int_t pixel = PixelIndexOf(x, y);
      if (pixel < 0)
        continue;
      sum->Fill(energy);
      pixels[pixel]->Fill(energy);
    }
    tree->ResetBranchAddresses();
  }
  file->Close();
  delete file;
  return kTRUE;
}

// A fit's status is good when Minuit converged; with IMPROVE ("M") a good
// fit reports 4000.
inline Bool_t Converged(const TFitResultPtr &r) {
  return r.Get() != nullptr && (r->Status() == 0 || r->Status() == 4000);
}

inline TF1 *NewLineFit(const TString &name, Int_t line) {
  TF1 *f = new TF1(name, LineModel, WINDOW_LO[line], WINDOW_HI[line], 9);
  f->FixParameter(8, Constants::BIN_WIDTH_KEV);
  return f;
}

// The summed spectrum's line, every shape parameter free. The main line is
// the window's highest peak; where it sits in the uncalibrated energies
// changes from day to day.
inline Bool_t FitSum(TH1D *h, RunAlignment &a) {
  TF1 *f = NewLineFit(Form("line_sum_%s_%d", a.run.Data(), a.line), a.line);
  Double_t top = h->GetBinCenter(h->GetMaximumBin());
  Double_t height = h->GetMaximum();
  const Double_t sigma0 = 0.6;
  f->SetParameter(0, top);
  f->SetParLimits(0, top - 1.0, top + 1.0);
  f->SetParameter(1, sigma0);
  f->SetParLimits(1, 0.2, 1.5);
  f->SetParameter(2, height * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 /
                         Constants::BIN_WIDTH_KEV);
  f->SetParLimits(2, 0, 1e9);
  if (a.line == LINE_PB_KA) {
    f->SetParameter(3, 0.6);
    f->SetParLimits(3, 0.3, 1.0);
    f->SetParameter(4,
                    0.99 * (ReferenceLines::PB_KA1 - ReferenceLines::PB_KA2));
    f->SetParLimits(4, 2.0, 2.3);
  } else {
    f->FixParameter(3, 0);
    f->FixParameter(4, 2.0);
  }
  f->SetParameter(5, 0.1);
  f->SetParLimits(5, 0, 2.0);
  f->SetParameter(6, 3.0);
  f->SetParLimits(6, 0.5, 20.0);
  f->SetParameter(7, TMath::Max(h->GetBinContent(h->GetNbinsX()), 1.0));
  f->SetParLimits(7, 0, 1e9);
  h->Fit(f, "QLRN0");
  TFitResultPtr r = h->Fit(f, "QLRN0MS");
  a.valid = Converged(r);
  a.mu = f->GetParameter(0);
  a.mu_err = f->GetParError(0);
  a.sigma = f->GetParameter(1);
  a.ka2_ratio = f->GetParameter(3);
  a.spacing = f->GetParameter(4);
  a.tail_amp = f->GetParameter(5);
  a.tail_decay = f->GetParameter(6);
  delete f;
  return a.valid;
}

// A pixel's line with the summed shape held.
inline Bool_t FitPixel(TH1D *h, const RunAlignment &a, Double_t &mu,
                       Double_t &mu_err) {
  TF1 *f = NewLineFit(TString("line_pixel_") + h->GetName(), a.line);
  f->SetParameter(0, a.mu);
  f->SetParLimits(0, a.mu - 1.0, a.mu + 1.0);
  f->FixParameter(1, a.sigma);
  f->SetParameter(2, 0.5 * h->Integral());
  f->SetParLimits(2, 0, 1e9);
  f->FixParameter(3, a.ka2_ratio);
  f->FixParameter(4, a.spacing);
  f->FixParameter(5, a.tail_amp);
  f->FixParameter(6, a.tail_decay);
  f->SetParameter(7, TMath::Max(h->GetBinContent(h->GetNbinsX()), 0.5));
  f->SetParLimits(7, 0, 1e9);
  TFitResultPtr r = h->Fit(f, "QLRN0S");
  mu = f->GetParameter(0);
  mu_err = f->GetParError(0);
  Bool_t ok = (Converged(r) && mu_err > 0 && TMath::Abs(mu - a.mu) < 0.9);
  delete f;
  return ok;
}

inline RunAlignment Measure(const TString &run, Int_t line) {
  RunAlignment a;
  a.run = run;
  a.line = line;
  a.pixel_mu.assign(N_PIXELS, 0);
  a.pixel_mu_err.assign(N_PIXELS, 0);
  a.pixel_events.assign(N_PIXELS, 0);
  a.measured.assign(N_PIXELS, kFALSE);
  TH1D *sum = NewWindowHist(Form("window_sum_%s_%d", run.Data(), line), line);
  std::vector<TH1D *> pixels;
  for (Int_t p = 0; p < N_PIXELS; p++)
    pixels.push_back(
        NewWindowHist(Form("window_%s_%d_%d", run.Data(), line, p), line));
  if (FillWindow(run, line, sum, pixels) && FitSum(sum, a)) {
    for (Int_t p = 0; p < N_PIXELS; p++) {
      a.pixel_events[p] = (Int_t)pixels[p]->GetEntries();
      if (a.pixel_events[p] < MIN_EVENTS)
        continue;
      a.measured[p] = FitPixel(pixels[p], a, a.pixel_mu[p], a.pixel_mu_err[p]);
    }
  }
  Int_t n_measured = 0, n_applied = 0;
  Double_t s = 0, s2 = 0;
  std::vector<Double_t> errors;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (!a.measured[p])
      continue;
    n_measured++;
    if (a.OffsetErrEv(p) <= MAX_OFFSET_ERR_EV)
      n_applied++;
    s += a.OffsetEv(p);
    s2 += a.OffsetEv(p) * a.OffsetEv(p);
    errors.push_back(a.OffsetErrEv(p));
  }
  std::sort(errors.begin(), errors.end());
  Double_t mean = n_measured ? s / n_measured : 0;
  std::cout << "ALIGN " << run << ", " << LINE_NAMES[line]
            << ": summed line at " << std::fixed << std::setprecision(4) << a.mu
            << " (sigma " << std::setprecision(3) << a.sigma << ", tail "
            << a.tail_amp << " x " << a.tail_decay << " sigma)"
            << (a.valid ? "" : " FIT FAILED") << "; " << n_measured
            << " pixels measured, " << n_applied << " to " << MAX_OFFSET_ERR_EV
            << " eV; offsets rms " << std::setprecision(1)
            << (n_measured ? TMath::Sqrt(s2 / n_measured - mean * mean) : 0)
            << " eV, median error "
            << (errors.empty() ? 0 : errors[errors.size() / 2])
            << " eV (at 68.75 keV)" << std::endl;
  delete sum;
  for (Int_t p = 0; p < N_PIXELS; p++)
    delete pixels[p];
  return a;
}

// Each run and line measured once per process.
inline const RunAlignment &AlignmentOf(const TString &run, Int_t line) {
  static std::map<TString, RunAlignment> measured_runs;
  TString key = Form("%s_%d", run.Data(), line);
  std::map<TString, RunAlignment>::iterator found = measured_runs.find(key);
  if (found != measured_runs.end())
    return found->second;
  measured_runs[key] = Measure(run, line);
  return measured_runs[key];
}

// A day's map: each pixel's line against its run's, averaged over the day's
// runs with their errors. Its mu is 1 and its pixel_mu the pixels' ratios, so
// Factor and OffsetEv read as a run's do.
inline RunAlignment CombineRuns(const TString &label,
                                const std::vector<TString> &runs, Int_t line) {
  RunAlignment day;
  day.run = label;
  day.line = line;
  day.mu = 1.0;
  day.pixel_mu.assign(N_PIXELS, 1.0);
  day.pixel_mu_err.assign(N_PIXELS, 0);
  day.pixel_events.assign(N_PIXELS, 0);
  day.measured.assign(N_PIXELS, kFALSE);
  std::vector<Double_t> sw(N_PIXELS, 0), swx(N_PIXELS, 0);
  for (size_t r = 0; r < runs.size(); r++) {
    const RunAlignment &a = AlignmentOf(runs[r], line);
    if (!a.valid)
      continue;
    day.valid = kTRUE;
    for (Int_t p = 0; p < N_PIXELS; p++) {
      if (!a.measured[p])
        continue;
      Double_t ratio = a.pixel_mu[p] / a.mu;
      Double_t ratio_err = a.pixel_mu_err[p] / a.mu;
      sw[p] += 1.0 / (ratio_err * ratio_err);
      swx[p] += ratio / (ratio_err * ratio_err);
      day.pixel_events[p] += a.pixel_events[p];
    }
  }
  Int_t n_applied = 0;
  for (Int_t p = 0; p < N_PIXELS; p++) {
    if (sw[p] <= 0)
      continue;
    day.measured[p] = kTRUE;
    day.pixel_mu[p] = swx[p] / sw[p];
    day.pixel_mu_err[p] = 1.0 / TMath::Sqrt(sw[p]);
    if (day.OffsetErrEv(p) <= MAX_OFFSET_ERR_EV)
      n_applied++;
  }
  std::cout << "ALIGN day " << label << ", " << LINE_NAMES[line] << ", from "
            << runs.size() << " runs: " << n_applied << " pixels measured to "
            << MAX_OFFSET_ERR_EV << " eV" << std::endl;
  return day;
}

inline const RunAlignment &AlignmentOfDay(const TString &label,
                                          const std::vector<TString> &runs,
                                          Int_t line) {
  static std::map<TString, RunAlignment> combined_days;
  TString key = Form("%s_%d", label.Data(), line);
  std::map<TString, RunAlignment>::iterator found = combined_days.find(key);
  if (found != combined_days.end())
    return found->second;
  combined_days[key] = CombineRuns(label, runs, line);
  return combined_days[key];
}

} // namespace PixelAlignment

#endif
