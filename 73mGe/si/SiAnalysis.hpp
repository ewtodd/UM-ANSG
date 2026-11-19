#ifndef SI_ANALYSIS_HPP
#define SI_ANALYSIS_HPP

// The PSI Si analysis's shared pieces (Amptek Si on DT5780 channel 0): the
// raw and processed files, the waveform processing (Analysis-Utilities'
// PulseHeightProcessor), the runs, spectra from the processed trees in
// either of two energies, and a binned-likelihood fitter of Gaussian lines
// on a polynomial background with a linear calibration.
//
// The two energies are the trapezoid of the processed waveform and the
// digitiser's own (on-board DPP) energy of the same record. Spectra are in
// channels of the digitiser's scale, the trapezoid divided by its gain to it
// (TRAPEZOID_PER_CHANNEL), so one set of line windows serves both.

#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"
#include "WaveformProcessingUtils.hpp"

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Si {

// ---------------------------------------------------------------- files

const char *const RAW_DIR = "/labdata/PSI/Si/";
const char *const FILE_PREFIX = "Data_CH0@DT5780P_14068_";
const char *const FEATURES_DIR = "/labdata/ANSG/PSI/processed/Si/features";
const Double_t SAMPLE_NS = 10.0;

// A run is the raw file's name without the prefix and ".root", e.g.
// "PSI_Beam_Ge_20251011_1715_2".
inline TString RawFile(const TString &run) {
  return TString(RAW_DIR) + FILE_PREFIX + run + ".root";
}

inline TString FeaturesFile(const TString &run) {
  return TString(FEATURES_DIR) + "/" + FILE_PREFIX + run + ".root";
}

// --------------------------------------------------------- processing

// Records are 2000 samples long up to 10/11 and 4500 from 10/12 on. The
// trigger (15% of the pulse's peak) sits at sample 675 or 1925; above 50 keV
// a few pulses in a thousand trigger earlier, and the records that trigger
// far earlier (about a tenth, on an earlier pulse's tail or noise) are low
// in energy and rejected by the crop.
const Int_t SHORT_RECORD = 2000, LONG_RECORD = 4500;
inline Int_t NominalTrigger(Int_t record_length) {
  return record_length > SHORT_RECORD ? 1925 : 675;
}

// The baseline is the record's first samples, ending this far before the
// nominal trigger, so the baseline under the pulse is measured as close to
// it as the early-triggering pulses allow.
const Int_t BASELINE_GAP = 75;
inline Int_t BaselineSamples(Int_t record_length) {
  return NominalTrigger(record_length) - BASELINE_GAP;
}

// The record length of a raw file, from its first record.
inline Int_t RecordLength(const TString &run) {
  TFile *file = TFile::Open(RawFile(run));
  if (!file || file->IsZombie())
    return -1;
  TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
  Int_t length = -1;
  if (tree && tree->GetEntries() > 0) {
    TArrayS *samples = nullptr;
    tree->SetBranchAddress("Samples", &samples);
    tree->GetEntry(0);
    length = samples->GetSize();
    tree->ResetBranchAddresses();
    delete samples;
  }
  file->Close();
  delete file;
  return length;
}

inline FileProcessingConfig Config(Int_t record_length) {
  FileProcessingConfig config;
  config.polarity = 1;
  config.trigger_threshold = 0.15;
  config.num_samples_baseline = BaselineSamples(record_length);
  config.pre_samples = 500;
  config.post_samples = 1000;
  config.sample_waveforms_to_save = 0;
  config.verbose = kTRUE;
  config.store_waveforms = kFALSE;
  config.input_format = InputFormat::kCOMPASS;
  return config;
}

// Preamplifier decay 324 samples (the source runs' average pulse; the beam
// runs' steps are flat with it). The shaping is SiShapingScan.cpp's choice,
// one for every run so that the calibrations never mix shapings: rise 4 us
// is best in the 2000-sample records and close to it in the 4500-sample
// ones, and the 1 us flat top keeps the reading point well inside the flat
// top of the slower of the two rise bands (180-200 and 310-330 ns). The
// height is read at 90% of the flat top.
inline PulseHeightSettings Filter() {
  PulseHeightSettings filter;
  filter.decay_time_samples = 324.0;
  filter.rise_samples = 400;
  filter.flat_top_samples = 100;
  filter.peaking_fraction = 0.9;
  filter.peak_average_samples = 1;
  filter.rise_plateau_start = 150;
  filter.rise_plateau_samples = 100;
  filter.rise_search_back = 50;
  return filter;
}

// ------------------------------------------------------------- energies

// Records whose raw baseline sits at the bottom of the ADC range have their
// noise clipped (baseline RMS ~2.7 ADC, against ~33) and carry a fake step
// near 11.4 keV; they are dropped below this baseline RMS.
const Float_t FLOOR_BASELINE_RMS = 8.0;
// Digitiser channels per keV, a first guess (Am 59.54 keV near channel 1428).
const Double_t DPP_PER_KEV = 23.98;
// Trapezoid ADC per digitiser channel in the source runs (SiWaveformLook).
const Double_t TRAPEZOID_PER_CHANNEL = 4.163;

// Optional event selection: a time window in seconds from the run's first
// record, isolation (no other accepted record within min_gap_us), and a
// window on the 10-90% rise time in ns.
struct Selection {
  Double_t t_min = -1, t_max = -1;
  Double_t min_gap_us = -1;
  Double_t rise_min_ns = -1, rise_max_ns = -1;
};

const Int_t N_CHANNELS = 4096;

inline TH1D *NewSpectrum(const TString &name) {
  TH1D *h = new TH1D(name, ";Energy [channel];Counts / channel", N_CHANNELS,
                     -0.5, N_CHANNELS - 0.5);
  h->SetDirectory(nullptr);
  return h;
}

// One run's trapezoid spectrum added to h, in digitiser channels.
inline void FillSpectrum(TH1D *h, const TString &run,
                         const Selection &selection) {
  TFile file(FeaturesFile(run));
  TTree *tree = static_cast<TTree *>(file.Get("features"));
  if (!tree) {
    std::cout << "SI no features tree for " << run << std::endl;
    return;
  }
  ULong64_t timestamp = 0;
  Float_t trapezoid_height = 0, baseline_rms = 0, rise_time = 0;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("timestamp", 1);
  tree->SetBranchStatus("trapezoid_height", 1);
  tree->SetBranchStatus("baseline_rms", 1);
  tree->SetBranchStatus("rise_time", 1);
  tree->SetBranchAddress("timestamp", &timestamp);
  tree->SetBranchAddress("trapezoid_height", &trapezoid_height);
  tree->SetBranchAddress("baseline_rms", &baseline_rms);
  tree->SetBranchAddress("rise_time", &rise_time);
  Long64_t n = tree->GetEntries();
  Bool_t need_gaps = selection.min_gap_us > 0;
  std::vector<ULong64_t> times;
  if (need_gaps) {
    times.resize(n);
    for (Long64_t i = 0; i < n; i++) {
      tree->GetEntry(i);
      times[i] = timestamp;
    }
  }
  ULong64_t first_time = 0;
  for (Long64_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    if (i == 0)
      first_time = timestamp;
    if (baseline_rms < FLOOR_BASELINE_RMS)
      continue;
    Double_t t_s = (timestamp - first_time) * 1e-12;
    if (selection.t_min >= 0 && t_s < selection.t_min)
      continue;
    if (selection.t_max >= 0 && t_s >= selection.t_max)
      continue;
    Double_t rise_ns = SAMPLE_NS * rise_time;
    if (selection.rise_min_ns >= 0 && rise_ns < selection.rise_min_ns)
      continue;
    if (selection.rise_max_ns >= 0 && rise_ns >= selection.rise_max_ns)
      continue;
    if (need_gaps) {
      Double_t gap_before = i > 0 ? (times[i] - times[i - 1]) * 1e-6 : 1e30;
      Double_t gap_after = i + 1 < n ? (times[i + 1] - times[i]) * 1e-6 : 1e30;
      if (gap_before < selection.min_gap_us || gap_after < selection.min_gap_us)
        continue;
    }
    h->Fill(trapezoid_height / TRAPEZOID_PER_CHANNEL);
  }
}

inline TH1D *SumSpectrum(const std::vector<TString> &runs, const TString &name,
                         const Selection &selection = Selection()) {
  TH1D *h = NewSpectrum(name);
  for (size_t i = 0; i < runs.size(); i++)
    FillSpectrum(h, runs[i], selection);
  return h;
}

// Seconds from a run's first to its last record, and its record count.
inline void RunSpan(const TString &run, Double_t &seconds, Long64_t &entries) {
  seconds = 0;
  entries = 0;
  TFile file(FeaturesFile(run));
  TTree *tree = static_cast<TTree *>(file.Get("features"));
  if (!tree || tree->GetEntries() == 0)
    return;
  ULong64_t timestamp = 0;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("timestamp", 1);
  tree->SetBranchAddress("timestamp", &timestamp);
  entries = tree->GetEntries();
  tree->GetEntry(0);
  ULong64_t first = timestamp;
  tree->GetEntry(entries - 1);
  seconds = (timestamp - first) * 1e-12;
}

// --------------------------------------------------------------- fitter

// Gaussian lines on a polynomial background, in channels. Parameters: [0..2]
// background c0 + c1 (x - x_ref) + c2 (x - x_ref)^2 per channel, then
// n_sigma widths, then per line its area (counts) and mean.
struct LineModel {
  Int_t n_lines = 1;
  Int_t n_sigma = 1;
  std::vector<Int_t> sigma_slot;
  Double_t x_ref = 0;
  Int_t NPar() const { return 3 + n_sigma + 2 * n_lines; }
  Int_t SigmaPar(Int_t s) const { return 3 + s; }
  Int_t AreaPar(Int_t k) const { return 3 + n_sigma + 2 * k; }
  Int_t MeanPar(Int_t k) const { return 4 + n_sigma + 2 * k; }
  Double_t Background(Double_t x, const Double_t *p) const {
    Double_t d = x - x_ref;
    return p[0] + p[1] * d + p[2] * d * d;
  }
  Double_t Line(Double_t x, const Double_t *p, Int_t k) const {
    Double_t sigma = p[SigmaPar(sigma_slot[k])];
    Double_t z = (x - p[MeanPar(k)]) / sigma;
    return p[AreaPar(k)] * TMath::Exp(-0.5 * z * z) /
           (TMath::Sqrt(2.0 * TMath::Pi()) * sigma);
  }
  Double_t operator()(Double_t *x, Double_t *p) {
    Double_t v = Background(x[0], p);
    for (Int_t k = 0; k < n_lines; k++)
      v += Line(x[0], p, k);
    return v;
  }
};

struct LineSpec {
  Double_t mean_guess = 0;
  Int_t sigma_slot = 0;
  // Channels either side of the guess the mean may move.
  Double_t mean_range = 15;
};

struct FitOutput {
  Bool_t converged = kFALSE;
  std::vector<Double_t> mean, mean_err, area, area_err;
  std::vector<Double_t> sigma, sigma_err;
  Double_t c0 = 0, c1 = 0, c2 = 0;
  // Pearson chi2 against the fitted model.
  Double_t chi2 = 0;
  Int_t ndf = 0;
  Double_t lo = 0, hi = 0;
  Int_t status = -1;
};

inline Int_t &FitCounter() {
  static Int_t counter = 0;
  return counter;
}

inline Int_t &UnconvergedFitCounter() {
  static Int_t counter = 0;
  return counter;
}

// Fit option "M" runs TMinuit's MIGRAD then IMPROVE, and adds 1000 times
// IMPROVE's return code to the status; IMPROVE returns 4 when it finds no
// lower minimum than MIGRAD's, so a good fit reads 4000.
inline Bool_t FitConverged(const TFitResultPtr &result) {
  if (!result.Get())
    return kFALSE;
  Int_t migrad_and_hesse_status = result->Status() % 1000;
  return migrad_and_hesse_status == 0 && result->CovMatrixStatus() == 3;
}

// Binned Poisson likelihood fit of the lines on a background of order
// background_order (0 flat, 1 linear, 2 quadratic) over [lo, hi] channels.
inline FitOutput FitLines(TH1D *h, Double_t lo, Double_t hi,
                          const std::vector<LineSpec> &lines, Int_t n_sigma,
                          Int_t background_order, Double_t sigma_guess = 6.0,
                          Bool_t fix_sigma = kFALSE) {
  FitOutput out;
  out.lo = lo;
  out.hi = hi;
  LineModel model;
  model.n_lines = (Int_t)lines.size();
  model.n_sigma = n_sigma;
  model.x_ref = 0.5 * (lo + hi);
  for (size_t k = 0; k < lines.size(); k++)
    model.sigma_slot.push_back(lines[k].sigma_slot);
  TF1 f(Form("si_fit_%d", FitCounter()++), model, lo, hi, model.NPar());
  f.SetNpx(1000);
  Int_t bin_lo = h->FindBin(lo), bin_hi = h->FindBin(hi);
  Double_t edge = 0;
  Int_t n_edge = 0;
  for (Int_t b = bin_lo; b < bin_lo + 5; b++) {
    edge += h->GetBinContent(b);
    n_edge++;
  }
  for (Int_t b = bin_hi - 4; b <= bin_hi; b++) {
    edge += h->GetBinContent(b);
    n_edge++;
  }
  Double_t level = TMath::Max(edge / n_edge, 0.1);
  f.SetParameter(0, level);
  f.SetParameter(1, 0);
  f.SetParameter(2, 0);
  if (background_order < 2)
    f.FixParameter(2, 0);
  if (background_order < 1)
    f.FixParameter(1, 0);
  for (Int_t s = 0; s < n_sigma; s++) {
    f.SetParameter(model.SigmaPar(s), sigma_guess);
    f.SetParLimits(model.SigmaPar(s), 1.0, 40.0);
    if (fix_sigma)
      f.FixParameter(model.SigmaPar(s), sigma_guess);
  }
  for (size_t k = 0; k < lines.size(); k++) {
    Int_t b = h->FindBin(lines[k].mean_guess);
    Double_t height = TMath::Max(h->GetBinContent(b) - level, 1.0);
    f.SetParameter(model.AreaPar(k), height * 2.5 * sigma_guess);
    f.SetParLimits(model.AreaPar(k), 0.0, 1e9);
    f.SetParameter(model.MeanPar(k), lines[k].mean_guess);
    f.SetParLimits(model.MeanPar(k), lines[k].mean_guess - lines[k].mean_range,
                   lines[k].mean_guess + lines[k].mean_range);
  }
  h->Fit(&f, "QLRN0");
  TFitResultPtr result = h->Fit(&f, "QLRN0SM");
  out.status = result.Get() ? result->Status() : -1;
  out.converged = FitConverged(result);
  if (!out.converged) {
    UnconvergedFitCounter()++;
    std::cout << "SI FIT NOT CONVERGED: " << h->GetName() << " over " << lo
              << "-" << hi << " ch, " << lines.size() << " lines, first at "
              << lines[0].mean_guess << " ch, background order "
              << background_order << ", status " << out.status
              << ", covariance status "
              << (result.Get() ? result->CovMatrixStatus() : -1) << std::endl;
  }
  for (size_t k = 0; k < lines.size(); k++) {
    out.mean.push_back(f.GetParameter(model.MeanPar(k)));
    out.mean_err.push_back(f.GetParError(model.MeanPar(k)));
    out.area.push_back(f.GetParameter(model.AreaPar(k)));
    out.area_err.push_back(f.GetParError(model.AreaPar(k)));
  }
  for (Int_t s = 0; s < n_sigma; s++) {
    out.sigma.push_back(f.GetParameter(model.SigmaPar(s)));
    out.sigma_err.push_back(f.GetParError(model.SigmaPar(s)));
  }
  out.c0 = f.GetParameter(0);
  out.c1 = f.GetParameter(1);
  out.c2 = f.GetParameter(2);
  Double_t chi2 = 0;
  Int_t n_bins = 0;
  for (Int_t b = bin_lo; b <= bin_hi; b++) {
    Double_t expected = f.Eval(h->GetBinCenter(b));
    if (expected <= 0)
      continue;
    Double_t d = h->GetBinContent(b) - expected;
    chi2 += d * d / expected;
    n_bins++;
  }
  out.chi2 = chi2;
  out.ndf = n_bins - f.GetNumberFreeParameters();
  return out;
}

// ----------------------------------------------------------- calibration

// E = c0 + c1 (ch - pivot): a weighted least-squares line, the channel
// errors turned into energy errors by the slope, iterated.
struct Calibration {
  Double_t c0 = 0, c1 = 1.0 / DPP_PER_KEV, pivot = 1650;
  Double_t v00 = 0, v11 = 0, v01 = 0;
  Double_t chi2 = 0;
  Int_t ndf = 0;
  Double_t Energy(Double_t ch) const { return c0 + c1 * (ch - pivot); }
  Double_t Channel(Double_t e) const { return pivot + (e - c0) / c1; }
  // Calibration error at channel ch, and its offset (c0) and slope (c1)
  // parts.
  Double_t Error(Double_t ch) const {
    Double_t d = ch - pivot;
    return TMath::Sqrt(TMath::Max(v00 + d * d * v11 + 2 * d * v01, 0.0));
  }
  Double_t OffsetError() const { return TMath::Sqrt(v00); }
  Double_t SlopeError(Double_t ch) const {
    return TMath::Abs(ch - pivot) * TMath::Sqrt(v11);
  }
};

// A negative pivot puts it at the points' weighted centroid, where offset
// and slope are uncorrelated.
inline Calibration FitCalibration(const std::vector<Double_t> &ch,
                                  const std::vector<Double_t> &ch_err,
                                  const std::vector<Double_t> &energy,
                                  const std::vector<Double_t> &energy_err,
                                  Double_t pivot) {
  Calibration cal;
  if (pivot < 0) {
    Double_t sw = 0, swx = 0;
    for (size_t i = 0; i < ch.size(); i++) {
      Double_t w = 1.0 / (ch_err[i] * ch_err[i]);
      sw += w;
      swx += w * ch[i];
    }
    pivot = swx / sw;
  }
  cal.pivot = pivot;
  for (Int_t iteration = 0; iteration < 5; iteration++) {
    Double_t s = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < ch.size(); i++) {
      Double_t err = TMath::Sqrt(TMath::Power(cal.c1 * ch_err[i], 2) +
                                 TMath::Power(energy_err[i], 2));
      Double_t w = 1.0 / (err * err);
      Double_t x = ch[i] - pivot;
      s += w;
      sx += w * x;
      sy += w * energy[i];
      sxx += w * x * x;
      sxy += w * x * energy[i];
    }
    Double_t det = s * sxx - sx * sx;
    cal.c0 = (sxx * sy - sx * sxy) / det;
    cal.c1 = (s * sxy - sx * sy) / det;
    cal.v00 = sxx / det;
    cal.v11 = s / det;
    cal.v01 = -sx / det;
  }
  cal.chi2 = 0;
  for (size_t i = 0; i < ch.size(); i++) {
    Double_t err = TMath::Sqrt(TMath::Power(cal.c1 * ch_err[i], 2) +
                               TMath::Power(energy_err[i], 2));
    cal.chi2 += TMath::Power((cal.Energy(ch[i]) - energy[i]) / err, 2);
  }
  cal.ndf = (Int_t)ch.size() - 2;
  return cal;
}

} // namespace Si

#endif
