#ifndef GAINTRANSFER_FIT_HPP
#define GAINTRANSFER_FIT_HPP

// The fitting core of the CZT 73Ge analysis.
//
// CalibrateAndFit.cpp (the measurement) and InjectionTest.cpp (its testbed
// closure) both include this, so the testbed exercises the same in-situ
// (Method 1) and direct (Method 2) fits that produce the result.

#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "PixelAlignment.hpp"
#include "PlottingUtils.hpp"
#include "ReferenceLines.hpp"
#include "RooFitPhotopeakPdfs.hpp"
#include "RooFitUtils.hpp"

#include <RooRealVar.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------- constants

const Float_t E_AM241 = ReferenceLines::AM241;
const Float_t E_BA133_79 = ReferenceLines::BA133_79;
const Float_t E_BA133_81 = ReferenceLines::BA133_81;
// ENSDF intensities, 2.65(5) and 32.9(3) per 100 decays.
const Double_t BA133_79_OVER_81 = 2.65 / 32.9;
const Float_t E_PB_KA2 = ReferenceLines::PB_KA2;
const Float_t E_PB_KA1 = ReferenceLines::PB_KA1;

const Float_t E_NEAR_GE = 69;

// Anchor reference-energy uncertainties (keV), the cal fit's y-errors (energy
// is its y axis), so they propagate into the pol1 parameter covariance and
// thus the Ge energy. dE_PB_KA rounds SRD 128's 0.24 and 0.17 eV up to 1 eV;
// dE_AM241 is DDEP's.
//
// dE_LINE is a conservative ~10 eV lump for the Ba-133 lines of the 01/15
// source calibration.
const Float_t dE_PB_KA = 0.001;
const Float_t dE_AM241 = 0.0001;
const Float_t dE_LINE = 0.010;

const Float_t CAL_RANGE_LOW = 65;
const Float_t CAL_RANGE_HIGH = 80;
const Float_t RATESUB_LO = 62.0;
const Float_t RATESUB_HI = 71.5;

// The summed spectra carry a comb with this period in raw energy, the same on
// 01/13 and 01/15 (czt/CombPeriod.cpp measures it, results/comb_period.result).
// Bins one period wide (times the day's slope) hold one tooth each wherever
// they start; 255 eV bins beat against it. The fit plots and their chi2, and
// the 95.902 keV line's fits (averaged over N_RATESUB_PHASES starting points),
// use them.
const Double_t COMB_PERIOD_RAW_KEV = 0.264;
const Int_t N_RATESUB_PHASES = 5;

// TRANSFER=offset: the Cd pairs inherit the Cu 10% calibration shifted by
// their Pb K-alpha offset from it instead of scaled by the Pb ratio. At Ge the
// two differ by ~7% of the Pb shift and no line inside the calibrated range
// resolves them (Pb K-beta lies above it), so the difference is a systematic
// (CombineGeResult's transfer model).
inline Bool_t OffsetTransferRequested() {
  const char *env = getenv("TRANSFER");
  return env && TString(env) == "offset";
}

// The 95.902 keV line's K X-ray escapes enter both methods' fits as fixed
// lines (see "95.902 keV escapes" below). ESCAPES=off leaves them out, for
// CombineGeResult's escape systematic.
inline Bool_t EscapeLinesRequested() {
  const char *env = getenv("ESCAPES");
  return !(env && TString(env) == "off");
}

// ESCAPE_FRACTIONS=up or down: every run's escape fractions moved up or down
// by their measured uncertainty (ESCAPES_INSIDE_REL_ERR), for
// CombineGeResult's escape systematic. +1, -1, or 0 for the measured values.
inline Int_t EscapeFractionShift() {
  const char *env = getenv("ESCAPE_FRACTIONS");
  if (env && TString(env) == "up")
    return 1;
  if (env && TString(env) == "down")
    return -1;
  return 0;
}

// For the dissertation's account of what the pixel gain matching does to the
// results, not for a systematic: PIXEL_MAP=off leaves out the day pixel map
// below, and PIXEL_GAINS=off also undoes the 558 keV gains the filter
// applied (each energy over its pixel's gain), which implies PIXEL_MAP=off.
inline Bool_t PixelGainsRequested() {
  const char *env = getenv("PIXEL_GAINS");
  return !(env && TString(env) == "off");
}

inline Bool_t PixelMapRequested() {
  const char *env = getenv("PIXEL_MAP");
  return PixelGainsRequested() && !(env && TString(env) == "off");
}

// Read once, since every event's energy checks them.
const Bool_t PIXEL_GAINS_ON = PixelGainsRequested();
const Bool_t PIXEL_MAP_ON = PixelMapRequested();

// The 558 keV gain table the filter applied (PixelCalibration.cpp), indexed
// as PixelAlignment::PixelIndexOf, for PIXEL_GAINS=off.
inline const std::vector<Double_t> &PixelGains558() {
  static std::vector<Double_t> gains;
  if (!gains.empty())
    return gains;
  gains.assign(PixelAlignment::N_PIXELS, 0.0);
  TFile *file = IO::OpenForReading("pixel_calibration.root");
  TH2D *table = file ? static_cast<TH2D *>(file->Get("pixel_gain")) : nullptr;
  if (!table) {
    std::cerr << "ERROR: no pixel_gain in root_files/pixel_calibration.root"
              << std::endl;
    return gains;
  }
  for (Int_t ix = 0; ix < PixelAlignment::N_SIDE; ix++)
    for (Int_t iy = 0; iy < PixelAlignment::N_SIDE; iy++)
      gains[ix * PixelAlignment::N_SIDE + iy] =
          table->GetBinContent(ix + 1, iy + 1);
  file->Close();
  delete file;
  return gains;
}

// ---------------------------------------------------------------- pixel map

// Every event is aligned on its day's pixel map (PixelAlignment.hpp) before
// the calibration: each pixel's line position against its run's, averaged
// over the day's beam runs, on Pb K-alpha on 01/13 (behind the Cu or Cd
// shield) and on the 95.902 keV capture line on 01/15 (Pb K-alpha stands a
// few percent above the continuum there). A day's source runs take its map.

inline TString RunDate(const TString &run) {
  if (run.Contains("20260113"))
    return "20260113";
  if (run.Contains("20260115"))
    return "20260115";
  return "";
}

// The beam runs whose line makes a day's map.
inline std::vector<TString> PixelMapRuns(const TString &date) {
  std::vector<TString> runs;
  if (date == "20260113") {
    runs.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);
    runs.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
    runs.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
    runs.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
    runs.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
    runs.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
    runs.push_back(Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113);
    runs.push_back(Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113);
  }
  if (date == "20260115") {
    runs.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115);
    runs.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_20260115);
  }
  return runs;
}

// The map run's pixels take, or null for a run outside the two days.
inline const PixelAlignment::RunAlignment *PixelMapFor(const TString &run) {
  const TString date = RunDate(run);
  if (date == "" || !PIXEL_MAP_ON)
    return nullptr;
  const Int_t line = (date == "20260113") ? PixelAlignment::LINE_PB_KA
                                          : PixelAlignment::LINE_CAPTURE_96;
  return &PixelAlignment::AlignmentOfDay(date, PixelMapRuns(date), line);
}

// An event's energy with its pixel's factor (or, with PIXEL_GAINS=off,
// without its 558 keV gain).
inline Double_t AlignedEnergy(const PixelAlignment::RunAlignment *map,
                              Float_t energy, Float_t x_mm, Float_t y_mm) {
  if (!PIXEL_GAINS_ON) {
    Int_t pixel = PixelAlignment::PixelIndexOf(x_mm, y_mm);
    Double_t gain = (pixel >= 0) ? PixelGains558()[pixel] : 0.0;
    return (gain > 0) ? energy / gain : energy;
  }
  if (!map)
    return energy;
  return energy * map->Factor(PixelAlignment::PixelIndexOf(x_mm, y_mm));
}

// At 8 the low-exp tail railed on every line. Up to 20 it settles inside
// (tau/sigma 9.5 Am, 10.5 Ba, 10.7 Pb in-situ) and the Pb Ka chi2 drops from
// 2.6 to 1.1 per degree of freedom.
const Double_t TAIL_RATIO_MAX = 20.0;

// Fixed when the macro loads, so production and its systematic partners
// each write their own calibrated trees, plots and results. The fixed parts
// of the name are where the saved fit states live: tau20 the tail bound,
// sigpb the signal-run Pb anchor, pxday the day pixel map (pxoff without it,
// nogain without the 558 keV gains either; their states are copies of
// production's, as are escup's and escdn's, the escape fractions moved).
const TString RESULT_TAG =
    TString("tau20_sigpb") + (EscapeLinesRequested() ? "_esc" : "") +
    (OffsetTransferRequested() ? "_oxfer" : "") +
    (EscapeFractionShift() > 0 ? "_escup"
                               : (EscapeFractionShift() < 0 ? "_escdn" : "")) +
    (!PIXEL_GAINS_ON ? "_nogain" : (PIXEL_MAP_ON ? "_pxday" : "_pxoff"));
const TString OUT_SUBDIR = "calibrated_" + RESULT_TAG;

// ------------------------------------------------------------------- events

inline std::vector<Double_t> LoadCalibrated(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading(OUT_SUBDIR + "/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << OUT_SUBDIR << "/" << input_name
              << ".root" << std::endl;
    return events;
  }
  TTree *tree = static_cast<TTree *>(file->Get("calibrated_tree"));
  if (tree)
    events = RooFitUtils::LoadEventsFromTree(tree, "depositedEnergykeV");
  file->Close();
  delete file;
  return events;
}

// Rate histograms of an already-calibrated event list over its live time, in
// bins of width bin at n_phases starting points bin / n_phases apart: the
// event-list counterpart of CalibrateAndFit.cpp's BuildRateHists. That one
// divides each crystal by its own live time; an event list pools the crystals
// over their summed live time, and within a run the crystals' live times agree
// to under 1%. Sumw2 so the propagated errors survive a subtraction.
inline std::vector<TH1D *>
RateHistsFromEvents(const std::vector<Double_t> &events, Double_t live_time,
                    const TString &name, Double_t lo, Double_t hi, Double_t bin,
                    Int_t n_phases) {
  std::vector<TH1D *> hists;
  Int_t nbins = (Int_t)std::lround((hi - lo) / bin);
  for (Int_t k = 0; k < n_phases; k++) {
    Double_t start = lo + k * bin / n_phases;
    TH1D *h = new TH1D(name + Form("_phase%d", k),
                       ";Deposited Energy [keV];Rate [counts/s]", nbins, start,
                       start + nbins * bin);
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (size_t i = 0; i < events.size(); i++)
      h->Fill(events[i]);
    if (live_time > 0)
      h->Scale(1.0 / live_time);
    hists.push_back(h);
  }
  return hists;
}

// Tail/Gaussian HEIGHT ratio, which the TF1 models need, from the tail's
// full-range yield ratio R = k * 2 rho exp(1/(2 rho^2)) / sqrt(2 pi). R is what
// RooFitUtils reports (low_exp_tail_amplitude / gaus_amplitude) and does not
// depend on the fit window.
inline Double_t TailHeightRatio(const PeakFitResult &pk) {
  if (!(pk.gaus_amplitude > 0) || !(pk.low_exp_tail_amplitude > 0) ||
      !(pk.low_exp_tail_ratio > 0))
    return 0.0;
  Double_t rho = pk.low_exp_tail_ratio;
  Double_t full_range_ratio = pk.low_exp_tail_amplitude / pk.gaus_amplitude;
  return full_range_ratio * TMath::Sqrt(2.0 * TMath::Pi()) /
         (2.0 * rho * TMath::Exp(1.0 / (2.0 * rho * rho)));
}

// The inverse of TailHeightRatio: the full-range tail yield ratio (RooFitUtils'
// LowExpTailAmplitude) of a low exponential tail with this height ratio to
// its Gaussian and this tau/sigma.
inline Double_t FullRangeTailRatio(Double_t tail_height_ratio,
                                   Double_t tau_over_sigma) {
  Double_t rho = tau_over_sigma;
  return tail_height_ratio * 2.0 * rho * TMath::Exp(1.0 / (2.0 * rho * rho)) /
         TMath::Sqrt(2.0 * TMath::Pi());
}

// ------------------------------------------------------- 95.902 keV escapes

// 113Cd(n,g) feeds a 95.902 keV transition in 114Cd, converted in the K shell
// about half the time. The Cd K X-rays of the conversions, and the Cd and Te
// K X-rays of the photon's absorption, escape the crystal and leave lines at
// 95.902 keV minus the X-ray energy: Cd K-alpha under Pb K-alpha2, Cd K-beta
// at 69.3-69.8 keV and Te K-alpha at 68.4-68.7 keV, beside and under Ge.
// Both methods fit them as fixed lines at the run's own 95.902 keV centroid
// minus the X-ray energy, which is where the beam background runs put them
// (the line itself reads up to 270 eV low), with the run's 95.902 keV
// Gaussian yield times the fractions below and the lineshape of the lines
// they sit among. The X-ray energies and series are ReferenceLines.hpp's,
// each line's share of its series Geant4's (below).
using ReferenceLines::ESCAPE_SERIES;
using ReferenceLines::ESCAPE_XRAY_KEV;
using ReferenceLines::N_ESCAPE_LINES;

// Escapes of each series over the 95.902 keV line's Gaussian yield.
struct EscapeFractions {
  Double_t series[4];
};

// In every run, shielded or not, the line comes from captures in the CZT's
// own cathode layer: its depth and pixel-ring profiles show no dependence on
// the Cd shield beyond what the unshielded runs show (Capture96Positions.cpp).
// So every run takes the escape fractions of captures inside the crystal,
// which CaptureEscapeFraction.cpp measures in the 01/15 beam background run
// (no shield, no Pb): an unbinned fit of both Cd escape pairs as one yield, Cd
// K-alpha held to Cd K-beta and the Te escapes to Cd K-alpha by the Geant4
// ratios below. The macro prints its values against these and writes them to
// results/escape_fractions.result.
const EscapeFractions ESCAPES_INSIDE = {{0.05154, 0.01404, 0.00275, 0.00083}};
// Their relative uncertainty: the fit's statistics (0.0016 on Cd K-beta) and
// the largest shift over its alternatives, the photopeak width and a quadratic
// continuum (0.0008).
const Double_t ESCAPES_INSIDE_REL_ERR = 0.125;
// Geant4's ratios for captures inside the crystal, families of the
// czt96_surface run (geant4/Escape96Plots.cpp prints them).
const Double_t GEANT4_INSIDE_CD_KA_OVER_CD_KB = 3.671;
const Double_t GEANT4_INSIDE_TE_KA_OVER_CD_KA = 0.0534;
const Double_t GEANT4_INSIDE_TE_KB_OVER_TE_KA = 0.300;
// Each line's share of its series, czt96_surface's ("shares in family"): Cd
// Ka1, Ka2, Kb1', Kb2', Te Ka1, Ka2, Kb1', Kb2'.
const Double_t GEANT4_INSIDE_SERIES_SHARE[N_ESCAPE_LINES] = {
    0.6569, 0.3431, 0.8387, 0.1613, 0.6584, 0.3416, 0.8206, 0.1794};

// The set every fit uses: ESCAPES_INSIDE, or under ESCAPE_FRACTIONS=up or
// down that set moved by its uncertainty.
inline EscapeFractions ShiftedInsideEscapes() {
  EscapeFractions shifted = ESCAPES_INSIDE;
  for (Int_t s = 0; s < 4; s++)
    shifted.series[s] *= 1.0 + EscapeFractionShift() * ESCAPES_INSIDE_REL_ERR;
  return shifted;
}
const EscapeFractions ESCAPES_USED = ShiftedInsideEscapes();

// Photons from captures in a Cd shield, onto the crystal's face (Geant4). No
// run uses them: the Cd-shield background runs' 95.902 keV line is too weak
// to tell them from ESCAPES_INSIDE (CaptureEscapeFraction.cpp compares the
// two). The testbed still injects them.
const EscapeFractions ESCAPES_OUTSIDE = {{0.0115, 0.0033, 0.0057, 0.0017}};

// keV per raw unit on 01/13 (the pairs' calibrations give 0.9882-0.9902). It
// only spaces the escapes from the 95.902 keV line in the raw
// pre-calibration fit, where 50 eV of spacing moves Ge after the signal-run
// anchor by 0.1 eV.
const Double_t RAW_KEV_PER_UNIT_0113 = 0.989;

// The 95.902 keV line with a fixed low tail on a linear continuum, per bin.
//   [0] Gaussian yield  [1] mu  [2] sigma  [3] tail/Gaussian height ratio
//   [4] tau/sigma  [5] continuum at 96  [6] continuum slope  [7] bin width
inline Double_t CaptureLineModel(Double_t *x, Double_t *p) {
  Double_t height = p[0] * p[7] / (p[2] * TMath::Sqrt(2.0 * TMath::Pi()));
  Double_t gaus[3] = {p[1], p[2], height};
  Double_t tail[6] = {p[1], p[2], height * p[3], p[4], 0.0, 0.0};
  return FittingFunctions::Gaussian(x, gaus) +
         FittingFunctions::LowTail(x, tail) + p[5] + p[6] * (x[0] - 96.0);
}

struct CaptureLine {
  // In the histogram's units times keV per bin width: events for a counts
  // histogram, counts/s for a rate histogram.
  Double_t gaus_yield = 0;
  // The fit's statistical error on gaus_yield; over phases, their mean.
  Double_t gaus_yield_err = 0;
  // The fitted Gaussian width; over phases, their mean.
  Double_t sigma = 0;
  Double_t mu = 0;
  Bool_t valid = kFALSE;
};

// The window the 95.902 keV line is fitted in, in the events' units. It holds
// the line in raw (96.3) and calibrated (95.6-95.9) energies alike, and stays
// put: a window that follows the highest bin moves with the binning.
const Double_t CAPTURE_FIT_LO = 91.5;
const Double_t CAPTURE_FIT_HI = 100.5;

// AU's valid flag also passes Minuit status 3-4 at EDM < 10 and a forced
// positive-definite covariance, and a refit of a saved state is always valid,
// so a result is used only from a clean minimum with an accurate covariance.
inline Bool_t FitConverged(const FitResult &r, const TString &label) {
  if (r.fit_status == 0 && r.cov_qual == 3)
    return kTRUE;
  std::cerr << "ERROR: " << label << " did not converge cleanly (status "
            << r.fit_status << ", covariance quality " << r.cov_qual << ", edm "
            << r.edm << "); the fit is not used" << std::endl;
  return kFALSE;
}

// A fitted value pinned at its limit is no measurement, whatever the fit
// status says.
inline Bool_t ParameterAtLimit(const TF1 &f, Int_t i) {
  Double_t lo = 0, hi = 0;
  f.GetParLimits(i, lo, hi);
  if (!(hi > lo))
    return kFALSE;
  Double_t margin = 1e-4 * (hi - lo);
  Double_t value = f.GetParameter(i);
  return value <= lo + margin || value >= hi - margin;
}

// The 95.902 keV line in h, with the tail of the lines the escapes will sit
// among. Counts take a Poisson likelihood; a rate-sub residual, which goes
// negative, a chi2.
inline CaptureLine FitCaptureLine(TH1 *h, Double_t tail_height_ratio,
                                  Double_t tau_over_sigma, Bool_t poisson,
                                  const TString &tag) {
  CaptureLine out;
  const Double_t lo = CAPTURE_FIT_LO, hi = CAPTURE_FIT_HI;
  Int_t first = h->FindBin(lo), last = h->FindBin(hi);
  Int_t peak_bin = first;
  for (Int_t b = first; b <= last; b++)
    if (h->GetBinContent(b) > h->GetBinContent(peak_bin))
      peak_bin = b;
  Double_t mu0 = h->GetBinCenter(peak_bin);
  Double_t bin = h->GetBinWidth(1);
  Double_t level =
      0.5 * (h->GetBinContent(first + 1) + h->GetBinContent(last - 1));
  Double_t sigma0 = 0.9;
  Double_t yield0 = (h->GetBinContent(peak_bin) - level) * sigma0 *
                    TMath::Sqrt(2.0 * TMath::Pi()) / bin;
  TF1 model("capture_line_" + tag, CaptureLineModel, lo, hi, 8);
  model.SetParameters(yield0, mu0, sigma0, tail_height_ratio, tau_over_sigma,
                      level, 0.0, bin);
  model.FixParameter(3, tail_height_ratio);
  model.FixParameter(4, tau_over_sigma);
  model.FixParameter(7, bin);
  model.SetParLimits(1, mu0 - 1.0, mu0 + 1.0);
  model.SetParLimits(2, 0.4, 1.5);
  TFitResultPtr fr = h->Fit(&model, poisson ? "S L Q R N" : "S Q R N");
  if (fr.Get() && fr->IsValid() && model.GetParameter(0) > 0 &&
      !ParameterAtLimit(model, 1) && !ParameterAtLimit(model, 2)) {
    out.gaus_yield = model.GetParameter(0);
    out.gaus_yield_err = model.GetParError(0);
    out.sigma = model.GetParameter(2);
    out.mu = model.GetParameter(1);
    out.valid = kTRUE;
  }
  return out;
}

// The line averaged over its fits in bins one comb period wide at
// N_RATESUB_PHASES starting points, which beat neither against the comb nor
// against where they start: in 255 eV bins at one starting point the Cu 01/13
// background run's yield moved by 20%.
inline CaptureLine AverageCaptureLine(const std::vector<TH1D *> &hists,
                                      Double_t tail_height_ratio,
                                      Double_t tau_over_sigma, Bool_t poisson,
                                      const TString &tag) {
  CaptureLine out;
  Double_t sum_yield = 0, sum_yield2 = 0, sum_mu = 0, sum_mu2 = 0;
  Double_t sum_yield_err = 0;
  Double_t sum_sigma = 0;
  Int_t n_valid = 0;
  for (size_t k = 0; k < hists.size(); k++) {
    CaptureLine line =
        FitCaptureLine(hists[k], tail_height_ratio, tau_over_sigma, poisson,
                       tag + TString::Format("_phase%d", (Int_t)k));
    if (!line.valid)
      continue;
    sum_yield += line.gaus_yield;
    sum_yield_err += line.gaus_yield_err;
    sum_sigma += line.sigma;
    sum_yield2 += line.gaus_yield * line.gaus_yield;
    sum_mu += line.mu;
    sum_mu2 += line.mu * line.mu;
    n_valid++;
  }
  if (n_valid > 0) {
    out.gaus_yield = sum_yield / n_valid;
    out.gaus_yield_err = sum_yield_err / n_valid;
    out.sigma = sum_sigma / n_valid;
    out.mu = sum_mu / n_valid;
    out.valid = kTRUE;
  }
  Double_t yield_rms = 0, mu_rms = 0;
  if (n_valid > 0) {
    yield_rms = TMath::Sqrt(TMath::Max(
        sum_yield2 / n_valid - out.gaus_yield * out.gaus_yield, 0.0));
    mu_rms = TMath::Sqrt(TMath::Max(sum_mu2 / n_valid - out.mu * out.mu, 0.0));
  }
  std::cout << "CAPTURE-LINE " << tag << ": " << std::fixed
            << std::setprecision(4) << out.mu << " (phase rms "
            << std::setprecision(1) << 1000.0 * mu_rms
            << " eV), Gaussian yield " << std::setprecision(2) << out.gaus_yield
            << " (phase rms " << yield_rms << "), " << n_valid << " of "
            << hists.size() << " phases" << std::endl;
  return out;
}

// The escapes of a capture line: centroids and Gaussian yields.
// units_per_kev turns an X-ray energy into a spacing in the line's units: 1
// for calibrated energies.
inline void EscapeLinesOf(const CaptureLine &line,
                          const EscapeFractions &fractions,
                          Double_t units_per_kev, std::vector<Double_t> &mus,
                          std::vector<Double_t> &gaus_yields) {
  mus.clear();
  gaus_yields.clear();
  for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
    mus.push_back(line.mu - ESCAPE_XRAY_KEV[i] * units_per_kev);
    gaus_yields.push_back(line.gaus_yield * fractions.series[ESCAPE_SERIES[i]] *
                          GEANT4_INSIDE_SERIES_SHARE[i]);
  }
}

// How a pair's escapes enter one of its in-situ fits.
struct InSituEscapes {
  const EscapeFractions *fractions = nullptr;
  // The events' units per keV, which spaces the escapes from the line, and one
  // comb period in those units, the bin width of the line's fits.
  Double_t units_per_kev = 1.0;
  Double_t comb_period = COMB_PERIOD_RAW_KEV;
};

// Fit a channel's capture line in its events with the tail of tail_source
// and add its escapes to that channel as fixed lines with the shape of
// shape_peak. False when the line cannot be measured: a fit asked for escapes
// must not go on without them.
inline Bool_t AddEscapeLinesToChannel(RooFitUtils &sim, const TString &channel,
                                      Int_t shape_peak,
                                      const std::vector<Double_t> &events,
                                      const PeakFitResult &tail_source,
                                      const InSituEscapes &escapes,
                                      const TString &tag) {
  const Double_t bin = escapes.comb_period;
  const Int_t n_bins = (Int_t)std::lround(16.0 / bin);
  std::vector<TH1D *> hists;
  for (Int_t k = 0; k < N_RATESUB_PHASES; k++) {
    Double_t start = 88.0 + k * bin / N_RATESUB_PHASES;
    TH1D *h = new TH1D(TString::Format("capture_hist_%s_%d", tag.Data(), k), "",
                       n_bins, start, start + n_bins * bin);
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < events.size(); i++)
      h->Fill(events[i]);
    hists.push_back(h);
  }
  Double_t tau_over_sigma =
      tail_source.low_exp_tail_ratio > 0 ? tail_source.low_exp_tail_ratio : 1.0;
  CaptureLine line = AverageCaptureLine(hists, TailHeightRatio(tail_source),
                                        tau_over_sigma, kTRUE, tag);
  for (size_t k = 0; k < hists.size(); k++)
    delete hists[k];
  if (!line.valid) {
    std::cerr << "ERROR: no 95.902 keV line for " << tag
              << "; the fit is abandoned" << std::endl;
    return kFALSE;
  }
  std::vector<Double_t> mus, gaus_yields;
  EscapeLinesOf(line, *escapes.fractions, escapes.units_per_kev, mus,
                gaus_yields);
  Double_t total = 0;
  for (size_t i = 0; i < gaus_yields.size(); i++)
    total += gaus_yields[i];
  std::cout << "ESCAPES " << tag << ": " << std::fixed << std::setprecision(1)
            << total << " events in 8 fixed lines, Cd K-alpha at "
            << std::setprecision(4) << mus[0] << " / " << mus[1]
            << ", Te K-alpha at " << mus[4] << " / " << mus[5] << std::endl;
  sim.AddFixedLines(channel, shape_peak, mus, gaus_yields);
  return kTRUE;
}

// The escapes per unit capture-line Gaussian yield, per keV at x: each line a
// Gaussian of area equal to its fraction, plus the low tail.
inline Double_t EscapeDensity(Double_t x, Double_t capture_mu, Double_t sigma,
                              Double_t tail_height_ratio,
                              Double_t tau_over_sigma,
                              const EscapeFractions &fractions) {
  Double_t xx[1] = {x};
  Double_t density = 0;
  for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
    Double_t mu = capture_mu - ESCAPE_XRAY_KEV[i];
    Double_t gaus[3] = {mu, sigma, 1.0};
    Double_t tail[6] = {mu, sigma, tail_height_ratio, tau_over_sigma, 0.0, 0.0};
    Double_t unit_line = (FittingFunctions::Gaussian(xx, gaus) +
                          FittingFunctions::LowTail(xx, tail)) /
                         (sigma * TMath::Sqrt(2.0 * TMath::Pi()));
    density += fractions.series[ESCAPE_SERIES[i]] *
               GEANT4_INSIDE_SERIES_SHARE[i] * unit_line;
  }
  return density;
}

// ------------------------------------------------------------ simultaneous

struct BkgGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid = kFALSE;
};

// The in-situ model's components. No high-side exponential tail: it never
// takes an interior value across the hand-fitted runs, either switching off or
// railing at both bounds at once.
// The continuum is flat: a sloped one takes up shape that belongs to the
// lines' tails, and its fits end with a forced covariance.
const Bool_t INSITU_FLAT_BKG = kTRUE;
const Bool_t INSITU_STEP = kFALSE;
const Bool_t INSITU_LOW_EXP = kTRUE;
const Bool_t INSITU_LOW_LIN = kFALSE;
const Bool_t INSITU_HIGH_EXP = kFALSE;

// One simultaneous fit of a pair, from its background-run seed. With escapes,
// each channel also carries its own run's 95.902 keV escapes as fixed lines
// with the Pb Ka1 shape, the line measured with tail_source's tail. The fit
// is unbinned; bin_kev only bins its plots and chi2, one comb period in the
// events' units so that each bin holds one tooth of the comb.
inline BkgGeSimResult FitBkgGeSimFromSeed(
    const std::vector<Double_t> &bkg_events,
    const std::vector<Double_t> &sig_events,
    const std::vector<Double_t> &bkg_peak_mus, const TString &sig_label,
    const TString &fit_label, Double_t bin_kev, const FitResult &seed,
    const FitResult *precomputed_sig_seed, const InSituEscapes *escapes,
    const PeakFitResult *tail_source) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();

  // The saved .simroofits carries the hand-tuned starting point a cold start
  // cannot reach; SetRefitAfterLoad makes it SEED a real fit rather than be
  // returned as the result, which once produced post-cal rows holding day-one
  // constants at chi2 7-9.
  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  sim.SetInteractive();
  sim.SetRefitAfterLoad();

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back((Double_t)E_NEAR_GE);
  std::vector<Bool_t> bkg_fixed(n_bkg, kFALSE);
  std::vector<Bool_t> sig_fixed(n_bkg + 1, kFALSE);

  // The slope is held only for the flat model; a linear one fits it.
  sim.AddChannel("bkg", bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH, bin_kev,
                 n_bkg, bkg_peak_mus, INSITU_FLAT_BKG, INSITU_STEP,
                 INSITU_LOW_EXP, INSITU_LOW_LIN, INSITU_HIGH_EXP, bkg_fixed,
                 kFALSE, INSITU_FLAT_BKG, kFALSE);
  sim.AddChannel("sig", sig_events, CAL_RANGE_LOW, CAL_RANGE_HIGH, bin_kev,
                 n_bkg + 1, sig_mus, INSITU_FLAT_BKG, INSITU_STEP,
                 INSITU_LOW_EXP, INSITU_LOW_LIN, INSITU_HIGH_EXP, sig_fixed,
                 kFALSE, INSITU_FLAT_BKG, kFALSE);
  // Share ONE low-side tail across both Pb-Ka lines and the Ge peak, tied to
  // bkg:Ka1. Fit independently, Ka1 inflates its own tail to absorb the
  // inter-peak fill so neither Ka tail is trustworthy, and the Ge tail rails
  // to zero, biasing its centroid low. Ge mu/sigma/yield stay free.
  // Pushed BEFORE LinkPeakShape so they win first-match on the tail params.
  // Route every target straight to bkg:Ka1 -- chaining via Ka2 hits a
  // ResolveOrCreate quirk that orphans sig:Ka2's tail into a free parameter.
  const Int_t n_tp = 2;
  const char *tp[n_tp] = {"LowExpTailAmplitude", "LowExpTailRatio"};
  TString ge_idx = TString::Format("%d", n_bkg + 1);
  for (Int_t t = 0; t < n_tp; t++) {
    sim.LinkParameter(TString("bkg:") + tp[t] + "2",
                      TString("bkg:") + tp[t] + "1");
    for (Int_t i = 1; i < n_bkg; i++)
      sim.LinkParameter(TString("sig:") + tp[t] + TString::Format("%d", i + 1),
                        TString("bkg:") + tp[t] + "1");
    sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                      TString("bkg:") + tp[t] + "1");
  }

  sim.LinkParameter("bkg:Sigma2", "bkg:Sigma1");

  // The signal run's Pb lines keep their own centroids: the two runs of a pair
  // sit at different rates, and the reading moves with rate by tens of eV.
  for (Int_t i = 0; i < n_bkg; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  if (escapes && tail_source) {
    if (!AddEscapeLinesToChannel(sim, "bkg", n_bkg - 1, bkg_events,
                                 *tail_source, *escapes, fit_label + "_bkg"))
      return out;
    if (!AddEscapeLinesToChannel(sim, "sig", n_bkg - 1, sig_events,
                                 *tail_source, *escapes, fit_label + "_sig"))
      return out;
  }

  sim.SeedChannel("bkg", seed);
  if (precomputed_sig_seed)
    sim.SeedChannel("sig", *precomputed_sig_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_label, fit_label);
  if (results.size() < 2)
    return out;
  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid && FitConverged(results[1], fit_label);
  return out;
}

// bkg channel: the Pb Ka2 and Ka1 peaks. sig channel: the same peaks plus the
// Ge line at E_NEAR_GE. One low-side tail and one Pb shape are shared across
// the two channels; see the comments in FitBkgGeSimFromSeed for why.
//
// With escapes, the 95.902 keV line is measured with a converged shared tail.
// A calibrated fit takes it from the raw fit it is seeded with. The raw fit's
// own seed is the background run alone, whose saved tail can sit at its
// bounds (at tau/sigma 1 it takes 30% of the line's Gaussian yield), so the
// raw fit first runs without the escapes and the line is measured with its
// tail.
inline BkgGeSimResult
RunBkgGeSim(const std::vector<Double_t> &bkg_events,
            const std::vector<Double_t> &sig_events,
            const std::vector<Double_t> &bkg_peak_mus, const TString &bkg_label,
            const TString &sig_label, const TString &fit_label,
            Double_t bin_kev, const FitResult *precomputed_bkg_seed = nullptr,
            const FitResult *precomputed_sig_seed = nullptr,
            const InSituEscapes &escapes = InSituEscapes()) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();
  if (n_bkg != 2 || bkg_events.empty() || sig_events.empty())
    return out;

  FitResult seed;
  if (precomputed_bkg_seed) {
    seed = *precomputed_bkg_seed;
  } else {
    RooFitUtils bkg_fitter(bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH, bin_kev,
                           INSITU_FLAT_BKG, INSITU_STEP, INSITU_LOW_EXP,
                           INSITU_LOW_LIN, INSITU_HIGH_EXP);
    bkg_fitter.SetTailRatioMax(TAIL_RATIO_MAX);
    bkg_fitter.SetInteractive();
    TString seed_label = fit_label + "_BkgSeed";
    seed = bkg_fitter.FitDoublePeak(bkg_label, seed_label, bkg_peak_mus[0],
                                    bkg_peak_mus[1]);
    if (!seed.valid) {
      std::cerr << "ERROR: bkg-only seed fit failed for " << bkg_label
                << std::endl;
      return out;
    }
  }

  if (!escapes.fractions)
    return FitBkgGeSimFromSeed(bkg_events, sig_events, bkg_peak_mus, sig_label,
                               fit_label, bin_kev, seed, precomputed_sig_seed,
                               nullptr, nullptr);
  if (precomputed_bkg_seed && !seed.peaks.empty())
    return FitBkgGeSimFromSeed(bkg_events, sig_events, bkg_peak_mus, sig_label,
                               fit_label, bin_kev, seed, precomputed_sig_seed,
                               &escapes, &seed.peaks.back());
  BkgGeSimResult without = FitBkgGeSimFromSeed(
      bkg_events, sig_events, bkg_peak_mus, sig_label, fit_label, bin_kev, seed,
      precomputed_sig_seed, nullptr, nullptr);
  if (!without.valid || without.bkg_channel.peaks.empty())
    return without;
  return FitBkgGeSimFromSeed(bkg_events, sig_events, bkg_peak_mus, sig_label,
                             fit_label, bin_kev, seed, precomputed_sig_seed,
                             &escapes, &without.bkg_channel.peaks.back());
}

// ------------------------------------------------------------- rate-sub fit

struct RateSubResult {
  Float_t mu = -1, mu_err = -1, chi2 = -1;
  Float_t sigma = -1, bkg0 = 0, slope = 0;
  Float_t tail_amp = -1, tail_ratio = -1; // -1 when no tail was fitted
  // Production's phase-averaged fit only: the spread over bin phases, the
  // shift to a quadratic residual background, and the two in quadrature.
  Float_t phase_rms = 0, shape_shift = 0, sys_err = 0;
  Bool_t valid = kFALSE;
};

// A Gaussian on a linear background fitted to an already-subtracted residual,
// the seed of the shaped fit below. Callers own the histogram; pass out_model
// to keep the TF1, which is otherwise deleted here.
//
// chi2 fit, not FittingUtils: its Poisson option assumes integer counts, but
// this histogram holds counts/second and goes negative after subtraction. The
// histogram must have Sumw2 set so propagated errors exist.
inline RateSubResult FitRateSubHist(TH1D *res, Double_t lo, Double_t hi,
                                    const TString &tag,
                                    TF1 **out_model = nullptr) {
  RateSubResult out;
  if (!res)
    return out;

  Double_t edge_lo = res->Integral(res->FindBin(lo + 0.2), res->FindBin(66.0));
  Double_t nedge = res->FindBin(66.0) - res->FindBin(lo + 0.2) + 1;
  Double_t bkg0 = (nedge > 0) ? edge_lo / nedge : 0.0;
  Int_t pk_bin = res->GetMaximumBin();
  Double_t pk_val = res->GetBinContent(pk_bin) - bkg0;
  if (pk_val < 0)
    pk_val = 0;
  Double_t sigma0 = 0.85;
  static Int_t model_counter = 0;
  TString mname = Form("ratesub_model_%s_%d", tag.Data(), ++model_counter);
  //   [0] area   [1] mu   [2] sigma   [3] bkg0   [4] bkg slope
  const char *kGausLinBkg =
      "[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*(x-68.75)";
  TF1 *model = new TF1(mname, kGausLinBkg, lo, hi);
  model->SetParNames("Area", "Mu", "Sigma", "Bkg0", "BkgSlope");
  model->SetParameters(pk_val * sigma0 * 2.5066, E_NEAR_GE, sigma0, bkg0, 0.0);
  model->SetParLimits(1, lo + 1.0, hi - 1.0);
  model->SetParLimits(2, 0.2, 3.0);
  TFitResultPtr fr = res->Fit(model, "S Q R N");
  Bool_t ok = fr.Get() && fr->IsValid() && model->GetNDF() > 0;
  if (ok) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.sigma = model->GetParameter(2);
    out.bkg0 = model->GetParameter(3);
    out.slope = model->GetParameter(4);
    out.valid = kTRUE;
  }
  if (out_model)
    *out_model = model;
  else
    delete model;
  return out;
}

// ------------------------------------------- rate-sub fit, shared lineshape

// The residual alone cannot fix the Ge tail: it is degenerate with the linear
// background. So the rate-sub fit borrows the tail from the same day's
// reference lines, which puts its centroid on the same definition as the
// calibration anchors and the in-situ fit.

// f_T / f_G: the fraction of a low exponential tail inside [lo, hi] over the
// fraction of its Gaussian, evaluated by RooFitUtils' own factor. A full-range
// ratio R times this is the tail's in-window yield ratio for that peak.
inline Double_t LowExpTailWindowFactor(Double_t lo, Double_t hi, Double_t mu,
                                       Double_t sigma, Double_t rho) {
  RooRealVar m("window_factor_mu", "", mu);
  RooRealVar s("window_factor_sigma", "", sigma);
  RooRealVar r("window_factor_rho", "", rho);
  RooRealVar a("window_factor_lo", "", lo);
  RooRealVar b("window_factor_hi", "", hi);
  RooExpTailWindowFactor factor("window_factor", "", m, s, r, a, b, kFALSE);
  return factor.getVal();
}

// What the rate-sub fit takes as given: the Ge tail and the residual's
// 95.902 keV line, whose escapes it carries with the fractions of where the
// pair's captures happen (none while escapes is null or capture_rate is 0).
struct RateSubShape {
  Double_t ge_tail_height_ratio = 0;
  Double_t ge_tau_over_sigma = 1;
  const EscapeFractions *escapes = nullptr;
  Double_t capture_rate = 0; // Gaussian yield, counts/s
  Double_t capture_mu = 0;
};

// Line-cal day: the Ge tail interpolated linearly in energy between the
// Am-241 and Ba-133 80.98 anchor fits. They are in raw units, but tau/sigma
// and the height ratio do not depend on the energy scale.
inline RateSubShape RateSubShapeFromLines(const PeakFitResult &am,
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

const Int_t RATESUB_N_PARAMS = 14;

// Ge Gaussian + tail on a linear (optionally quadratic) background, plus the
// fixed escapes of the residual's 95.902 keV line, with Ge's width and tail.
// Parameters:
//   [0] Ge height  [1] Ge mu  [2] Ge sigma  [3] Ge tail/gaus height ratio
//   [4] Ge tau/sigma  [5] bkg0  [6] bkg slope  [7] bkg curvature
//   [8] 95.902 keV Gaussian yield x bin width  [9] its centroid
//   [10-13] the escape fractions of the four series
inline Double_t RateSubEscapes(Double_t x, const Double_t *p) {
  if (p[8] == 0)
    return 0.0;
  EscapeFractions fractions = {{p[10], p[11], p[12], p[13]}};
  return p[8] * EscapeDensity(x, p[9], p[2], p[3], p[4], fractions);
}

inline Double_t RateSubShapedModel(Double_t *x, Double_t *p) {
  Double_t ge_gaus[3] = {p[1], p[2], p[0]};
  Double_t ge_tail[6] = {p[1], p[2], p[0] * p[3], p[4], 0.0, 0.0};
  Double_t u = x[0] - 68.75;
  return FittingFunctions::Gaussian(x, ge_gaus) +
         FittingFunctions::LowTail(x, ge_tail) + p[5] + p[6] * u +
         p[7] * u * u + RateSubEscapes(x[0], p);
}

// Fit the already-subtracted residual over the production window with the
// borrowed shape, seeded from the plain Gaussian fit. The residual
// background is linear unless quadratic is set.
inline RateSubResult FitRateSubHistShaped(TH1D *res, const RateSubShape &shape,
                                          const TString &tag,
                                          TF1 **out_model = nullptr,
                                          Bool_t quadratic = kFALSE) {
  RateSubResult out;
  if (!res)
    return out;
  TF1 *seed = nullptr;
  RateSubResult gaus =
      FitRateSubHist(res, RATESUB_LO, RATESUB_HI, tag + "_seed", &seed);
  if (!gaus.valid || !seed) {
    delete seed;
    return out;
  }
  const char *names[RATESUB_N_PARAMS] = {
      "GeHeight",         "GeMu",      "GeSigma",    "GeTailRatio",
      "GeTauRatio",       "Bkg0",      "BkgSlope",   "BkgCurvature",
      "CaptureYieldBins", "CaptureMu", "EscapeCdKa", "EscapeCdKb",
      "EscapeTeKa",       "EscapeTeKb"};
  static Int_t model_counter = 0;
  TString mname = Form("ratesub_shaped_%s_%d", tag.Data(), ++model_counter);
  TF1 *model = new TF1(mname, RateSubShapedModel, RATESUB_LO, RATESUB_HI,
                       RATESUB_N_PARAMS);
  for (Int_t i = 0; i < RATESUB_N_PARAMS; i++)
    model->SetParName(i, names[i]);
  model->SetParameter(7, 0.0);
  if (!quadratic)
    model->FixParameter(7, 0.0);
  Double_t capture_yield_bins =
      shape.escapes ? shape.capture_rate * res->GetBinWidth(1) : 0.0;
  model->FixParameter(8, capture_yield_bins);
  model->FixParameter(9, shape.capture_mu);
  for (Int_t s = 0; s < 4; s++)
    model->FixParameter(10 + s, shape.escapes ? shape.escapes->series[s] : 0.0);
  Double_t sigma0 = seed->GetParameter(2);
  model->SetParameter(0, seed->GetParameter(0) /
                             (sigma0 * TMath::Sqrt(2.0 * TMath::Pi())));
  model->SetParameter(1, seed->GetParameter(1));
  model->SetParameter(2, sigma0);
  model->FixParameter(3, shape.ge_tail_height_ratio);
  model->FixParameter(4, shape.ge_tau_over_sigma);
  model->SetParameter(5, seed->GetParameter(3));
  model->SetParameter(6, seed->GetParameter(4));
  model->SetParLimits(1, RATESUB_LO + 1.0, RATESUB_HI - 1.0);
  model->SetParLimits(2, 0.2, 3.0);
  delete seed;

  TFitResultPtr fr = res->Fit(model, "S Q R N");
  if (fr.Get() && fr->IsValid() && model->GetNDF() > 0 &&
      !ParameterAtLimit(*model, 1) && !ParameterAtLimit(*model, 2)) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.sigma = model->GetParameter(2);
    out.tail_amp = shape.ge_tail_height_ratio;
    out.tail_ratio = shape.ge_tau_over_sigma;
    out.bkg0 = model->GetParameter(5);
    out.slope = model->GetParameter(6);
    out.valid = kTRUE;
  }
  if (out_model)
    *out_model = model;
  else
    delete model;
  return out;
}

// The 95.902 keV line is 1.08 times stronger in the 01/15 signal run than live
// time predicts, so the residual keeps a share of it and of its escapes. Its
// Gaussian yield and centroid in the residual of rate histograms over 88-104
// keV at the rate-sub fit's starting points, with the Ge tail, set shape's
// capture line. False when the line cannot be measured.
inline Bool_t SetResidualCaptureLine(const std::vector<TH1D *> &sigs,
                                     const std::vector<TH1D *> &bkgs,
                                     const TString &tag, RateSubShape &shape) {
  if (sigs.size() != bkgs.size() || sigs.empty())
    return kFALSE;
  std::vector<TH1D *> residuals;
  for (size_t k = 0; k < sigs.size(); k++) {
    TH1D *res = static_cast<TH1D *>(sigs[k]->Clone(
        TString::Format("capture_residual_%s_%d", tag.Data(), (Int_t)k)));
    res->SetDirectory(nullptr);
    res->Add(bkgs[k], -1.0);
    residuals.push_back(res);
  }
  CaptureLine line =
      AverageCaptureLine(residuals, shape.ge_tail_height_ratio,
                         shape.ge_tau_over_sigma, kFALSE, "ratesub_" + tag);
  for (size_t k = 0; k < residuals.size(); k++)
    delete residuals[k];
  if (!line.valid) {
    std::cerr << "ERROR: no residual 95.902 keV line for " << tag
              << "; the rate-sub fit is abandoned" << std::endl;
    return kFALSE;
  }
  shape.capture_rate = line.gaus_yield;
  shape.capture_mu = line.mu;
  return kTRUE;
}

// The rate-sub fit of the residuals of pairs of rate histograms binned one comb
// period wide at N_RATESUB_PHASES starting points, averaged over the starting
// points. The linear residual background is the central model and the
// quadratic one its systematic, the mean shift over the starting points where
// both fits converge. The residual is overdispersed (chi2/ndf
// 1.4-1.9 on 01/15-16: the pixels' discrete energy grids, which the two runs
// fill differently), so the fit error is scaled by sqrt(chi2/ndf). The first
// converged starting point's model and residual are handed back for drawing
// when the pointers are given (the caller then owns them).
inline RateSubResult FitRateSubPhases(const std::vector<TH1D *> &sigs,
                                      const std::vector<TH1D *> &bkgs,
                                      const RateSubShape &shape,
                                      const TString &tag,
                                      TF1 **out_model = nullptr,
                                      TH1D **out_residual = nullptr) {
  RateSubResult out;
  if (sigs.size() != bkgs.size() || sigs.empty())
    return out;
  Double_t sum_err = 0, sum_chi2 = 0, sum_quadratic_shift = 0;
  std::vector<Double_t> linear_mus;
  Int_t n_quadratic = 0;
  RateSubResult fit;
  TF1 *model = nullptr;
  TH1D *res = nullptr;
  for (Int_t k = 0; k < (Int_t)sigs.size(); k++) {
    TString phase_tag = tag + Form("_phase%d", k);
    TH1D *phase_res =
        static_cast<TH1D *>(sigs[k]->Clone("rate_residual_" + phase_tag));
    phase_res->SetDirectory(nullptr);
    phase_res->Add(bkgs[k], -1.0);
    TF1 *phase_model = nullptr;
    RateSubResult linear =
        FitRateSubHistShaped(phase_res, shape, phase_tag, &phase_model, kFALSE);
    RateSubResult quadratic = FitRateSubHistShaped(
        phase_res, shape, phase_tag + "_quadratic", nullptr, kTRUE);
    std::cout << "RATESUB-PHASE " << tag << " " << k << ": linear "
              << std::fixed << std::setprecision(4) << linear.mu << " +/- "
              << linear.mu_err << " (chi2/ndf " << std::setprecision(2)
              << linear.chi2 << ")   quadratic " << std::setprecision(4)
              << quadratic.mu << " +/- " << quadratic.mu_err << std::endl;
    if (linear.valid) {
      linear_mus.push_back(linear.mu);
      sum_err += linear.mu_err;
      sum_chi2 += linear.chi2;
    }
    if (linear.valid && quadratic.valid) {
      sum_quadratic_shift += (Double_t)quadratic.mu - (Double_t)linear.mu;
      n_quadratic++;
    }
    if (!model && linear.valid && phase_model) {
      fit = linear;
      model = phase_model;
      res = phase_res;
    } else {
      delete phase_model;
      delete phase_res;
    }
  }
  const Int_t n_linear = (Int_t)linear_mus.size();
  if (n_linear == 0 || !model || n_quadratic == 0) {
    std::cerr << "ERROR: rate-sub fit failed for " << tag << " (" << n_linear
              << " linear, " << n_quadratic << " quadratic of " << sigs.size()
              << " starting points)" << std::endl;
    delete model;
    delete res;
    return out;
  }

  Double_t mean = 0;
  for (Int_t k = 0; k < n_linear; k++)
    mean += linear_mus[k];
  mean /= n_linear;
  // The spread from the deviations: E^2 - mean^2 at 68.77 keV loses the
  // few-eV spread to rounding.
  Double_t sum_dev2 = 0;
  for (Int_t k = 0; k < n_linear; k++)
    sum_dev2 += (linear_mus[k] - mean) * (linear_mus[k] - mean);
  Double_t chi2 = sum_chi2 / n_linear;
  out = fit;
  out.mu = mean;
  out.chi2 = chi2;
  out.mu_err = (sum_err / n_linear) * TMath::Sqrt(TMath::Max(chi2, 1.0));
  out.phase_rms = TMath::Sqrt(sum_dev2 / n_linear);
  out.shape_shift =
      (n_quadratic > 0) ? TMath::Abs(sum_quadratic_shift / n_quadratic) : 0.0;
  out.sys_err = TMath::Sqrt(out.phase_rms * out.phase_rms +
                            out.shape_shift * out.shape_shift);
  std::cout << "RATESUB-FIT " << tag << ": mu = " << std::fixed
            << std::setprecision(4) << out.mu << " +/- " << out.mu_err
            << " (fit, scaled by sqrt(chi2/ndf " << std::setprecision(2) << chi2
            << ")) +/- " << std::setprecision(4) << out.sys_err
            << " (sys: bin phase " << out.phase_rms << ", quadratic "
            << out.shape_shift << ") keV   sigma = " << out.sigma << std::endl;
  if (out_model)
    *out_model = model;
  else
    delete model;
  if (out_residual)
    *out_residual = res;
  else
    delete res;
  return out;
}

// ----------------------------------------------- Method 2: the no-shield run

// The Ge line fitted in the 01/15 no-shield signal run alone, without its
// background run. The structure near Ge that once made a background
// subtraction necessary is the 95.902 keV line's Cd and Te K X-ray escapes,
// which the fit carries: an unbinned likelihood over DIRECT_WINDOW_LO-HI of the
// run's calibrated events, the Ge line on a linear continuum (the no-shield
// continuum is not flat here) with the escapes as fixed lines of the Ge peak's
// shape. The Ge tail is held at the day's source-run tail, the Am-241 and
// Ba-133 80.998 keV lines interpolated to Ge. The escapes sit at the 95.902
// keV line's centroid minus the SRD 128 X-ray energies, with yields of the
// day's escape fractions times its Gaussian yield, the line measured in the
// same run on comb-period bins with the Ge tail. The same fit with Ge's tail
// free, started from the held one, gives the tail systematic: both the source
// lines' tail and Ge's own can be defended. The fit starts cold; it needs no
// saved seed.
const Double_t DIRECT_WINDOW_LO = 62.0;
const Double_t DIRECT_WINDOW_HI = 78.0;
// The 95.902 keV line's window, as AverageCaptureLine's callers bin it.
const Double_t DIRECT_CAPTURE_LO = 88.0;
const Double_t DIRECT_CAPTURE_WIDTH = 16.0;
// RooFitUtils draws 0.85 times a window's low edge to 1.15 times its high
// edge; the events loaded cover that and the capture line's window.
const Double_t DIRECT_DISPLAY_LOW_FACTOR = 0.85;
const Double_t DIRECT_DISPLAY_HIGH_FACTOR = 1.15;

// The run's events the direct fit is given: its window as RooFitUtils draws
// it, and the capture line's window.
inline Bool_t InDirectFitRange(Double_t energy) {
  const Double_t lo = DIRECT_DISPLAY_LOW_FACTOR * DIRECT_WINDOW_LO;
  const Double_t hi = TMath::Max(DIRECT_DISPLAY_HIGH_FACTOR * DIRECT_WINDOW_HI,
                                 DIRECT_CAPTURE_LO + DIRECT_CAPTURE_WIDTH);
  return energy >= lo && energy <= hi;
}

// The run's 95.902 keV line on comb-period bins at N_RATESUB_PHASES starting
// points, with the Ge tail.
inline CaptureLine MeasureDirectCaptureLine(const std::vector<Double_t> &events,
                                            Double_t bin_kev,
                                            const RateSubShape &shape,
                                            const TString &tag) {
  const Int_t n_bins = (Int_t)std::lround(DIRECT_CAPTURE_WIDTH / bin_kev);
  std::vector<TH1D *> hists;
  for (Int_t k = 0; k < N_RATESUB_PHASES; k++) {
    Double_t start = DIRECT_CAPTURE_LO + k * bin_kev / N_RATESUB_PHASES;
    TH1D *h = new TH1D(Form("direct_capture_%s_%d", tag.Data(), k), "", n_bins,
                       start, start + n_bins * bin_kev);
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < events.size(); i++)
      h->Fill(events[i]);
    hists.push_back(h);
  }
  CaptureLine line =
      AverageCaptureLine(hists, shape.ge_tail_height_ratio,
                         shape.ge_tau_over_sigma, kTRUE, "direct_" + tag);
  for (size_t k = 0; k < hists.size(); k++)
    delete hists[k];
  return line;
}

// One direct fit: the Ge tail held at the source-run tail, or free. A seed,
// when given, starts the fit.
inline FitResult FitDirectGeOnce(const std::vector<Double_t> &events,
                                 Double_t bin_kev, const RateSubShape &shape,
                                 const std::vector<Double_t> &escape_mus,
                                 const std::vector<Double_t> &escape_yields,
                                 Bool_t hold_tail, const FitResult *seed,
                                 const TString &input_name,
                                 const TString &label) {
  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  std::vector<Double_t> mus(1, (Double_t)E_NEAR_GE);
  const Bool_t flat_continuum = kFALSE;
  sim.AddChannel("direct", events, DIRECT_WINDOW_LO, DIRECT_WINDOW_HI, bin_kev,
                 1, mus, flat_continuum, kFALSE, kTRUE, kFALSE, kFALSE);
  if (hold_tail) {
    sim.HoldParameter("direct:LowExpTailAmplitude1",
                      FullRangeTailRatio(shape.ge_tail_height_ratio,
                                         shape.ge_tau_over_sigma));
    sim.HoldParameter("direct:LowExpTailRatio1", shape.ge_tau_over_sigma);
  }
  if (!escape_mus.empty())
    sim.AddFixedLines("direct", 0, escape_mus, escape_yields);
  if (seed && seed->valid)
    sim.SeedChannel("direct", *seed);
  std::vector<FitResult> results = sim.FitSimultaneous(input_name, label);
  if (results.empty())
    return FitResult();
  // Per-bin data, model and pull, for the residual check and the figure.
  sim.DumpChannelCSV("direct",
                     PlottingUtils::GetPlotsBaseDir() + "/fits/" + label);
  FitResult fit = results[0];
  fit.valid = fit.valid && !fit.peaks.empty() && fit.fit_status == 0 &&
              fit.cov_qual == 3;
  return fit;
}

inline void PrintDirectFit(const TString &label, const FitResult &fit) {
  if (!fit.valid) {
    std::cout << "DIRECT-FIT " << label << ": failed (status " << fit.fit_status
              << ", covQual " << fit.cov_qual << ")" << std::endl;
    return;
  }
  const PeakFitResult &ge = fit.peaks[0];
  // The tail as RooFitUtils' full-range ratio to the Gaussian yield.
  Double_t tail_ratio = (ge.gaus_amplitude > 0)
                            ? ge.low_exp_tail_amplitude / ge.gaus_amplitude
                            : 0.0;
  std::cout << "DIRECT-FIT " << label << ": Ge " << std::fixed
            << std::setprecision(4) << ge.mu << " +/- " << ge.mu_error
            << " keV, sigma " << ge.sigma << ", Gaussian yield "
            << std::setprecision(0) << ge.gaus_amplitude << ", tail ratio "
            << std::setprecision(4) << tail_ratio << " tau/sigma "
            << std::setprecision(3) << ge.low_exp_tail_ratio << "   chi2/ndf "
            << fit.reduced_chi2 << std::endl;
}

// Where production leaves the day's held Ge tail for the testbed, which has no
// source-run fits of its own: beside the direct fit's other files.
inline TString DirectTailPath(const TString &plots_dir,
                              const TString &date_label) {
  return plots_dir + "/fits/Direct_" + date_label + "_held_tail.txt";
}

inline Bool_t WriteDirectTail(const TString &path, const RateSubShape &shape) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    std::cerr << "ERROR: cannot write " << path << std::endl;
    return kFALSE;
  }
  out << std::setprecision(17);
  out << "ge_tail_height_ratio " << shape.ge_tail_height_ratio << std::endl;
  out << "ge_tau_over_sigma " << shape.ge_tau_over_sigma << std::endl;
  return kTRUE;
}

inline Bool_t ReadDirectTail(const TString &path, RateSubShape &shape) {
  std::ifstream in(path.Data());
  std::string name;
  Double_t value = 0;
  Int_t found = 0;
  while (in >> name >> value) {
    if (name == "ge_tail_height_ratio") {
      shape.ge_tail_height_ratio = value;
      found++;
    } else if (name == "ge_tau_over_sigma") {
      shape.ge_tau_over_sigma = value;
      found++;
    }
  }
  return found == 2;
}

#endif // GAINTRANSFER_FIT_HPP
