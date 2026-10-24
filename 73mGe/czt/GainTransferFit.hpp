#ifndef GAINTRANSFER_FIT_HPP
#define GAINTRANSFER_FIT_HPP

// The fitting core of the CZT 73mGe analysis.
//
// CalibrateAndFit.cpp (the measurement) and InjectionTest.cpp (its testbed
// closure) both include this, so the testbed exercises the same in-situ and
// rate-subtracted fits that produce the result.

#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "RooFitPhotopeakPdfs.hpp"
#include "RooFitUtils.hpp"

#include <RooRealVar.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------- constants

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_79 = 79.6142;
const Float_t E_BA133_81 = 80.9979;
// ENSDF intensities, 2.65(5) and 32.9(3) per 100 decays.
const Double_t BA133_79_OVER_81 = 2.65 / 32.9;
const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;

const Float_t E_NEAR_GE = 69;

// Anchor reference-energy uncertainties (keV), folded into the cal-fit x-errors
// so they propagate into the pol1 parameter covariance and thus the Ge energy.
// Pb Ka are known to <1 eV (Deslattes).
//
// dE_LINE remains the conservative ~10 eV lump, now used ONLY for the Ba-133
// lines of the 15th/16th line-cal days.
// TODO: replace with per-line evaluated Ba-133 uncertainties.
const Float_t dE_PB_KA = 0.001;
const Float_t dE_AM241 = 0.0001;
const Float_t dE_LINE = 0.010;

const Float_t CAL_RANGE_LOW = 65;
const Float_t CAL_RANGE_HIGH = 80;
const Float_t RATESUB_LO = 62.0;
const Float_t RATESUB_HI = 71.5;

// The summed spectra carry a comb with this period in raw energy, the same on
// 01/13, 01/15 and 01/16. Rate-sub bins one period wide (times the day's
// slope) hold one tooth each wherever they start; 255 eV bins beat against it
// and move Ge by 13-16 eV rms with where the bins start. The fit is averaged
// over N_RATESUB_PHASES starting points, and their spread is a systematic.
const Double_t COMB_PERIOD_RAW_KEV = 0.264;
const Int_t N_RATESUB_PHASES = 5;

// The 01/15-16 calibration comes from source runs and is carried to the beam
// runs, and no line appears in both to check it. On 01/13 the same step, Am
// source run to beam run, holds to about 20 eV: the beam run's own Pb Ka1 ->
// Kb1 slope matches the Am + Pb Ka slope.
const Double_t RATESUB_TRANSFER_SYS_KEV = 0.020;

const Bool_t LINK_PB_KA_SIGMA = kTRUE;

// BKG=linear: a sloped continuum in the in-situ fit instead of the flat one, a
// check of the continuum shape.
inline Bool_t LinearBkgRequested() {
  const char *env = getenv("BKG");
  return env && TString(env) == "linear";
}

const Bool_t USE_FLAT_BKG = !LinearBkgRequested();

// FIT_HARD=1: Minuit strategy 2, a final HESSE and 200k calls in RunBkgGeSim.
inline Bool_t HardFitRequested() {
  const char *env = getenv("FIT_HARD");
  return env && *env && atoi(env) != 0;
}

// Each in-situ signal run is calibrated on its own Pb Ka lines rather than
// inheriting its background run's: Ge follows its own run's Pb, and the Cd 10%
// signal run reads +36 +/- 11 eV above its background run within the same
// crystal, pixel ring and depth slice, a shift of the whole run.
// PB_ANCHOR=background restores the background-run anchor.
inline Bool_t SignalPbAnchorRequested() {
  const char *env = getenv("PB_ANCHOR");
  return !(env && TString(env) == "background");
}

// TRANSFER=offset: the pairs without their own Am run inherit the reference
// calibration shifted by their Pb K-alpha offset from the reference run, not
// scaled by the Pb ratio. At Ge the two differ by ~7% of the Pb shift. No line
// inside the calibrated range resolves them (Pb K-beta lies above it), so the
// choice is a systematic, not a measurement.
inline Bool_t OffsetTransferRequested() {
  const char *env = getenv("TRANSFER");
  return env && TString(env) == "offset";
}

// FIDUCIAL=interior: every run, calibration and source runs included, keeps
// only the 9 x 9 interior pixels of each crystal, dropping its outer ring and
// the ring along the gap between crystals. Pb K-alpha1 reads 7-155 eV away
// from its all-pixel value in those rings.
inline Bool_t InteriorPixelsRequested() {
  const char *env = getenv("FIDUCIAL");
  return env && TString(env) == "interior";
}

// Pixel centres sit at |x|, |y| = 1.743 ... 20.769 mm on a 1.9 mm pitch.
inline Bool_t InInteriorPixels(Float_t x_mm, Float_t y_mm) {
  const Float_t inner_edge = Constants::PIXEL_CENTERS_X_MM[11] + 0.5;
  const Float_t outer_edge = Constants::PIXEL_CENTERS_X_MM[21] - 0.5;
  Float_t ax = TMath::Abs(x_mm), ay = TMath::Abs(y_mm);
  return ax > inner_edge && ax < outer_edge && ay > inner_edge &&
         ay < outer_edge;
}

// The 95.902 keV line's K X-ray escapes enter both methods' fits as fixed
// lines (see "95.902 keV escapes" below). ESCAPES=off leaves them out.
inline Bool_t EscapeLinesRequested() {
  const char *env = getenv("ESCAPES");
  return !(env && TString(env) == "off");
}

// At 8 the low-exp tail railed on every line. Up to 20 it settles inside
// (tau/sigma 9.5 Am, 10.5 Ba, 10.7 Pb in-situ) and the Pb Ka chi2 drops from
// 2.6 to 1.1 per degree of freedom.
const Double_t TAIL_RATIO_MAX = 20.0;

// Fixed when the macro loads, from the switches above, so each
// configuration writes its own calibrated trees, plots and results. The tag
// names the Pb anchor and the escapes even though both are the default:
// production is tau20_sigpb_esc, tau20_sigpb is it without the escapes, and
// tau20 is the background-run anchor.
const TString RESULT_TAG = TString("tau20") +
                           (HardFitRequested() ? "_hard" : "") +
                           (SignalPbAnchorRequested() ? "_sigpb" : "") +
                           (EscapeLinesRequested() ? "_esc" : "") +
                           (OffsetTransferRequested() ? "_oxfer" : "") +
                           (InteriorPixelsRequested() ? "_int" : "") +
                           (LinearBkgRequested() ? "_linbkg" : "");
const TString OUT_SUBDIR = "calibrated_" + RESULT_TAG;

// Which result set this run writes. RESULT_TAG names the default one;
// RESULT_TAG_OVERRIDE lets the same macro produce a second set (e.g. a
// harder-convergence re-run) instead of keeping a duplicate source file.
// It changes OUTPUT paths only -- the calibrated inputs always come from
// OUT_SUBDIR, because overriding the tag must not change which calibrated
// trees are read, only where the results land.
inline TString RuntimeResultTag() {
  const char *env = getenv("RESULT_TAG_OVERRIDE");
  if (env && *env)
    return TString(env);
  return RESULT_TAG;
}

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
// event-list counterpart of CalibrateAndFit.cpp's BuildRateHists. Sumw2 so the
// propagated errors survive a subtraction.
inline std::vector<TH1F *>
RateHistsFromEvents(const std::vector<Double_t> &events, Double_t live_time,
                    const TString &name, Double_t lo, Double_t hi, Double_t bin,
                    Int_t n_phases) {
  std::vector<TH1F *> hists;
  Int_t nbins = (Int_t)std::lround((hi - lo) / bin);
  for (Int_t k = 0; k < n_phases; k++) {
    Double_t start = lo + k * bin / n_phases;
    TH1F *h = new TH1F(name + Form("_phase%d", k),
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
// they sit among.
const Int_t N_ESCAPE_LINES = 8;
// Cd Ka1, Ka2, Kb1', Kb2', Te Ka1, Ka2, Kb1', Kb2', in keV; Kb1' and Kb2' at
// their groups' intensity-weighted energies.
const Double_t ESCAPE_XRAY_KEV[N_ESCAPE_LINES] = {
    23.1739, 22.9843, 26.092, 26.650, 27.4726, 27.202, 30.988, 31.720};
// Series of each line: 0 Cd K-alpha, 1 Cd K-beta, 2 Te K-alpha, 3 Te K-beta.
const Int_t ESCAPE_SERIES[N_ESCAPE_LINES] = {0, 0, 1, 1, 2, 2, 3, 3};
// Each line's share of its series: DDEP emission intensities times the
// probability of escaping through the face of a slab.
const Double_t ESCAPE_SERIES_SHARE[N_ESCAPE_LINES] = {
    1.0 / 1.5242, 0.5242 / 1.5242, 0.3394 / 0.4035, 0.0641 / 0.4035,
    1.0 / 1.5275, 0.5275 / 1.5275, 0.3548 / 0.4349, 0.0801 / 0.4349};

// Escapes of each series over the 95.902 keV line's Gaussian yield.
struct EscapeFractions {
  Double_t series[4];
};

// Captures inside the crystal (no shield, Cu shield). Cd K-beta, the one
// escape clear of other lines, is measured in the 01/15 and 01/16 beam
// background runs (0.0122 and 0.0115, each +/- ~0.0017), and Cd K-alpha
// follows from the shares. Te K-alpha takes Geant4's 0.053 of Cd K-alpha for
// captures inside the crystal (those runs give 0.0006 and 0.0021, +/-
// 0.002-0.004), Te K-beta Geant4's 0.30 of Te K-alpha.
const EscapeFractions ESCAPES_INSIDE = {{0.0450, 0.0119, 0.0024, 0.0007}};

// Photons from captures outside the crystal, in the Cd shield, which stops
// the thermal neutrons: Geant4, the line's photons onto the crystal's face.
// The 01/13 Cd background runs show no Cd K-beta escape (0 +/- 0.004 at 10%,
// +/- 0.006 at 25%) where captures inside would give 0.012.
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
  Double_t mu = 0;
  Bool_t valid = kFALSE;
};

// The window the 95.902 keV line is fitted in, in the events' units. It holds
// the line in raw (96.3) and calibrated (95.6-95.9) energies alike, and stays
// put: a window that follows the highest bin moves with the binning.
const Double_t CAPTURE_FIT_LO = 91.5;
const Double_t CAPTURE_FIT_HI = 100.5;

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
  if (fr.Get() && fr->IsValid() && model.GetParameter(0) > 0) {
    out.gaus_yield = model.GetParameter(0);
    out.mu = model.GetParameter(1);
    out.valid = kTRUE;
  }
  return out;
}

// The line averaged over its fits in bins one comb period wide at
// N_RATESUB_PHASES starting points, which beat neither against the comb nor
// against where they start: in 255 eV bins at one starting point the Cu 01/13
// background run's yield moved by 20%.
inline CaptureLine AverageCaptureLine(const std::vector<TH1F *> &hists,
                                      Double_t tail_height_ratio,
                                      Double_t tau_over_sigma, Bool_t poisson,
                                      const TString &tag) {
  CaptureLine out;
  Double_t sum_yield = 0, sum_yield2 = 0, sum_mu = 0, sum_mu2 = 0;
  Int_t n_valid = 0;
  for (size_t k = 0; k < hists.size(); k++) {
    CaptureLine line =
        FitCaptureLine(hists[k], tail_height_ratio, tau_over_sigma, poisson,
                       tag + TString::Format("_phase%d", (Int_t)k));
    if (!line.valid)
      continue;
    sum_yield += line.gaus_yield;
    sum_yield2 += line.gaus_yield * line.gaus_yield;
    sum_mu += line.mu;
    sum_mu2 += line.mu * line.mu;
    n_valid++;
  }
  if (n_valid > 0) {
    out.gaus_yield = sum_yield / n_valid;
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
                          ESCAPE_SERIES_SHARE[i]);
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
// shape_peak.
inline void AddEscapeLinesToChannel(RooFitUtils &sim, const TString &channel,
                                    Int_t shape_peak,
                                    const std::vector<Double_t> &events,
                                    const PeakFitResult &tail_source,
                                    const InSituEscapes &escapes,
                                    const TString &tag) {
  const Double_t bin = escapes.comb_period;
  const Int_t n_bins = (Int_t)std::lround(16.0 / bin);
  std::vector<TH1F *> hists;
  for (Int_t k = 0; k < N_RATESUB_PHASES; k++) {
    Double_t start = 88.0 + k * bin / N_RATESUB_PHASES;
    TH1F *h = new TH1F(TString::Format("capture_hist_%s_%d", tag.Data(), k), "",
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
    std::cerr << "WARNING: no 95.902 keV line for " << tag
              << "; its escapes are left out" << std::endl;
    return;
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
    density +=
        fractions.series[ESCAPE_SERIES[i]] * ESCAPE_SERIES_SHARE[i] * unit_line;
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
const Bool_t INSITU_STEP = kFALSE;
const Bool_t INSITU_LOW_EXP = kTRUE;
const Bool_t INSITU_LOW_LIN = kFALSE;
const Bool_t INSITU_HIGH_EXP = kFALSE;

// One simultaneous fit of a pair, from its background-run seed. With escapes,
// each channel also carries its own run's 95.902 keV escapes as fixed lines
// with the Pb Ka1 shape, the line measured with tail_source's tail.
inline BkgGeSimResult FitBkgGeSimFromSeed(
    const std::vector<Double_t> &bkg_events,
    const std::vector<Double_t> &sig_events,
    const std::vector<Double_t> &bkg_peak_mus, const TString &sig_label,
    const TString &fit_label, const FitResult &seed,
    const FitResult *precomputed_sig_seed, const InSituEscapes *escapes,
    const PeakFitResult *tail_source) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();
  const Bool_t kFlatBkg = USE_FLAT_BKG;

  // The saved .simroofits carries the hand-tuned starting point a cold start
  // cannot reach; SetRefitAfterLoad makes it SEED a real fit rather than be
  // returned as the result, which once produced post-cal rows holding day-one
  // constants at chi2 7-9.
  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  sim.SetInteractive();
  sim.SetRefitAfterLoad();
  if (HardFitRequested()) {
    sim.SetFitStrategy(2);
    sim.SetFitFinalHesse(kTRUE);
    sim.SetFitMaxCalls(200000);
  }

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back((Double_t)E_NEAR_GE);
  std::vector<Bool_t> bkg_fixed(n_bkg, kFALSE);
  std::vector<Bool_t> sig_fixed(n_bkg + 1, kFALSE);

  std::vector<Bool_t> sig_step(n_bkg + 1, kTRUE);
  sig_step[n_bkg] = kFALSE;

  // The slope is held only for the flat model; a linear one fits it.
  sim.AddChannel("bkg", bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                 Constants::BIN_WIDTH_KEV, n_bkg, bkg_peak_mus, kFlatBkg,
                 INSITU_STEP, INSITU_LOW_EXP, INSITU_LOW_LIN, INSITU_HIGH_EXP,
                 bkg_fixed, kFALSE, kFlatBkg, kFALSE);
  sim.AddChannel("sig", sig_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                 Constants::BIN_WIDTH_KEV, n_bkg + 1, sig_mus, kFlatBkg,
                 INSITU_STEP, INSITU_LOW_EXP, INSITU_LOW_LIN, INSITU_HIGH_EXP,
                 sig_fixed, kFALSE, kFlatBkg, kFALSE, sig_step);
  // Share ONE low-side tail (LowExp + LowLin) across both Pb-Ka lines and the
  // Ge peak, tied to bkg:Ka1. Fit independently, Ka1 inflates its own tail to
  // absorb the inter-peak fill so neither Ka tail is trustworthy, and the Ge
  // tail rails to zero, biasing its centroid low. Ge mu/sigma/yield stay free.
  // Pushed BEFORE LinkPeakShape so they win first-match on the tail params.
  // Route every target straight to bkg:Ka1 -- chaining via Ka2 hits a
  // ResolveOrCreate quirk that orphans sig:Ka2's tail into a free parameter.
  // Only link components that exist: with kLowLin off the LowLinTail params
  // are never built and linking to them aborts the channel build.
  const char *tp[4] = {"LowExpTailAmplitude", "LowExpTailRatio",
                       "LowLinTailAmplitude", "LowLinTailSlope"};
  const Int_t n_tp = INSITU_LOW_LIN ? 4 : 2;
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

  if (LINK_PB_KA_SIGMA)
    sim.LinkParameter("bkg:Sigma2", "bkg:Sigma1");

  // The signal run's Pb lines keep their own centroids: the two runs of a pair
  // sit at different rates, and the reading moves with rate by tens of eV.
  for (Int_t i = 0; i < n_bkg; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  if (escapes && tail_source) {
    AddEscapeLinesToChannel(sim, "bkg", n_bkg - 1, bkg_events, *tail_source,
                            *escapes, fit_label + "_bkg");
    AddEscapeLinesToChannel(sim, "sig", n_bkg - 1, sig_events, *tail_source,
                            *escapes, fit_label + "_sig");
  }

  sim.SeedChannel("bkg", seed);
  if (precomputed_sig_seed)
    sim.SeedChannel("sig", *precomputed_sig_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_label, fit_label);
  if (results.size() < 2)
    return out;
  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;
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
            const FitResult *precomputed_bkg_seed = nullptr,
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
    RooFitUtils bkg_fitter(bkg_events, CAL_RANGE_LOW, CAL_RANGE_HIGH,
                           Constants::BIN_WIDTH_KEV, USE_FLAT_BKG, INSITU_STEP,
                           INSITU_LOW_EXP, INSITU_LOW_LIN, INSITU_HIGH_EXP);
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
                               fit_label, seed, precomputed_sig_seed, nullptr,
                               nullptr);
  if (precomputed_bkg_seed && !seed.peaks.empty())
    return FitBkgGeSimFromSeed(bkg_events, sig_events, bkg_peak_mus, sig_label,
                               fit_label, seed, precomputed_sig_seed, &escapes,
                               &seed.peaks.back());
  BkgGeSimResult without = FitBkgGeSimFromSeed(
      bkg_events, sig_events, bkg_peak_mus, sig_label, fit_label, seed,
      precomputed_sig_seed, nullptr, nullptr);
  if (!without.valid || without.bkg_channel.peaks.empty())
    return without;
  return FitBkgGeSimFromSeed(bkg_events, sig_events, bkg_peak_mus, sig_label,
                             fit_label, seed, precomputed_sig_seed, &escapes,
                             &without.bkg_channel.peaks.back());
}

// ------------------------------------------------------------- rate-sub fit

struct RateSubResult {
  Float_t mu = -1, mu_err = -1, chi2 = -1;
  Float_t sigma = -1, bkg0 = 0, slope = 0;
  Float_t tail_amp = -1, tail_ratio = -1; // -1 when no tail was fitted
  // Production's phase-averaged fit only: the spread over bin phases, the
  // shift to a quadratic residual background, and both with the calibration
  // transfer in quadrature.
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
inline RateSubResult FitRateSubHist(TH1F *res, Double_t lo, Double_t hi,
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
// 95.902 keV line, whose escapes it carries (none while capture_rate is 0).
struct RateSubShape {
  Double_t ge_tail_height_ratio = 0;
  Double_t ge_tau_over_sigma = 1;
  Double_t capture_rate = 0; // Gaussian yield, counts/s
  Double_t capture_mu = 0;
};

// Ge Gaussian + tail on a linear (optionally quadratic) background, plus the
// fixed escapes of the residual's 95.902 keV line (no shield: captures inside
// the crystal), with Ge's width and tail. Parameters:
//   [0] Ge height  [1] Ge mu  [2] Ge sigma  [3] Ge tail/gaus height ratio
//   [4] Ge tau/sigma  [5] bkg0  [6] bkg slope  [7] bkg curvature
//   [8] 95.902 keV Gaussian yield x bin width  [9] its centroid
inline Double_t RateSubShapedModel(Double_t *x, Double_t *p) {
  Double_t ge_gaus[3] = {p[1], p[2], p[0]};
  Double_t ge_tail[6] = {p[1], p[2], p[0] * p[3], p[4], 0.0, 0.0};
  Double_t u = x[0] - 68.75;
  Double_t escapes = 0;
  if (p[8] != 0)
    escapes =
        p[8] * EscapeDensity(x[0], p[9], p[2], p[3], p[4], ESCAPES_INSIDE);
  return FittingFunctions::Gaussian(x, ge_gaus) +
         FittingFunctions::LowTail(x, ge_tail) + p[5] + p[6] * u +
         p[7] * u * u + escapes;
}

// Fit the already-subtracted residual over the production window with the
// borrowed shape, seeded from the plain Gaussian fit. The residual
// background is linear unless quadratic is set.
inline RateSubResult FitRateSubHistShaped(TH1F *res, const RateSubShape &shape,
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
  const char *names[10] = {
      "GeHeight", "GeMu",     "GeSigma",      "GeTailRatio",      "GeTauRatio",
      "Bkg0",     "BkgSlope", "BkgCurvature", "CaptureYieldBins", "CaptureMu"};
  static Int_t model_counter = 0;
  TString mname = Form("ratesub_shaped_%s_%d", tag.Data(), ++model_counter);
  TF1 *model = new TF1(mname, RateSubShapedModel, RATESUB_LO, RATESUB_HI, 10);
  for (Int_t i = 0; i < 10; i++)
    model->SetParName(i, names[i]);
  model->SetParameter(7, 0.0);
  if (!quadratic)
    model->FixParameter(7, 0.0);
  model->FixParameter(8, shape.capture_rate * res->GetBinWidth(1));
  model->FixParameter(9, shape.capture_mu);
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
  if (fr.Get() && fr->IsValid() && model->GetNDF() > 0) {
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

// The 95.902 keV line is 1.2 times stronger in the signal runs than live time
// predicts, so the residual keeps a share of it and of its escapes. Its
// Gaussian yield and centroid in the residual of rate histograms over 88-104
// keV at the rate-sub fit's starting points, with the Ge tail, set shape's
// capture line; without one the fit carries no escapes.
inline void SetResidualCaptureLine(const std::vector<TH1F *> &sigs,
                                   const std::vector<TH1F *> &bkgs,
                                   const TString &tag, RateSubShape &shape) {
  if (sigs.size() != bkgs.size() || sigs.empty())
    return;
  std::vector<TH1F *> residuals;
  for (size_t k = 0; k < sigs.size(); k++) {
    TH1F *res = static_cast<TH1F *>(sigs[k]->Clone(
        TString::Format("capture_residual_%s_%d", tag.Data(), (Int_t)k)));
    res->SetDirectory(nullptr);
    res->Add(bkgs[k], -1.0);
    residuals.push_back(res);
  }
  CaptureLine line =
      AverageCaptureLine(residuals, shape.ge_tail_height_ratio,
                         shape.ge_tau_over_sigma, kFALSE, "ratesub_" + tag);
  if (line.valid) {
    shape.capture_rate = line.gaus_yield;
    shape.capture_mu = line.mu;
  } else {
    std::cerr << "WARNING: no residual 95.902 keV line for " << tag
              << "; the rate-sub fit carries no escapes" << std::endl;
  }
  for (size_t k = 0; k < residuals.size(); k++)
    delete residuals[k];
}

// The rate-sub fit of the residuals of pairs of rate histograms binned one comb
// period wide at N_RATESUB_PHASES starting points, averaged over the starting
// points. The linear residual background is the central model and the
// quadratic one its systematic. The residual is overdispersed (chi2/ndf
// 1.4-1.9 on 01/15-16: the pixels' discrete energy grids, which the two runs
// fill differently), so the fit error is scaled by sqrt(chi2/ndf). The first
// starting point's model and residual are handed back for drawing when the
// pointers are given (the caller then owns them).
inline RateSubResult FitRateSubPhases(const std::vector<TH1F *> &sigs,
                                      const std::vector<TH1F *> &bkgs,
                                      const RateSubShape &shape,
                                      const TString &tag,
                                      TF1 **out_model = nullptr,
                                      TH1F **out_residual = nullptr) {
  RateSubResult out;
  if (sigs.size() != bkgs.size() || sigs.empty())
    return out;
  Double_t sum_mu = 0, sum_mu2 = 0, sum_err = 0, sum_chi2 = 0, sum_quad = 0;
  Int_t n_linear = 0, n_quadratic = 0;
  RateSubResult fit;
  TF1 *model = nullptr;
  TH1F *res = nullptr;
  for (Int_t k = 0; k < (Int_t)sigs.size(); k++) {
    TString phase_tag = tag + Form("_phase%d", k);
    TH1F *phase_res =
        static_cast<TH1F *>(sigs[k]->Clone("rate_residual_" + phase_tag));
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
      sum_mu += linear.mu;
      sum_mu2 += linear.mu * linear.mu;
      sum_err += linear.mu_err;
      sum_chi2 += linear.chi2;
      n_linear++;
    }
    if (quadratic.valid) {
      sum_quad += quadratic.mu;
      n_quadratic++;
    }
    if (k == 0 && linear.valid && phase_model) {
      fit = linear;
      model = phase_model;
      res = phase_res;
    } else {
      delete phase_model;
      delete phase_res;
    }
  }
  if (n_linear == 0 || !model) {
    std::cerr << "WARNING: rate-sub chi2 fit failed for " << tag << std::endl;
    delete model;
    delete res;
    return out;
  }

  Double_t mean = sum_mu / n_linear;
  Double_t chi2 = sum_chi2 / n_linear;
  out = fit;
  out.mu = mean;
  out.chi2 = chi2;
  out.mu_err = (sum_err / n_linear) * TMath::Sqrt(TMath::Max(chi2, 1.0));
  out.phase_rms =
      TMath::Sqrt(TMath::Max(sum_mu2 / n_linear - mean * mean, 0.0));
  out.shape_shift =
      (n_quadratic > 0) ? TMath::Abs(sum_quad / n_quadratic - mean) : 0.0;
  out.sys_err = TMath::Sqrt(
      out.phase_rms * out.phase_rms + out.shape_shift * out.shape_shift +
      RATESUB_TRANSFER_SYS_KEV * RATESUB_TRANSFER_SYS_KEV);
  std::cout << "RATESUB-FIT " << tag << ": mu = " << std::fixed
            << std::setprecision(4) << out.mu << " +/- " << out.mu_err
            << " (fit, scaled by sqrt(chi2/ndf " << std::setprecision(2) << chi2
            << ")) +/- " << std::setprecision(4) << out.sys_err
            << " (sys: bin phase " << out.phase_rms << ", quadratic "
            << out.shape_shift << ", transfer " << RATESUB_TRANSFER_SYS_KEV
            << ") keV   sigma = " << out.sigma << std::endl;
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

#endif // GAINTRANSFER_FIT_HPP
