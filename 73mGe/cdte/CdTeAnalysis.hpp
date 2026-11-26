#ifndef CDTE_ANALYSIS_HPP
#define CDTE_ANALYSIS_HPP

// The 73mGe line in the PSI CdTe (Amptek, DT5780 channel 1, BOA beamline,
// October 2025), from our own trapezoid energies (ProcessWaveforms.cpp, through
// Analysis-Utilities' PulseHeightProcessor), calibrated as the CZT is: by the
// source runs' Am-241 and Ba-133 lines, carried to the beam runs through
// their own Pb K-alpha lines. What every cdte/ macro shares: the paths and
// runs, the processed trees and their rise-gain correction, the lineshape,
// the scale, the in-run fit and the source-run calibration.
//
// The in-run fit places every known line at its true energy on a free scale
// x = a + b (E - 70 + kappa (E - 70)^2), the relative curvature kappa fixed at
// the source runs' value (SourceCurvature.cpp): Pb Ka2/Ka1 (free ratio), the Pb
// Kb1' and Kb2' groups (sub-lines at fixed relative intensities, free group
// ratio), the 95.902 keV line's Cd and Te K X-ray escapes at the capture
// line's measured position minus the X-ray energies, and Ge, whose energy is
// a fit parameter. One lineshape for all lines (PeakShape, the source fit's
// as well), its width linear in energy, on a linear continuum whose slope is
// held at the line-free side bands'; the beam fits hold the shape's step at
// the source runs' value and leave its tails free. The fit gives the
// positions of Ge and the Pb K-alpha lines with their correlations; Ge's
// energy is then read on the source runs' line (GainTransfer).

#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"

#include <RtypesCore.h>
#include <TChain.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TString.h>

#include <iomanip>
#include <iostream>
#include <vector>

namespace CdTe {

// ------------------------------------------------------------ data and runs

// The raw CoMPASS files, read only.
const char *const RAW_DIR = "/labdata/PSI/CdTe/";
// ProcessWaveforms.cpp's trees (`features`), one per raw file, under the
// same name.
const char *const FEATURES_DIR = "/labdata/ANSG/PSI/processed/CdTe/features/";
const char *const FILE_PREFIX = "Data_CH1@DT5780P_14068_";

const Double_t SAMPLE_NS = 10.0;
// Spectra are filled on the Am-only scale, 59.5409 keV at 5577 trapezoid
// units; the in-run fit's own scale maps it to keV.
const Double_t KEV_PER_UNIT = 59.5409 / 5577.0;
const Double_t BIN = 0.05; // keV

inline TString FeaturesFile(const TString &run) {
  return TString(FEATURES_DIR) + FILE_PREFIX + run + ".root";
}

// The Ge beam runs: the CdTe behind the Ge target, on the days with source
// runs. 10/09's BeamTurnOn_4 has none, so no calibration can be carried to
// it; AllGeRuns adds it, for the comparisons of whole data sets.
inline std::vector<TString> GeRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Beam_GeTarget_20251010_1525_1");
  runs.push_back("PSI_Beam_GeTarget_20251010_1525_2");
  runs.push_back("PSI_Beam_Ge_20251011_1715_1");
  runs.push_back("PSI_Beam_Ge_20251011_1715_2");
  return runs;
}

inline std::vector<TString> AllGeRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Running_BeamTurnOn_20251009_1914_4");
  std::vector<TString> with_sources = GeRuns();
  runs.insert(runs.end(), with_sources.begin(), with_sources.end());
  return runs;
}

// Beam runs without Ge in front of the CdTe, with Pb: 10/10 NOGe, and
// 10/11-12 with the Ge target in front of the Si only.
inline std::vector<TString> NoGePbRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Beam_NOGe_20251010_1106_2");
  runs.push_back("PSI_Beam_NOGe_20251010_1106_3");
  runs.push_back("PSI_Beam_CdTe_NoGe_Si_YesGe_20251011_1056_1");
  runs.push_back("PSI_Beam_LongWf_CdTe_NoGe_Si_YesGe_20251012_1517_1");
  return runs;
}

// Am-241 and Ba-133 source runs, every day.
inline std::vector<TString> AmSourceRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_20251010_1422_1");
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_20251010_1422_2");
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_20251011_1617_1");
  runs.push_back("PSI_Calib_CdTe_Am_Si_Am_20251011_0920_1");
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_NoGe_20251010_1014_1");
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_NoGe_20251010_1014_2");
  runs.push_back("PSI_Calib_Long_CdTe_Am_Si_Sm_20251012_0917_1");
  runs.push_back("PSI_Calibration_CdTe_Am_Si_Am_20251013_0652_1");
  return runs;
}

inline std::vector<TString> BaSourceRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Calibration_CdTe_Ba_Si_Co_20251010_1433_1");
  runs.push_back("PSI_Calibration_CdTe_Ba_Si_Co_20251011_1617_1");
  runs.push_back("PSI_Calibration_CdTe_Ba_Si_Co_NoGe_20251010_1043_1");
  runs.push_back("PSI_Calibration_CdTe_Ba_Si_Co_NoGe_20251010_1043_2");
  runs.push_back("PSI_Calib_CdTe_ba_Si_am_20251011_1043_1");
  runs.push_back("PSI_Calib_CdTe_ba_Si_am_20251011_1043_2");
  runs.push_back("PSI_Calib_Long_CdTe_Ba_Si_Co_20251012_0917_1");
  runs.push_back("PSI_Calibration_CdTe_Ba_Si_Co_20251013_0652_1");
  return runs;
}

// The 10-90% rise-time slices, in ns. Slow events lose charge to hole
// trapping, so each slice is calibrated on its own lines. Only rise times
// below 281 ns are used: above about 282 ns every line grows a heavy low
// tail, however narrow the slice (SliceWidthTails.cpp), because there the
// charge lost at a given rise varies. Below it, each slice holds the same
// number of the Ge runs' events over 62-89 keV (heights over the rise gain;
// RiseTime.cpp prints the edges from EqualCountEdges).
const Int_t N_RISE_SLICES = 4;
const Double_t RESULT_MAX_RISE_NS = 281.0;
const Double_t RISE_SLICE_LO[N_RISE_SLICES] = {0, 269.24, 277.30, 279.44};
const Double_t RISE_SLICE_HI[N_RISE_SLICES] = {269.24, 277.30, 279.44, 281.0};

// ---------------------------------------------------------------- rise gain

// The pulse height changes with the 10-90% rise inside the slices: within
// 0.1% up to 295 ns, then down 1.3% by 330 ns. The lines do not share one
// rise distribution (it follows the depth they convert at), so the lines of
// one slice would sample different gains, which a fit reads as a curvature of
// the scale that changes from slice to slice (SourceCurvature.cpp). Every
// height is therefore divided by the gain at its rise: the Ba-133 81.00 keV
// core's position over 81.00 keV in the Ba source runs, in rise bins as
// narrow as the counts allow (each measures the core to 20 eV, from 0.25 ns
// where the rise distribution peaks to 10 ns in its tails), linear between
// the bins and held beyond the measured ones. The slices then only group
// events of similar width.
Bool_t g_rise_gain_correction = kTRUE;
TGraph *g_rise_gain = nullptr;

inline Double_t RiseGain(Double_t rise_ns);

// Energy (on the Am-only scale) against 10-90% rise, in 0.25 ns and BIN-wide
// bins, the heights divided by the rise gain or as measured.
inline TH2D *FillRiseHeight(const std::vector<TString> &runs,
                            const TString &name, Bool_t corrected) {
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  TH2D *h = new TH2D(name, ";10-90% rise [ns];Energy [keV]", 800, 200, 400,
                     (Int_t)(50.0 / BIN + 0.5), 50, 100);
  h->SetDirectory(nullptr);
  Float_t height = 0, rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("trapezoid_height", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("rise_time", &rise);
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    if (height <= 0 || rise <= 0)
      continue;
    Double_t rise_ns = SAMPLE_NS * rise;
    Double_t gain = corrected ? RiseGain(rise_ns) : 1.0;
    h->Fill(rise_ns, height * KEV_PER_UNIT / gain);
  }
  return h;
}

// The core of the line at true energy e in rise bins [first, last] of h: a
// Gaussian over +-0.3 keV around the maximum. False when the bins hold fewer
// than 150 counts within 0.5 keV of it or the core is not measured to
// max_err keV.
inline Bool_t CorePosition(TH2D *h, Int_t first, Int_t last, Double_t e,
                           Double_t max_err, Double_t &mu, Double_t &err) {
  TH1D *p = h->ProjectionY("core_projection", first, last);
  p->GetXaxis()->SetRangeUser(0.975 * e, 1.004 * e);
  Double_t mode = p->GetBinCenter(p->GetMaximumBin());
  Double_t height = p->GetMaximum();
  p->GetXaxis()->SetRange(0, 0);
  Double_t n = p->Integral(p->FindBin(mode - 0.5), p->FindBin(mode + 0.5));
  Bool_t measured = kFALSE;
  if (n >= 150) {
    TF1 core("core", "gaus", mode - 0.3, mode + 0.3);
    core.SetParameters(height, mode, 0.3);
    p->Fit(&core, "QLRN0");
    mu = core.GetParameter(1);
    err = core.GetParError(1);
    measured = err > 0 && err < max_err;
  }
  delete p;
  return measured;
}

// The core position over true energy e against rise, in bins of equal
// precision: h's rise bins merged until the core is measured to target_err
// keV, or, where 10 ns of rise do not reach it, to fallback_err; each point
// at its core events' mean rise, its x error spanning its merged bins.
inline TGraphAsymmErrors *CoreRatioByRise(TH2D *h, Double_t e,
                                          Double_t target_err,
                                          Double_t fallback_err) {
  const Double_t max_width_ns = 10.0;
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  TAxis *rise_axis = h->GetXaxis();
  const Int_t n_bins = h->GetNbinsX();
  Int_t first = 1;
  while (first <= n_bins) {
    Int_t last = first;
    Double_t mu = 0, err = 0;
    Bool_t measured = CorePosition(h, first, last, e, target_err, mu, err);
    while (!measured && last < n_bins &&
           rise_axis->GetBinUpEdge(last + 1) -
                   rise_axis->GetBinLowEdge(first) <=
               max_width_ns) {
      last++;
      measured = CorePosition(h, first, last, e, target_err, mu, err);
    }
    if (!measured)
      measured = CorePosition(h, first, last, e, fallback_err, mu, err);
    if (measured) {
      Int_t y_lo = h->GetYaxis()->FindBin(mu - 0.5);
      Int_t y_hi = h->GetYaxis()->FindBin(mu + 0.5);
      Double_t sw = 0, swx = 0;
      for (Int_t b = first; b <= last; b++) {
        Double_t w = h->Integral(b, b, y_lo, y_hi);
        sw += w;
        swx += w * h->GetXaxis()->GetBinCenter(b);
      }
      Double_t mean_rise = swx / sw;
      Int_t k = g->GetN();
      g->SetPoint(k, mean_rise, mu / e);
      g->SetPointError(k, mean_rise - rise_axis->GetBinLowEdge(first),
                       rise_axis->GetBinUpEdge(last) - mean_rise, err / e,
                       err / e);
    }
    first = last + 1;
  }
  return g;
}

inline TGraph *MeasureRiseGain() {
  TH2D *h = FillRiseHeight(BaSourceRuns(), "rise_gain_ba", kFALSE);
  TGraph *g = CoreRatioByRise(h, ReferenceLines::BA133_81, 0.02, 0.04);
  delete h;
  std::cout << "RISE GAIN from the Ba 81.00 keV core: " << g->GetN()
            << " bins, " << std::fixed << std::setprecision(1) << g->GetX()[0]
            << "-" << g->GetX()[g->GetN() - 1] << " ns:";
  for (Int_t k = 0; k < g->GetN(); k++)
    std::cout << " " << std::setprecision(1) << g->GetX()[k] << ":"
              << std::setprecision(5) << g->GetY()[k];
  std::cout << std::endl;
  return g;
}

inline Double_t RiseGain(Double_t rise_ns) {
  if (!g_rise_gain)
    g_rise_gain = MeasureRiseGain();
  const Int_t n = g_rise_gain->GetN();
  if (rise_ns <= g_rise_gain->GetX()[0])
    return g_rise_gain->GetY()[0];
  if (rise_ns >= g_rise_gain->GetX()[n - 1])
    return g_rise_gain->GetY()[n - 1];
  return g_rise_gain->Eval(rise_ns);
}

// The trapezoid spectrum of the given runs on the Am-only scale, each height
// divided by the rise gain (unless g_rise_gain_correction is off), events
// with a 10-90% rise in [rise_lo_ns, rise_hi_ns).
inline TH1D *FillRise(const std::vector<TString> &runs, const TString &name,
                      Double_t rise_lo_ns, Double_t rise_hi_ns) {
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  TH1D *h = new TH1D(name, Form(";Energy [keV];Counts / %.0f eV", 1000.0 * BIN),
                     (Int_t)(120.0 / BIN + 0.5), 0, 120);
  h->SetDirectory(nullptr);
  Float_t height = 0, rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("trapezoid_height", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("rise_time", &rise);
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    Double_t rise_ns = SAMPLE_NS * rise;
    if (height <= 0 || rise <= 0 || rise_ns < rise_lo_ns ||
        rise_ns >= rise_hi_ns)
      continue;
    Double_t gain = g_rise_gain_correction ? RiseGain(rise_ns) : 1.0;
    h->Fill(height * KEV_PER_UNIT / gain);
  }
  return h;
}

// The 10-90% rise distribution of the given runs' events with an energy (over
// the rise gain) in [e_lo, e_hi), in 0.01 ns bins up to 450 ns.
inline TH1D *FillRiseDistribution(const std::vector<TString> &runs,
                                  const TString &name, Double_t e_lo,
                                  Double_t e_hi) {
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  TH1D *h =
      new TH1D(name, ";Rise time 10-90% [ns];Counts / 0.01 ns", 45000, 0, 450);
  h->SetDirectory(nullptr);
  Float_t height = 0, rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("trapezoid_height", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("rise_time", &rise);
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    if (height <= 0 || rise <= 0)
      continue;
    Double_t rise_ns = SAMPLE_NS * rise;
    Double_t e = height * KEV_PER_UNIT / RiseGain(rise_ns);
    if (e >= e_lo && e < e_hi)
      h->Fill(rise_ns);
  }
  return h;
}

// The edges of n slices between 0 and max_ns holding equal counts of h: n + 1
// edges, the first 0 and the last max_ns.
inline std::vector<Double_t> EqualCountEdges(TH1D *h, Int_t n,
                                             Double_t max_ns) {
  const Int_t last_bin = h->FindBin(max_ns - 1e-6);
  const Double_t total = h->Integral(1, last_bin);
  std::vector<Double_t> edges;
  edges.push_back(0.0);
  Double_t cumulative = 0;
  for (Int_t b = 1; b <= last_bin && (Int_t)edges.size() < n; b++) {
    cumulative += h->GetBinContent(b);
    if (cumulative >= edges.size() * total / n)
      edges.push_back(h->GetXaxis()->GetBinUpEdge(b));
  }
  edges.push_back(max_ns);
  return edges;
}

// ------------------------------------------------------- lineshape and scale

// Every line in the CdTe fits, source and beam alike, has Analysis-Utilities'
// photopeak shape (FittingFunctions::PeakFunction): a Gaussian core with a
// step, a low exponential tail and a high exponential tail, each a ratio of
// the core's height, the decays in units of sigma. Its linear tail is left
// off: at zero slope it is the step, and its slope is undetermined when the
// tail is small. A line's position is its core's. shape[] holds step, low
// tail, its decay, high tail, its decay; the line is scaled so its core holds
// one count, in counts per BIN-wide bin.
const Int_t N_SHAPE = 5;

inline Double_t PeakShape(Double_t x, Double_t mu, Double_t sigma,
                          const Double_t *shape) {
  Double_t height = BIN / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma);
  Double_t par[12] = {mu, sigma, height,   shape[0], shape[1], shape[2],
                      0,  0,     shape[3], shape[4], 0,        0};
  return FittingFunctions::PeakFunction(&x, par);
}

// The shape's start in a staged fit, the step at parameter first and the
// rest after it in PeakShape's order: step and low tail free, high tail off.
inline void StartShape(TF1 *f, Int_t first) {
  f->SetParameter(first, 0.01);
  f->SetParLimits(first, 0, 0.5);
  f->SetParameter(first + 1, 0.1);
  f->SetParLimits(first + 1, 0, 2.0);
  f->SetParameter(first + 2, 2.0);
  f->SetParLimits(first + 2, 0.2, 20.0);
  f->FixParameter(first + 3, 0);
  f->FixParameter(first + 4, 1.0);
}

// The step's height over the core's in the source runs (SourceCurvature.cpp:
// 0.006-0.009 in the slices below 281 ns, 0.008 on average). A beam spectrum
// cannot tell a flat shelf under its lines from its continuum, so a free step
// there takes the continuum's counts; the beam fits hold it at this value
// unless g_free_step frees it (a check).
const Double_t STEP_SOURCE = 0.008;
Bool_t g_free_step = kFALSE;

inline void HoldStep(TF1 *f, Int_t first) {
  f->FixParameter(first, STEP_SOURCE);
}

// Where the beam fits take the shape's tails from. Below 281 ns the beam
// fits' own tails converge (low tail 0.1-0.7 of the core's height, decay
// 1.4-2.9 sigma), so they are free. Two alternatives are checks: the
// source runs' low tail per rise slice (with the rise-gain correction;
// SourceCurvature.cpp), and both tails of the beam spectrum's own 95.902
// keV capture line, fitted with them free (FitCaptureLine). Held in units of
// sigma, a tail keeps its share of a line's counts at any width.
enum { TAILS_FREE, TAILS_SOURCE, TAILS_CAPTURE };
Int_t g_tails = TAILS_FREE;
const Double_t LOW_EXP_SOURCE[N_RISE_SLICES] = {0.038, 0.017, 0.021, 0.025};
const Double_t LOW_EXP_DECAY_SOURCE[N_RISE_SLICES] = {5.036, 4.632, 5.003,
                                                      5.060};
// The source runs' high tail per slice, which the source-line fits hold
// with the step and the low tail.
const Double_t HIGH_EXP_SOURCE[N_RISE_SLICES] = {0.292, 0.078, 0.043, 0.020};
const Double_t HIGH_EXP_DECAY_SOURCE[N_RISE_SLICES] = {0.631, 0.546, 0.511,
                                                       0.452};
Double_t g_source_low_exp = LOW_EXP_SOURCE[0];
Double_t g_source_low_exp_decay = LOW_EXP_DECAY_SOURCE[0];
// The last capture-line fit's tails: low amplitude, its decay, high
// amplitude, its decay.
Double_t g_capture_tail[4] = {0, 1.0, 0, 1.0};

// The source low tail for a fit spanning the slices first to last: their
// mean.
inline void SetSourceLowTail(Int_t first, Int_t last) {
  g_source_low_exp = 0;
  g_source_low_exp_decay = 0;
  for (Int_t s = first; s <= last; s++) {
    g_source_low_exp += LOW_EXP_SOURCE[s];
    g_source_low_exp_decay += LOW_EXP_DECAY_SOURCE[s];
  }
  g_source_low_exp /= (last - first + 1);
  g_source_low_exp_decay /= (last - first + 1);
}

// For a fit over rise times 0 to rise_hi_ns.
inline void SetSourceLowTailBelow(Double_t rise_hi_ns) {
  Int_t last = 0;
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    if (RISE_SLICE_HI[s] <= rise_hi_ns)
      last = s;
  SetSourceLowTail(0, last);
}

inline void HoldTails(TF1 *f, Int_t first) {
  if (g_tails == TAILS_SOURCE) {
    f->FixParameter(first + 1, g_source_low_exp);
    f->FixParameter(first + 2, g_source_low_exp_decay);
  }
  if (g_tails == TAILS_CAPTURE)
    for (Int_t k = 0; k < 4; k++)
      f->FixParameter(first + 1 + k, g_capture_tail[k]);
}

// Whether the fits free the shape's high tail; off, it stays at zero (a
// check).
Bool_t g_high_tail = kTRUE;

inline void ReleaseHighTail(TF1 *f, Int_t first) {
  if (!g_high_tail)
    return;
  f->ReleaseParameter(first + 3);
  f->SetParameter(first + 3, 0.01);
  f->SetParLimits(first + 3, 0, 1.0);
  f->ReleaseParameter(first + 4);
  f->SetParLimits(first + 4, 0.2, 10.0);
}

// The high tail of a beam fit, unless it comes from the capture line.
inline void ReleaseBeamHighTail(TF1 *f, Int_t first) {
  if (g_tails != TAILS_CAPTURE)
    ReleaseHighTail(f, first);
}

// The source runs' relative curvature of the scale (SourceCurvature.cpp, the
// slices below 281 ns, heights corrected for the rise gain; errors from the
// likelihood profile, chi2 3.1/3 between the slices), 1/keV.
const Double_t KAPPA_SOURCE = 0.0534e-4, KAPPA_SOURCE_ERR = 0.92e-4;

// x = a + b (E - 70 + kappa (E - 70)^2). The offset a is the position of 70
// keV, among the lines, rather than of 0 keV, 70 keV below them, where it is
// almost fully correlated with the gain b.
inline Double_t CurvedScale(Double_t e, Double_t a, Double_t b,
                            Double_t kappa) {
  Double_t d = e - 70.0;
  return a + b * (d + kappa * d * d);
}

// The true energy at spectrum position x.
inline Double_t CurvedScaleInverse(Double_t x, Double_t a, Double_t b,
                                   Double_t kappa) {
  Double_t e = 70.0 + (x - a) / b;
  for (Int_t iteration = 0; iteration < 6; iteration++) {
    Double_t d = e - 70.0;
    e -= (CurvedScale(e, a, b, kappa) - x) / (b * (1.0 + 2.0 * kappa * d));
  }
  return e;
}

// ------------------------------------------------------------------ drawing

// A model and its components over the data, with pulls underneath; the y
// range from the window's contents.
inline void Plot(TH1D *h, TF1 *model, const std::vector<TGraph *> &components,
                 Double_t lo, Double_t hi, const TString &name,
                 const TString &label) {
  const Int_t n_pts = 800;
  TGraph *total = new TGraph(n_pts);
  for (Int_t i = 0; i < n_pts; i++) {
    Double_t x = lo + (hi - lo) * i / (n_pts - 1.0);
    total->SetPoint(i, x, model->Eval(x));
  }
  total->SetLineColor(kRed + 1);
  total->SetLineWidth(3);
  Double_t ymax = 0;
  for (Int_t b = h->FindBin(lo); b <= h->FindBin(hi); b++)
    ymax = TMath::Max(ymax, h->GetBinContent(b) + h->GetBinError(b));
  h->SetMinimum(0);
  h->SetMaximum(1.12 * ymax);
  PlottingUtils::PlotFitWithResiduals(h, total, components, lo, hi, name, "",
                                      label, kFALSE);
  h->GetXaxis()->SetRange(0, 0);
  h->SetMinimum();
  h->SetMaximum();
}

inline TGraph *Component(Double_t lo, Double_t hi, Int_t color, Int_t style,
                         Double_t (*eval)(Double_t, const Double_t *),
                         const Double_t *p) {
  const Int_t n_pts = 800;
  TGraph *g = new TGraph(n_pts);
  for (Int_t i = 0; i < n_pts; i++) {
    Double_t x = lo + (hi - lo) * i / (n_pts - 1.0);
    g->SetPoint(i, x, eval(x, p));
  }
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(2);
  return g;
}

// ------------------------------------------------------------- single lines

// One line in PeakShape, its shape free (the step held at the source runs'
// in a beam spectrum, and the low tail too with TAILS_SOURCE), on a linear
// continuum, on the spectrum's own scale: the 95.902 keV capture line (where
// its escapes sit, how strong it is, and the beam lines' tails) and the
// Am-241 check line of the source runs.
struct PeakFit {
  Double_t mu = 0, mu_err = 0, area = 0, area_err = 0, sigma = 0, chi2 = 0;
  Double_t tail[4] = {0, 1.0, 0, 1.0};
  Double_t tail_err[4] = {0, 0, 0, 0};
};

// [0] continuum at mu [1] its slope [2] sigma [3] the core's counts [4] mu
// [5-9] the shape
inline Double_t PeakModel(Double_t *x, Double_t *p) {
  return p[0] + p[1] * (x[0] - p[4]) +
         p[3] * PeakShape(x[0], p[4], p[2], &p[5]);
}

// Staged: the high tail off, then free.
inline PeakFit FitPeak(TH1D *h, Double_t e_line, Double_t lo, Double_t hi,
                       const TString &tag, Bool_t beam) {
  h->GetXaxis()->SetRangeUser(e_line - 1.0, e_line + 1.0);
  Double_t top = h->GetBinCenter(h->GetMaximumBin());
  Double_t height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  TF1 *f = new TF1("peak_" + tag, PeakModel, lo, hi, 5 + N_SHAPE);
  f->SetNpx(2000);
  const Double_t sigma0 = 0.3;
  f->SetParameter(0, h->GetBinContent(h->FindBin(hi - 0.1)));
  f->SetParameter(1, 0);
  f->SetParameter(2, sigma0);
  f->SetParLimits(2, 0.05, 1.5);
  f->SetParameter(3, height * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 / BIN);
  f->SetParLimits(3, 0, 1e10);
  f->SetParameter(4, top);
  f->SetParLimits(4, top - 1.0, top + 1.0);
  StartShape(f, 5);
  if (beam)
    HoldStep(f, 5);
  if (beam && g_tails == TAILS_SOURCE)
    HoldTails(f, 5);
  h->Fit(f, "QLRN0");
  ReleaseHighTail(f, 5);
  h->Fit(f, "QLRN0");
  h->Fit(f, "QLRN0M");
  PeakFit out;
  out.mu = f->GetParameter(4);
  out.mu_err = f->GetParError(4);
  out.area = f->GetParameter(3);
  out.area_err = f->GetParError(3);
  out.sigma = f->GetParameter(2);
  out.chi2 = f->GetChisquare() / f->GetNDF();
  for (Int_t k = 0; k < 4; k++) {
    out.tail[k] = f->GetParameter(6 + k);
    out.tail_err[k] = f->GetParError(6 + k);
  }
  delete f;
  return out;
}

// Also keeps its tails for the beam fits (TAILS_CAPTURE).
inline PeakFit FitCaptureLine(TH1D *h, const TString &tag) {
  PeakFit fit = FitPeak(h, ReferenceLines::CD113_CAPTURE, 90.0, 99.0,
                        "capture_" + tag, kTRUE);
  for (Int_t k = 0; k < 4; k++)
    g_capture_tail[k] = fit.tail[k];
  return fit;
}

// ------------------------------------------------------------ in-run model

enum {
  Q_A,
  Q_B,
  Q_KAPPA, // x = a + b (E - 70 + kappa (E - 70)^2)
  Q_SIG0,
  Q_SIG1, // sigma(E) = sig0 + sig1 (E - 70)
  Q_STEP, // the shape every line shares, in PeakShape's order
  Q_LOW_EXP,
  Q_LOW_EXP_DECAY,
  Q_HIGH_EXP,
  Q_HIGH_EXP_DECAY,
  Q_BKG,       // the continuum at the window's middle
  Q_BKG_SLOPE, // its slope per keV, held at the side bands' (StartContinuum)
  Q_GE_E,
  Q_GE_A,
  Q_KA1_A,
  Q_KA2_R,
  Q_KB1_A,
  Q_KB2_R,
  Q_ESC_A,  // the Cd K-beta escape's yield; the others tied to it
  Q_GE_SIG, // Ge's own width and shape, only with g_ge_own_shape
  Q_GE_STEP,
  Q_GE_LOW_EXP,
  Q_GE_LOW_EXP_DECAY,
  Q_GE_HIGH_EXP,
  Q_GE_HIGH_EXP_DECAY,
  N_Q
};

// The fit window.
Double_t g_lo = 62.0, g_hi = 89.0;
// Variants.
Bool_t g_ge_own_shape = kFALSE;
// Ge's own width with the shared shape.
Bool_t g_ge_own_width = kFALSE;
Bool_t g_free_kappa = kFALSE;
Double_t g_kappa = KAPPA_SOURCE;
// The last FitInRun's covariance matrix, every parameter (zero where fixed).
TMatrixDSym g_inrun_cov(N_Q);

// The Pb K-beta groups the fit uses: ReferenceLines' unless a variant
// changes the sub-lines' intensities.
ReferenceLines::GroupLine g_pb_kb1[ReferenceLines::N_PB_KB1];
ReferenceLines::GroupLine g_pb_kb2[ReferenceLines::N_PB_KB2];

inline void ResetKbGroups() {
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB1; k++)
    g_pb_kb1[k] = ReferenceLines::PB_KB1[k];
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB2; k++)
    g_pb_kb2[k] = ReferenceLines::PB_KB2[k];
}

// The escapes' weights per unit Cd K-beta yield, and the capture line's
// position on the spectrum's scale.
Double_t g_escape_weight[ReferenceLines::N_ESCAPE_LINES];
Double_t g_capture_meas = ReferenceLines::CD113_CAPTURE;

// Each escape line's share of its series, Geant4's for captures in CdTe
// (cdte96_surface, geant4/Escape96Plots.cpp's "shares in family"): Cd Ka1,
// Ka2, Kb1', Kb2', Te Ka1, Ka2, Kb1', Kb2'.
const Double_t GEANT4_SERIES_SHARE[ReferenceLines::N_ESCAPE_LINES] = {
    0.6568, 0.3432, 0.8380, 0.1620, 0.6576, 0.3424, 0.8244, 0.1756};

// The escape tie: Cd K-alpha from the Cd Kb/Ka ratio, Te K-alpha from its
// ratio to Cd K-beta, Te K-beta at 0.30 of Te K-alpha (Geant4, captures in
// CdTe; the runs do not resolve it).
inline void SetEscapeTie(Double_t cd_kb_over_ka, Double_t te_ka_over_cd_kb) {
  Double_t series_total[4] = {1.0 / cd_kb_over_ka, 1.0, te_ka_over_cd_kb,
                              0.30 * te_ka_over_cd_kb};
  for (Int_t i = 0; i < ReferenceLines::N_ESCAPE_LINES; i++)
    g_escape_weight[i] =
        series_total[ReferenceLines::ESCAPE_SERIES[i]] * GEANT4_SERIES_SHARE[i];
}

// The escape tie's values, for every macro that sets it. Te Ka over Cd Kb is
// measured by NoGeEscapes.cpp on the no-Ge Pb runs with Pb Ka2/Ka1 at DDEP's
// 0.595 (free, it goes to 0.30-0.54, trading with the Cd Ka escapes under Pb
// Ka2): rise < 281 ns, and per slice. Cd Kb/Ka is Geant4's 0.275 for captures
// in CdTe; the alternative is NoGeEscapes.cpp's own, rise < 281 ns. The
// continuum's slope is held at the side bands', the step at the source runs'
// and the tails are free. A slice whose own fit is unphysical keeps the
// rise < 281 ns value.
const Double_t TE_CDKB = 0.354, TE_CDKB_ERR = 0.040;
const Double_t TE_CDKB_SLICE[N_RISE_SLICES] = {0.247, 0.355, 0.395, 0.342};
const Double_t KB_KA_NOMINAL = 0.275, KB_KA_ALT = 0.298;

inline Double_t ScalePosition(Double_t e, const Double_t *p) {
  return CurvedScale(e, p[Q_A], p[Q_B], p[Q_KAPPA]);
}

inline Double_t ScaleEnergy(Double_t x, const Double_t *p) {
  return CurvedScaleInverse(x, p[Q_A], p[Q_B], p[Q_KAPPA]);
}

// A line at true energy e whose core holds one count, in counts per bin.
inline Double_t Shape(Double_t x, Double_t e, const Double_t *p) {
  Double_t sigma = p[Q_SIG0] + p[Q_SIG1] * (e - 70.0);
  return PeakShape(x, ScalePosition(e, p), sigma, &p[Q_STEP]);
}

inline Double_t Bkg(Double_t x, const Double_t *p) {
  return p[Q_BKG] + p[Q_BKG_SLOPE] * (x - 0.5 * (g_lo + g_hi));
}

inline Double_t PbLines(Double_t x, const Double_t *p) {
  Double_t v = p[Q_KA1_A] * (Shape(x, ReferenceLines::PB_KA1, p) +
                             p[Q_KA2_R] * Shape(x, ReferenceLines::PB_KA2, p));
  Double_t kb1 = 0, kb2 = 0;
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB1; k++)
    kb1 += g_pb_kb1[k].weight * Shape(x, g_pb_kb1[k].energy, p);
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB2; k++)
    kb2 += g_pb_kb2[k].weight * Shape(x, g_pb_kb2[k].energy, p);
  return v + p[Q_KB1_A] * (kb1 + p[Q_KB2_R] * kb2);
}

inline Double_t Escapes(Double_t x, const Double_t *p) {
  Double_t e_capture = ScaleEnergy(g_capture_meas, p);
  Double_t v = 0;
  for (Int_t k = 0; k < ReferenceLines::N_ESCAPE_LINES; k++)
    v += g_escape_weight[k] *
         Shape(x, e_capture - ReferenceLines::ESCAPE_XRAY_KEV[k], p);
  return p[Q_ESC_A] * v;
}

inline Double_t GePeak(Double_t x, const Double_t *p) {
  if (g_ge_own_width)
    return p[Q_GE_A] *
           PeakShape(x, ScalePosition(p[Q_GE_E], p), p[Q_GE_SIG], &p[Q_STEP]);
  if (!g_ge_own_shape)
    return p[Q_GE_A] * Shape(x, p[Q_GE_E], p);
  return p[Q_GE_A] *
         PeakShape(x, ScalePosition(p[Q_GE_E], p), p[Q_GE_SIG], &p[Q_GE_STEP]);
}

inline Double_t InRunModel(Double_t *x, Double_t *p) {
  return Bkg(x[0], p) + PbLines(x[0], p) + Escapes(x[0], p) + GePeak(x[0], p);
}

// Components for drawing, each on top of the continuum.
inline Double_t DrawBkg(Double_t x, const Double_t *p) { return Bkg(x, p); }
inline Double_t DrawPb(Double_t x, const Double_t *p) {
  return Bkg(x, p) + PbLines(x, p);
}
inline Double_t DrawEsc(Double_t x, const Double_t *p) {
  return Bkg(x, p) + Escapes(x, p);
}
inline Double_t DrawGe(Double_t x, const Double_t *p) {
  return Bkg(x, p) + GePeak(x, p);
}

// The continuum's slope comes from the spectrum's line-free side bands, not
// from the fit: left free, it trades counts with the lines' step and low
// tail, which add counts only below each line, and turns steep. A straight
// line is fitted to the side bands; the fit holds its slope and leaves the
// level free. The low band sits between the Pb K-beta groups' Cd K-alpha
// escapes, at 61.3-62.0 and 64.1-64.2 keV.
const Double_t SIDEBAND_LOW[2] = {62.8, 63.6};
const Double_t SIDEBAND_HIGH[2] = {78.5, 82.0};

// The slope variants of the background systematic: the side bands' slope
// shifted by g_slope_shift of its errors, or a flat continuum.
Double_t g_slope_shift = 0;
Bool_t g_flat_continuum = kFALSE;
// The side bands' slope and its error in the last StartContinuum.
Double_t g_sideband_slope = 0, g_sideband_slope_err = 0;

inline Bool_t InSideband(Double_t x) {
  return (x >= SIDEBAND_LOW[0] && x <= SIDEBAND_LOW[1]) ||
         (x >= SIDEBAND_HIGH[0] && x <= SIDEBAND_HIGH[1]);
}

// The side bands' straight line by weighted least squares (Poisson errors):
// its value at x0, and its slope with the slope's error.
inline void FitSidebands(TH1D *h, Double_t x0, Double_t &level, Double_t &slope,
                         Double_t &slope_err) {
  Double_t s = 0, sx = 0, sxx = 0, sy = 0, sxy = 0;
  for (Int_t b = 1; b <= h->GetNbinsX(); b++) {
    Double_t x = h->GetBinCenter(b);
    if (!InSideband(x))
      continue;
    Double_t y = h->GetBinContent(b);
    Double_t w = 1.0 / TMath::Max(y, 1.0);
    Double_t dx = x - x0;
    s += w;
    sx += w * dx;
    sxx += w * dx * dx;
    sy += w * y;
    sxy += w * dx * y;
  }
  Double_t det = s * sxx - sx * sx;
  slope = (s * sxy - sx * sy) / det;
  slope_err = TMath::Sqrt(s / det);
  level = (sy - slope * sx) / s;
}

// The slope held at the side bands' (or a variant's), the level started at
// the side bands' line and kept high enough that the continuum stays at zero
// or above across the window.
inline void StartContinuum(TF1 *f, TH1D *h) {
  const Double_t middle = 0.5 * (g_lo + g_hi);
  Double_t level = 0;
  FitSidebands(h, middle, level, g_sideband_slope, g_sideband_slope_err);
  Double_t slope = g_flat_continuum ? 0.0
                                    : g_sideband_slope +
                                          g_slope_shift * g_sideband_slope_err;
  f->FixParameter(Q_BKG_SLOPE, slope);
  Double_t floor = TMath::Abs(slope) * 0.5 * (g_hi - g_lo);
  f->SetParameter(Q_BKG, TMath::Max(level, floor + 1.0));
  f->SetParLimits(Q_BKG, floor, 10.0 * h->GetMaximum() + 10.0);
}

// The width, the shared shape (the step held, the tails held or the high
// tail off) and the curvature at their start, Ge's own shape held; for
// FitInRun and the no-Ge fit.
inline void StartInRunShape(TF1 *f, Double_t sigma0) {
  f->SetParameter(Q_SIG0, sigma0);
  f->SetParLimits(Q_SIG0, 0.1, 1.5);
  f->FixParameter(Q_SIG1, 0.0);
  StartShape(f, Q_STEP);
  if (!g_free_step)
    HoldStep(f, Q_STEP);
  HoldTails(f, Q_STEP);
  f->FixParameter(Q_KAPPA, g_kappa);
  f->FixParameter(Q_GE_SIG, sigma0);
  f->FixParameter(Q_GE_STEP, STEP_SOURCE);
  f->FixParameter(Q_GE_LOW_EXP, 0.1);
  f->FixParameter(Q_GE_LOW_EXP_DECAY, 2.0);
  f->FixParameter(Q_GE_HIGH_EXP, 0);
  f->FixParameter(Q_GE_HIGH_EXP_DECAY, 1.0);
}

// The width's slope, the high tail unless it is held, and the curvature or
// Ge's own width in those checks.
inline void ReleaseInRunShape(TF1 *f) {
  f->ReleaseParameter(Q_SIG1);
  f->SetParLimits(Q_SIG1, -0.05, 0.05);
  ReleaseBeamHighTail(f, Q_STEP);
  if (g_free_kappa) {
    f->ReleaseParameter(Q_KAPPA);
    f->SetParLimits(Q_KAPPA, -0.01, 0.01);
  }
  if (g_ge_own_width) {
    f->ReleaseParameter(Q_GE_SIG);
    f->SetParameter(Q_GE_SIG, f->GetParameter(Q_SIG0) +
                                  f->GetParameter(Q_SIG1) *
                                      (f->GetParameter(Q_GE_E) - 70.0));
    f->SetParLimits(Q_GE_SIG, 0.1, 1.5);
  }
  if (g_ge_own_shape) {
    f->ReleaseParameter(Q_GE_SIG);
    f->SetParameter(Q_GE_SIG, f->GetParameter(Q_SIG0) +
                                  f->GetParameter(Q_SIG1) *
                                      (f->GetParameter(Q_GE_E) - 70.0));
    f->SetParLimits(Q_GE_SIG, 0.1, 1.5);
    StartShape(f, Q_GE_STEP);
    HoldStep(f, Q_GE_STEP);
    for (Int_t k = 1; k < 3; k++)
      f->SetParameter(Q_GE_STEP + k, f->GetParameter(Q_STEP + k));
    HoldTails(f, Q_GE_STEP);
    ReleaseBeamHighTail(f, Q_GE_STEP);
  }
}

// Staged: (1) the scale fixed at the Pb Ka1 peak's offset, ratios fixed, no
// escapes, the high tail off; (2) the scale released; (3) the ratios,
// escapes, the width's slope and the high tail released (and Ge's own
// shape or the curvature in those checks). The caller sets the capture line
// and the escape tie first.
inline TF1 *FitInRun(TH1D *h, const TString &tag) {
  TF1 *f = new TF1("inrun_" + tag, InRunModel, g_lo, g_hi, N_Q);
  f->SetNpx(3000);
  h->GetXaxis()->SetRangeUser(74.2, 75.6);
  Double_t ka1_mode = h->GetBinCenter(h->GetMaximumBin());
  Double_t ka1_height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  Double_t level = h->GetBinContent(h->FindBin(65.0));
  Double_t a0 = ka1_mode - (ReferenceLines::PB_KA1 - 70.0);
  const Double_t sigma0 = 0.6;
  Double_t ka1 =
      (ka1_height - level) * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 / BIN;

  StartInRunShape(f, sigma0);
  StartContinuum(f, h);
  f->SetParameter(Q_GE_E, 68.75);
  f->SetParLimits(Q_GE_E, 68.0, 69.5);
  f->SetParameter(Q_GE_A, 0.4 * ka1);
  f->SetParLimits(Q_GE_A, 0, 1e9);
  f->SetParameter(Q_KA1_A, ka1);
  f->SetParLimits(Q_KA1_A, 0, 1e9);
  f->SetParameter(Q_KB1_A, 0.25 * ka1);
  f->SetParLimits(Q_KB1_A, 0, 1e9);

  f->FixParameter(Q_A, a0);
  f->FixParameter(Q_B, 1.0);
  f->FixParameter(Q_KA2_R, 0.6);
  f->FixParameter(Q_KB2_R, 0.3);
  f->FixParameter(Q_ESC_A, 0.0);
  h->Fit(f, "QLRN0");

  f->ReleaseParameter(Q_A);
  f->ReleaseParameter(Q_B);
  f->SetParLimits(Q_A, a0 - 1.0, a0 + 1.0);
  f->SetParLimits(Q_B, 0.97, 1.03);
  h->Fit(f, "QLRN0");

  f->ReleaseParameter(Q_KA2_R);
  f->SetParLimits(Q_KA2_R, 0.3, 1.2);
  f->ReleaseParameter(Q_KB2_R);
  f->SetParLimits(Q_KB2_R, 0.1, 0.6);
  f->ReleaseParameter(Q_ESC_A);
  f->SetParameter(Q_ESC_A, 0.03 * f->GetParameter(Q_KA1_A));
  f->SetParLimits(Q_ESC_A, 0, f->GetParameter(Q_KA1_A));
  ReleaseInRunShape(f);
  h->Fit(f, "QLRN0");
  TFitResultPtr result = h->Fit(f, "QLRN0MS");
  g_inrun_cov.Zero();
  if (result.Get() != nullptr && result->CovMatrixStatus() > 0) {
    g_inrun_cov = result->GetCovarianceMatrix();
  } else {
    for (Int_t k = 0; k < N_Q; k++)
      g_inrun_cov(k, k) = f->GetParError(k) * f->GetParError(k);
  }
  return f;
}

inline void DrawInRun(TH1D *h, TF1 *f, const TString &name,
                      const TString &label) {
  const Double_t *p = f->GetParameters();
  std::vector<TGraph *> components;
  components.push_back(Component(g_lo, g_hi, kGray + 2, 2, DrawBkg, p));
  components.push_back(Component(g_lo, g_hi, kOrange + 7, 1, DrawPb, p));
  components.push_back(Component(g_lo, g_hi, kBlue + 1, 1, DrawEsc, p));
  components.push_back(Component(g_lo, g_hi, kGreen + 2, 1, DrawGe, p));
  Plot(h, f, components, g_lo, g_hi, name, label);
}

// ---------------------------------------------------- source-run calibration

// Ge's energy is calibrated as the CZT's: in each rise slice and on each day,
// the source runs' Am-241 59.5409 keV and Ba-133 80.9979 keV lines fix the
// scale x = a + b (E - 70 + kappa (E - 70)^2), kappa the in-run fit's
// (g_kappa), which is carried to that day's Ge beam runs as a gain through
// their own Pb K-alpha lines, whose positions (with Ge's) come from the in-run
// fit. Am lies below Ge, so Ge is interpolated. The CdTe's bias was not changed
// between the runs, so the beam runs differ from the source runs by a gain
// (up to 1.1% in the fastest slice).
const Int_t N_GE_DAYS = 2;
const Int_t GE_DAYS[N_GE_DAYS] = {10, 11};
// Pb K-alpha2's weight in the transfer against K-alpha1's: their intensity.
const Double_t TRANSFER_KA2_WEIGHT = 0.595;
// Pb K-beta1's blend centroid (K-M3 with K-M2 at 0.52), a check the transfer
// does not use.
const Double_t PB_KB1_BLEND = (84.93908 + 0.5215 * 84.45045) / 1.5215;

// The runs of one day (the date in the run name).
inline std::vector<TString> OfDay(const std::vector<TString> &runs, Int_t day) {
  std::vector<TString> out;
  TString date = Form("_202510%02d_", day);
  for (size_t i = 0; i < runs.size(); i++)
    if (runs[i].Contains(date))
      out.push_back(runs[i]);
  return out;
}

// The source runs' scale x = a + b (E - 70 + kappa (E - 70)^2) through the Am
// and Ba positions.
struct SourceLine {
  Double_t x_am = 0, x_am_err = 0, x_ba = 0, x_ba_err = 0;
  Double_t B() const {
    return (x_ba - x_am) /
           (CurvedScale(ReferenceLines::BA133_81, 0, 1, g_kappa) -
            CurvedScale(ReferenceLines::AM241, 0, 1, g_kappa));
  }
  Double_t A() const {
    return x_am - B() * CurvedScale(ReferenceLines::AM241, 0, 1, g_kappa);
  }
  Double_t Energy(Double_t x) const {
    return CurvedScaleInverse(x, A(), B(), g_kappa);
  }
  Double_t Position(Double_t e) const {
    return CurvedScale(e, A(), B(), g_kappa);
  }
};

// A source line's position on one day in one slice: the line in PeakShape
// with the source runs' shape of that slice held (with its shape free, a
// day's spectrum in one slice cannot fix the tails and the fit rails), its
// width, yield and position free, on a linear continuum. The Ba-133 80.9979
// keV fit carries the 79.6142 keV line at ENSDF's 2.65 : 32.9, 1.38 keV
// below on the source runs' scale (b = 0.998, SourceCurvature.cpp).
const Double_t SOURCE_SCALE_B = 0.998;

// [0] continuum at mu [1] its slope [2] sigma [3] the core's counts [4] mu
// [5-9] the shape [10] the 79.61 keV line's yield over the core's
inline Double_t SourcePeakModel(Double_t *x, Double_t *p) {
  const Double_t spacing =
      SOURCE_SCALE_B * (ReferenceLines::BA133_81 - ReferenceLines::BA133_79);
  return p[0] + p[1] * (x[0] - p[4]) +
         p[3] * (PeakShape(x[0], p[4], p[2], &p[5]) +
                 p[10] * PeakShape(x[0], p[4] - spacing, p[2], &p[5]));
}

inline PeakFit FitSourcePeak(TH1D *h, Double_t e_line, Double_t lo, Double_t hi,
                             Int_t slice, Double_t ba79_ratio,
                             const TString &tag) {
  h->GetXaxis()->SetRangeUser(e_line - 1.0, e_line + 1.0);
  Double_t top = h->GetBinCenter(h->GetMaximumBin());
  Double_t height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  TF1 *f = new TF1("source_peak_" + tag, SourcePeakModel, lo, hi, 11);
  f->SetNpx(2000);
  const Double_t sigma0 = 0.3;
  f->SetParameter(0, h->GetBinContent(h->FindBin(hi - 0.1)));
  f->SetParameter(1, 0);
  f->SetParameter(2, sigma0);
  f->SetParLimits(2, 0.05, 1.5);
  f->SetParameter(3, height * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 / BIN);
  f->SetParLimits(3, 0, 1e10);
  f->SetParameter(4, top);
  f->SetParLimits(4, top - 1.0, top + 1.0);
  f->FixParameter(5, STEP_SOURCE);
  f->FixParameter(6, LOW_EXP_SOURCE[slice]);
  f->FixParameter(7, LOW_EXP_DECAY_SOURCE[slice]);
  f->FixParameter(8, HIGH_EXP_SOURCE[slice]);
  f->FixParameter(9, HIGH_EXP_DECAY_SOURCE[slice]);
  f->FixParameter(10, ba79_ratio);
  h->Fit(f, "QLRN0");
  h->Fit(f, "QLRN0M");
  PeakFit out;
  out.mu = f->GetParameter(4);
  out.mu_err = f->GetParError(4);
  out.area = f->GetParameter(3);
  out.area_err = f->GetParError(3);
  out.sigma = f->GetParameter(2);
  out.chi2 = f->GetChisquare() / f->GetNDF();
  delete f;
  return out;
}

inline SourceLine FitSourceLine(Int_t day, Int_t slice) {
  SourceLine line;
  const TString tag = Form("d%d_s%d", day, slice);
  TH1D *am = FillRise(OfDay(AmSourceRuns(), day), "am_" + tag,
                      RISE_SLICE_LO[slice], RISE_SLICE_HI[slice]);
  TH1D *ba = FillRise(OfDay(BaSourceRuns(), day), "ba_" + tag,
                      RISE_SLICE_LO[slice], RISE_SLICE_HI[slice]);
  PeakFit am_peak = FitSourcePeak(am, ReferenceLines::AM241, 56.0, 61.5, slice,
                                  0, "am_" + tag);
  PeakFit ba_peak = FitSourcePeak(ba, ReferenceLines::BA133_81, 78.9, 82.5,
                                  slice, 2.65 / 32.9, "ba_" + tag);
  line.x_am = am_peak.mu;
  line.x_am_err = am_peak.mu_err;
  line.x_ba = ba_peak.mu;
  line.x_ba_err = ba_peak.mu_err;
  std::cout << "SOURCE 10/" << day << ", rise " << std::fixed
            << std::setprecision(2) << RISE_SLICE_LO[slice] << "-"
            << RISE_SLICE_HI[slice] << " ns: Am at " << std::setprecision(4)
            << line.x_am << " +- " << std::setprecision(1)
            << 1000 * line.x_am_err << " eV (chi2/ndf " << std::setprecision(2)
            << am_peak.chi2 << "), Ba 81 at " << std::setprecision(4)
            << line.x_ba << " +- " << std::setprecision(1)
            << 1000 * line.x_ba_err << " eV (chi2/ndf " << std::setprecision(2)
            << ba_peak.chi2 << ") (spectrum keV)" << std::endl;
  delete am;
  delete ba;
  return line;
}

// Every Ge day's source line in every slice.
inline void FitSourceLines(SourceLine lines[N_GE_DAYS][N_RISE_SLICES]) {
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    for (Int_t s = 0; s < N_RISE_SLICES; s++)
      lines[d][s] = FitSourceLine(GE_DAYS[d], s);
}

// The beam-to-source gain from Pb K-alpha1 and K-alpha2 by least squares.
inline Double_t BeamGain(const SourceLine &line, const Double_t *p) {
  const Double_t e_pb[2] = {ReferenceLines::PB_KA1, ReferenceLines::PB_KA2};
  const Double_t w_pb[2] = {1.0, TRANSFER_KA2_WEIGHT};
  Double_t num = 0, den = 0;
  for (Int_t k = 0; k < 2; k++) {
    Double_t x_beam = ScalePosition(e_pb[k], p);
    Double_t x_source = line.Position(e_pb[k]);
    num += w_pb[k] * x_beam * x_source;
    den += w_pb[k] * x_source * x_source;
  }
  return num / den;
}

// The source scale's energy of the line the in-run fit puts at e_in_run.
inline Double_t TransferredEnergy(const SourceLine &line, const Double_t *p,
                                  Double_t e_in_run) {
  return line.Energy(ScalePosition(e_in_run, p) / BeamGain(line, p));
}

// Ge from an in-run fit (its parameters and the covariance g_inrun_cov) on a
// day's source line. err is the statistical error, from the in-run fit's
// scale and Ge position with their correlations; d_am and d_ba are Ge's
// derivatives with respect to the line's Am and Ba positions, the gain's
// response included; ka2_residual and kb1_residual are Pb K-alpha2 and
// K-beta1 on the transferred scale.
struct Transferred {
  Double_t e = 0, err = 0, d_am = 0, d_ba = 0, gain = 1;
  Double_t ka2_residual = 0, kb1_residual = 0;
};

inline Transferred GainTransfer(const SourceLine &line, TF1 *f) {
  std::vector<Double_t> p(f->GetParameters(), f->GetParameters() + N_Q);
  Transferred t;
  t.gain = BeamGain(line, &p[0]);
  t.e = TransferredEnergy(line, &p[0], p[Q_GE_E]);
  t.ka2_residual = TransferredEnergy(line, &p[0], ReferenceLines::PB_KA2) -
                   ReferenceLines::PB_KA2;
  t.kb1_residual = TransferredEnergy(line, &p[0], PB_KB1_BLEND) - PB_KB1_BLEND;

  // Numerical derivatives over one standard deviation; Ge's energy is close
  // to linear in all of these.
  const Int_t N_MOVED = 4;
  const Int_t moved[N_MOVED] = {Q_A, Q_B, Q_KAPPA, Q_GE_E};
  Double_t gradient[N_MOVED];
  for (Int_t i = 0; i < N_MOVED; i++) {
    gradient[i] = 0;
    Double_t step = TMath::Sqrt(g_inrun_cov(moved[i], moved[i]));
    if (!(step > 0))
      continue;
    const Double_t start = p[moved[i]];
    p[moved[i]] = start + step;
    Double_t e_up = TransferredEnergy(line, &p[0], p[Q_GE_E]);
    p[moved[i]] = start - step;
    Double_t e_down = TransferredEnergy(line, &p[0], p[Q_GE_E]);
    p[moved[i]] = start;
    gradient[i] = (e_up - e_down) / (2.0 * step);
  }
  Double_t variance = 0;
  for (Int_t i = 0; i < N_MOVED; i++)
    for (Int_t j = 0; j < N_MOVED; j++)
      variance += gradient[i] * gradient[j] * g_inrun_cov(moved[i], moved[j]);
  t.err = TMath::Sqrt(TMath::Max(variance, 0.0));

  SourceLine shifted = line;
  shifted.x_am = line.x_am + line.x_am_err;
  Double_t e_up = TransferredEnergy(shifted, &p[0], p[Q_GE_E]);
  shifted.x_am = line.x_am - line.x_am_err;
  Double_t e_down = TransferredEnergy(shifted, &p[0], p[Q_GE_E]);
  t.d_am = (e_up - e_down) / (2.0 * line.x_am_err);
  shifted = line;
  shifted.x_ba = line.x_ba + line.x_ba_err;
  e_up = TransferredEnergy(shifted, &p[0], p[Q_GE_E]);
  shifted.x_ba = line.x_ba - line.x_ba_err;
  e_down = TransferredEnergy(shifted, &p[0], p[Q_GE_E]);
  t.d_ba = (e_up - e_down) / (2.0 * line.x_ba_err);
  return t;
}

} // namespace CdTe

#endif
