#ifndef CDTE_ANALYSIS_HPP
#define CDTE_ANALYSIS_HPP

// The 73mGe line in the PSI CdTe (Amptek, DT5780 channel 1, BOA beamline,
// October 2025), from our own trapezoid energies (ProcessWaveforms.cpp, through
// Analysis-Utilities' PulseHeightProcessor), calibrated in each run by the
// lines around it. What every cdte/ macro shares: the paths and runs, the
// processed trees, the EMG lineshape, and the in-run self-calibrating fit.
//
// The in-run fit places every known line at its true energy on a free scale
// x = a + b E (+ c (E - 70)^2 in the extrapolation check): Pb Ka2/Ka1 (free
// ratio), the Pb Kb1' and Kb2' groups (sub-lines at fixed relative
// intensities, free group ratio), the 95.902 keV line's Cd and Te K X-ray
// escapes at the capture line's measured position minus the X-ray energies,
// and Ge, whose energy is a fit parameter. One EMG lineshape for all lines,
// its width linear in energy, on a linear continuum. Ge's error then carries
// the calibration's own uncertainty and correlations.

#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"

#include <RtypesCore.h>
#include <TChain.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMath.h>
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

// The Ge beam runs: the CdTe behind the Ge target.
inline std::vector<TString> GeRuns() {
  std::vector<TString> runs;
  runs.push_back("PSI_Running_BeamTurnOn_20251009_1914_4");
  runs.push_back("PSI_Beam_GeTarget_20251010_1525_1");
  runs.push_back("PSI_Beam_GeTarget_20251010_1525_2");
  runs.push_back("PSI_Beam_Ge_20251011_1715_1");
  runs.push_back("PSI_Beam_Ge_20251011_1715_2");
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
// trapping, so each slice is calibrated on its own lines.
const Int_t N_RISE_SLICES = 7;
const Double_t RISE_SLICE_LO[N_RISE_SLICES] = {0, 270, 280, 290, 300, 320, 350};
const Double_t RISE_SLICE_HI[N_RISE_SLICES] = {270, 280, 290, 300,
                                               320, 350, 450};

// The trapezoid spectrum of the given runs on the Am-only scale, events with
// a 10-90% rise in [rise_lo_ns, rise_hi_ns).
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
    h->Fill(height * KEV_PER_UNIT);
  }
  return h;
}

// ---------------------------------------------------------------- lineshape

// exp(b^2) erfc(b) for b >= 0, without overflow.
inline Double_t ErfcxPositive(Double_t b) {
  if (b < 20.0)
    return TMath::Exp(b * b) * TMath::Erfc(b);
  Double_t inv2 = 1.0 / (b * b);
  return (1.0 - 0.5 * inv2 + 0.75 * inv2 * inv2) /
         (b * TMath::Sqrt(TMath::Pi()));
}

// Low-side exponentially modified Gaussian of area A (counts), mean mu, width
// sigma and tail length tau (keV), in counts per BIN-wide bin.
inline Double_t Emg(Double_t x, Double_t area, Double_t mu, Double_t sigma,
                    Double_t tau) {
  Double_t z = x - mu;
  Double_t b = (z / sigma + sigma / tau) / TMath::Sqrt2();
  Double_t value;
  if (b >= 0)
    value = TMath::Exp(-0.5 * z * z / (sigma * sigma)) * ErfcxPositive(b);
  else
    value = TMath::Exp(z / tau + 0.5 * sigma * sigma / (tau * tau)) *
            TMath::Erfc(b);
  return area * BIN * value / (2.0 * tau);
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

// --------------------------------------------------------- the capture line

// The 95.902 keV line, one EMG on a linear continuum over 90-99 keV, on the
// spectrum's own scale: where its escapes sit, and how strong it is.
struct CaptureLine {
  Double_t mu = 0, mu_err = 0, area = 0, area_err = 0, sigma = 0, tau = 0;
  Double_t chi2 = 0;
};

// [0] continuum at mu [1] its slope [2] sigma [3] tau [4] area [5] mu
inline Double_t CaptureModel(Double_t *x, Double_t *p) {
  return p[0] + p[1] * (x[0] - p[5]) + Emg(x[0], p[4], p[5], p[2], p[3]);
}

inline Double_t CaptureContinuum(Double_t x, const Double_t *p) {
  return p[0] + p[1] * (x - p[5]);
}

inline CaptureLine FitCaptureLine(TH1D *h, const TString &tag,
                                  Bool_t draw = kFALSE) {
  const Double_t lo = 90.0, hi = 99.0;
  const Double_t e_line = ReferenceLines::CD113_CAPTURE;
  CaptureLine out;
  h->GetXaxis()->SetRangeUser(e_line - 0.6, e_line + 0.6);
  Double_t top = h->GetBinCenter(h->GetMaximumBin());
  Double_t height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  TF1 *f = new TF1("capture_" + tag, CaptureModel, lo, hi, 6);
  Double_t level = h->GetBinContent(h->FindBin(hi - 0.1));
  f->SetParameters(level, 0, 0.25, 0.4, height * 0.25 * 2.5 / BIN, top + 0.1);
  f->SetParLimits(2, 0.05, 0.9);
  f->SetParLimits(3, 0.02, 5.0);
  f->SetParLimits(4, 0, 1e9);
  h->Fit(f, "QLRN0");
  h->Fit(f, "QLRN0M");
  out.mu = f->GetParameter(5);
  out.mu_err = f->GetParError(5);
  out.area = f->GetParameter(4);
  out.area_err = f->GetParError(4);
  out.sigma = f->GetParameter(2);
  out.tau = f->GetParameter(3);
  out.chi2 = f->GetChisquare() / f->GetNDF();
  if (draw) {
    std::vector<TGraph *> components;
    components.push_back(
        Component(lo, hi, kGray + 2, 2, CaptureContinuum, f->GetParameters()));
    Plot(h, f, components, lo, hi, "psi_cdte_capture_" + tag, tag);
  }
  delete f;
  return out;
}

// ------------------------------------------------------------ in-run model

enum {
  Q_A,
  Q_B,
  Q_SCALE2, // x = a + b E + c (E - 70)^2; c fixed at 0 but in the check
  Q_SIG0,
  Q_SIG1, // sigma(E) = sig0 + sig1 (E - 70)
  Q_TAU,
  Q_C0,
  Q_C1, // linear continuum in (x - 70)
  Q_GE_E,
  Q_GE_A,
  Q_KA1_A,
  Q_KA2_R,
  Q_KB1_A,
  Q_KB2_R,
  Q_ESC_A,  // the Cd K-beta escape's yield; the others tied to it
  Q_TAU_GE, // Ge's own tail, only with g_ge_own_tail
  N_Q
};

// The fit window.
Double_t g_lo = 62.0, g_hi = 89.0;
// Variants.
Bool_t g_ge_own_tail = kFALSE;
Bool_t g_free_scale2 = kFALSE;

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

// The escape tie: Cd K-alpha from the Cd Kb/Ka ratio, Te K-alpha from its
// ratio to Cd K-beta, Te K-beta at 0.30 of Te K-alpha (Geant4, captures in
// CdTe; the runs do not resolve it).
inline void SetEscapeTie(Double_t cd_kb_over_ka, Double_t te_ka_over_cd_kb) {
  Double_t series_total[4] = {1.0 / cd_kb_over_ka, 1.0, te_ka_over_cd_kb,
                              0.30 * te_ka_over_cd_kb};
  for (Int_t i = 0; i < ReferenceLines::N_ESCAPE_LINES; i++)
    g_escape_weight[i] = series_total[ReferenceLines::ESCAPE_SERIES[i]] *
                         ReferenceLines::ESCAPE_SERIES_SHARE[i];
}

inline Double_t ScalePosition(Double_t e, const Double_t *p) {
  Double_t d = e - 70.0;
  return p[Q_A] + p[Q_B] * e + p[Q_SCALE2] * d * d;
}

// The true energy at spectrum position x, inverting the scale.
inline Double_t ScaleEnergy(Double_t x, const Double_t *p) {
  Double_t e = (x - p[Q_A]) / p[Q_B];
  if (p[Q_SCALE2] == 0)
    return e;
  for (Int_t iteration = 0; iteration < 4; iteration++) {
    Double_t d = e - 70.0;
    e -= (ScalePosition(e, p) - x) / (p[Q_B] + 2.0 * p[Q_SCALE2] * d);
  }
  return e;
}

// A unit-area line at true energy e, in counts per bin.
inline Double_t Shape(Double_t x, Double_t e, const Double_t *p) {
  Double_t sigma = p[Q_SIG0] + p[Q_SIG1] * (e - 70.0);
  return Emg(x, 1.0, ScalePosition(e, p), sigma, p[Q_TAU]);
}

inline Double_t Bkg(Double_t x, const Double_t *p) {
  return p[Q_C0] + p[Q_C1] * (x - 70.0);
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
  if (!g_ge_own_tail)
    return p[Q_GE_A] * Shape(x, p[Q_GE_E], p);
  Double_t sigma = p[Q_SIG0] + p[Q_SIG1] * (p[Q_GE_E] - 70.0);
  return p[Q_GE_A] *
         Emg(x, 1.0, ScalePosition(p[Q_GE_E], p), sigma, p[Q_TAU_GE]);
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

// Staged: (1) the scale fixed at the Pb Ka1 peak's offset, ratios fixed, no
// escapes; (2) the scale released; (3) the ratios, escapes and the width's
// slope released (and Ge's own tail or the scale's curvature in those
// checks). The caller sets the capture line and the escape tie first.
inline TF1 *FitInRun(TH1D *h, const TString &tag) {
  TF1 *f = new TF1("inrun_" + tag, InRunModel, g_lo, g_hi, N_Q);
  f->SetNpx(3000);
  h->GetXaxis()->SetRangeUser(74.2, 75.6);
  Double_t ka1_mode = h->GetBinCenter(h->GetMaximumBin());
  Double_t ka1_height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  Double_t level = h->GetBinContent(h->FindBin(65.0));
  Double_t a0 = ka1_mode - ReferenceLines::PB_KA1 + 0.3;
  Double_t ka1 = (ka1_height - level) * 0.6 * 2.5 / BIN;

  f->SetParameter(Q_SIG0, 0.6);
  f->SetParLimits(Q_SIG0, 0.1, 1.5);
  f->SetParameter(Q_TAU, 1.0);
  f->SetParLimits(Q_TAU, 0.02, 5.0);
  f->SetParameter(Q_C0, level);
  f->SetParameter(Q_C1, 0);
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
  f->FixParameter(Q_SCALE2, 0.0);
  f->FixParameter(Q_KA2_R, 0.6);
  f->FixParameter(Q_KB2_R, 0.3);
  f->FixParameter(Q_ESC_A, 0.0);
  f->FixParameter(Q_TAU_GE, 1.0);
  f->FixParameter(Q_SIG1, 0.0);
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
  if (g_ge_own_tail) {
    f->ReleaseParameter(Q_TAU_GE);
    f->SetParameter(Q_TAU_GE, f->GetParameter(Q_TAU));
    f->SetParLimits(Q_TAU_GE, 0.02, 5.0);
  }
  if (g_free_scale2) {
    f->ReleaseParameter(Q_SCALE2);
    f->SetParLimits(Q_SCALE2, -0.01, 0.01);
  }
  f->ReleaseParameter(Q_SIG1);
  f->SetParLimits(Q_SIG1, -0.05, 0.05);
  h->Fit(f, "QLRN0");
  h->Fit(f, "QLRN0M");
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

} // namespace CdTe

#endif
