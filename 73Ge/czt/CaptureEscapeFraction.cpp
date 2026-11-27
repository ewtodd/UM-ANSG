// The 95.902 keV line's K X-ray escapes, measured.
//
// Captures inside the crystal: the 01/15 beam background run has no shield,
// so its line comes from captures in the CZT's own cathode layer
// (Capture96Positions.cpp), and no Pb, so the Cd K-alpha escapes stand clear
// of Pb K-alpha2. Both Cd escape pairs are fitted as one yield, Cd K-alpha
// held to Cd K-beta by ReferenceLines' series ratio; that is the measurement
// ESCAPES_INSIDE takes. Its alternatives (the photopeak width, a quadratic
// continuum) and the two series fitted apart are printed beside it.
//
// The Cd-shield background runs of 01/13: does their spectrum carry the
// in-crystal escapes or the shield's (ESCAPES_OUTSIDE, Geant4)? Pb K-alpha
// lies over the Cd K-alpha escapes there, so the Cd K-beta escape is fitted
// beside free Pb K-alpha lines, and the two complete fraction sets are
// compared by their Baker-Cousins chi2.
//
// Every run's production energies (the calibrated trees of PRODUCTION_TAG)
// are binned one comb period wide at N_RATESUB_PHASES starting points, as
// production bins them, and the starting points are averaged. The 95.902 keV
// line is production's AverageCaptureLine. Escape lines sit at the line's
// centroid minus the SRD 128 X-ray energies with ReferenceLines' shares; the
// Te escapes, too weak to fit, follow Cd K-alpha by Geant4's ratios. Binned
// Poisson likelihood.
//
// The model is also fitted to the Geant4 run of captures in the crystal, its
// deposits smeared with the escape lines' measured width: the check that the
// lines, shares and ratios production holds describe the simulated escapes.
//
// Writes results/escape_fractions.result and the dissertation's figures of
// the 01/15 escapes and of the model on the simulation
// (plots/czt/figures/czt_escapes_0115, log_czt_escapes_geant4).
//
// Usage: root -b -q 'CaptureEscapeFraction.cpp+'
#include "Constants.hpp"
#include "FigureStyle.hpp"
#include "GainTransferFit.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace {

const char *const PRODUCTION_TAG = "tau20_sigpb_esc_pxday";
const Double_t INSIDE_WINDOW_LO_KEV = 66.0, INSIDE_WINDOW_HI_KEV = 74.5;
const Double_t SHIELD_WINDOW_LO_KEV = 66.0, SHIELD_WINDOW_HI_KEV = 77.0;
const Double_t PB_CHECK_LO_KEV = 73.8, PB_CHECK_HI_KEV = 77.0;
const Double_t CONTINUUM_PIVOT_KEV = 70.0;
// Cd K-alpha escapes per Cd K-beta escape, the ratio the tied fits hold:
// Geant4's for captures inside the crystal.
const Double_t CD_KA_OVER_KB = GEANT4_INSIDE_CD_KA_OVER_CD_KB;

// Parameters of EscapeWindowModel.
enum EscapeWindowParameter {
  kContinuum = 0, // counts per bin at CONTINUUM_PIVOT_KEV
  kSlope,         // per keV
  kCurvature,     // per keV^2, fixed at zero but for the alternative
  kSigma,
  kTailHeightRatio,
  kTauOverSigma,
  kBinWidth,
  kCaptureMu,
  kCdKalphaYield, // unused while kCdKalphaOverCdKbeta is set
  kCdKbetaYield,
  kCdKalphaOverCdKbeta,  // fixed; 0 frees Cd K-alpha
  kTeKalphaOverCdKalpha, // fixed
  kTeKbetaOverTeKalpha,  // fixed
  kPbKalpha2Yield,       // fixed at zero without Pb
  kPbKalpha1Yield,
  kPbKalpha2Mu,
  kPbKalpha1Mu,
  kNParameters
};

Double_t LineWithTail(Double_t x, Double_t mu, Double_t sigma,
                      Double_t gaus_yield, Double_t tail_height_ratio,
                      Double_t tau_over_sigma, Double_t bin_width) {
  Double_t height =
      gaus_yield * bin_width / (sigma * TMath::Sqrt(2.0 * TMath::Pi()));
  Double_t at[1] = {x};
  Double_t gaus[3] = {mu, sigma, height};
  Double_t tail[6] = {mu,  sigma, height * tail_height_ratio, tau_over_sigma,
                      0.0, 0.0};
  return FittingFunctions::Gaussian(at, gaus) +
         FittingFunctions::LowTail(at, tail);
}

Double_t Continuum(Double_t x, const Double_t *p) {
  Double_t dx = x - CONTINUUM_PIVOT_KEV;
  return p[kContinuum] + p[kSlope] * dx + p[kCurvature] * dx * dx;
}

// The Gaussian yield of escape series s, numbered as ESCAPE_SERIES: Cd
// K-alpha, Cd K-beta, Te K-alpha, Te K-beta.
Double_t SeriesYield(const Double_t *p, Int_t s) {
  Double_t cd_kalpha = p[kCdKalphaOverCdKbeta] > 0
                           ? p[kCdKalphaOverCdKbeta] * p[kCdKbetaYield]
                           : p[kCdKalphaYield];
  if (s == 0)
    return cd_kalpha;
  if (s == 1)
    return p[kCdKbetaYield];
  Double_t te_kalpha = cd_kalpha * p[kTeKalphaOverCdKalpha];
  if (s == 2)
    return te_kalpha;
  return te_kalpha * p[kTeKbetaOverTeKalpha];
}

// The escapes of one series (0-3), or of all four with series = -1.
Double_t Escapes(Double_t x, const Double_t *p, Int_t series) {
  Double_t sum = 0;
  for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
    if (series >= 0 && ESCAPE_SERIES[i] != series)
      continue;
    Double_t yield =
        SeriesYield(p, ESCAPE_SERIES[i]) * GEANT4_INSIDE_SERIES_SHARE[i];
    sum += LineWithTail(x, p[kCaptureMu] - ESCAPE_XRAY_KEV[i], p[kSigma], yield,
                        p[kTailHeightRatio], p[kTauOverSigma], p[kBinWidth]);
  }
  return sum;
}

Double_t PbLines(Double_t x, const Double_t *p) {
  return LineWithTail(x, p[kPbKalpha2Mu], p[kSigma], p[kPbKalpha2Yield],
                      p[kTailHeightRatio], p[kTauOverSigma], p[kBinWidth]) +
         LineWithTail(x, p[kPbKalpha1Mu], p[kSigma], p[kPbKalpha1Yield],
                      p[kTailHeightRatio], p[kTauOverSigma], p[kBinWidth]);
}

Double_t EscapeWindowModel(Double_t *x, Double_t *p) {
  return Continuum(x[0], p) + Escapes(x[0], p, -1) + PbLines(x[0], p);
}

// How the escapes enter a window fit.
enum EscapeMode {
  kCdTied,  // one free Cd yield, K-alpha held to K-beta by the series ratio
  kCdApart, // Cd K-alpha and K-beta free apart
  kFixedSet // every series at a given set of fractions
};

struct WindowFitConfig {
  Double_t lo = INSIDE_WINDOW_LO_KEV, hi = INSIDE_WINDOW_HI_KEV;
  Bool_t quadratic = kFALSE;
  Double_t fixed_sigma = -1.0; // <= 0 frees the width
  // Free Pb K-alpha lines in the window, with the escapes' width and tail.
  Bool_t with_pb = kFALSE;
  EscapeMode mode = kCdTied;
  Double_t te_kalpha_over_cd_kalpha = GEANT4_INSIDE_TE_KA_OVER_CD_KA;
  Double_t te_kbeta_over_te_kalpha = GEANT4_INSIDE_TE_KB_OVER_TE_KA;
  const EscapeFractions *fixed_set = nullptr;
};

struct EscapeFit {
  Bool_t valid = kFALSE;
  // The two fitted series, Cd K-alpha and Cd K-beta.
  Double_t series_yield[2] = {0, 0};
  Double_t series_yield_err[2] = {0, 0};
  Double_t sigma = 0;
  Double_t chi2 = 0; // Pearson, for the figure
  Int_t ndf = 0;
  Double_t baker_cousins = 0;
  std::vector<Double_t> parameters;
};

// Pearson chi2 and the Baker-Cousins likelihood chi2 of the model over the
// window's bins.
void WindowChi2(TH1D *h, TF1 &model, Double_t lo, Double_t hi, Int_t n_free,
                Double_t &chi2, Int_t &ndf, Double_t &baker_cousins) {
  chi2 = 0;
  baker_cousins = 0;
  Int_t n_bins = 0;
  for (Int_t b = h->FindBin(lo); b <= h->FindBin(hi); b++) {
    Double_t x = h->GetBinCenter(b);
    if (x < lo || x > hi)
      continue;
    Double_t expected = model.Eval(x);
    if (expected <= 0)
      continue;
    Double_t observed = h->GetBinContent(b);
    chi2 += (observed - expected) * (observed - expected) / expected;
    baker_cousins += 2.0 * (expected - observed);
    if (observed > 0)
      baker_cousins += 2.0 * observed * TMath::Log(observed / expected);
    n_bins++;
  }
  ndf = n_bins - n_free;
}

Double_t BinAt(TH1D *h, Double_t x) { return h->GetBinContent(h->FindBin(x)); }

EscapeFit FitEscapeWindow(TH1D *h, const CaptureLine &line,
                          const RateSubShape &shape, const WindowFitConfig &cfg,
                          const TString &tag) {
  EscapeFit out;
  TF1 model("escape_window_" + tag, EscapeWindowModel, cfg.lo, cfg.hi,
            kNParameters);
  Double_t bin = h->GetBinWidth(1);
  Double_t level_lo = BinAt(h, cfg.lo + 0.5);
  Double_t level_hi = BinAt(h, cfg.with_pb ? cfg.hi - 0.3 : cfg.hi - 2.5);
  Double_t span = (cfg.with_pb ? cfg.hi - 0.8 : cfg.hi - 3.0) - cfg.lo;
  model.SetParameter(kContinuum, 0.5 * (level_lo + level_hi));
  model.SetParameter(kSlope, (level_hi - level_lo) / span);
  model.SetParameter(kCurvature, 0.0);
  if (!cfg.quadratic)
    model.FixParameter(kCurvature, 0.0);
  model.SetParameter(kSigma, 0.85);
  model.SetParLimits(kSigma, 0.5, 1.3);
  if (cfg.fixed_sigma > 0)
    model.FixParameter(kSigma, cfg.fixed_sigma);
  model.FixParameter(kTailHeightRatio, shape.ge_tail_height_ratio);
  model.FixParameter(kTauOverSigma, shape.ge_tau_over_sigma);
  model.FixParameter(kBinWidth, bin);
  model.FixParameter(kCaptureMu, line.mu);

  model.FixParameter(kCdKalphaOverCdKbeta, 0.0);
  model.FixParameter(kTeKalphaOverCdKalpha, cfg.te_kalpha_over_cd_kalpha);
  model.FixParameter(kTeKbetaOverTeKalpha, cfg.te_kbeta_over_te_kalpha);
  model.SetParameter(kCdKalphaYield,
                     line.gaus_yield * ESCAPES_INSIDE.series[0]);
  model.SetParameter(kCdKbetaYield, line.gaus_yield * ESCAPES_INSIDE.series[1]);
  if (cfg.mode == kCdTied) {
    model.FixParameter(kCdKalphaYield, 0.0);
    model.FixParameter(kCdKalphaOverCdKbeta, CD_KA_OVER_KB);
  } else if (cfg.mode == kFixedSet) {
    const EscapeFractions &f = *cfg.fixed_set;
    model.FixParameter(kCdKalphaYield, line.gaus_yield * f.series[0]);
    model.FixParameter(kCdKbetaYield, line.gaus_yield * f.series[1]);
    model.FixParameter(kTeKalphaOverCdKalpha, f.series[2] / f.series[0]);
    model.FixParameter(kTeKbetaOverTeKalpha, f.series[3] / f.series[2]);
  }

  if (cfg.with_pb) {
    Double_t continuum0 = 0.5 * (level_lo + level_hi);
    Double_t to_yield = 0.85 * TMath::Sqrt(2.0 * TMath::Pi()) / bin;
    model.SetParameter(
        kPbKalpha2Yield,
        TMath::Max(BinAt(h, ReferenceLines::PB_KA2) - continuum0, 1.0) *
            to_yield);
    model.SetParameter(
        kPbKalpha1Yield,
        TMath::Max(BinAt(h, ReferenceLines::PB_KA1) - continuum0, 1.0) *
            to_yield);
    model.SetParameter(kPbKalpha2Mu, ReferenceLines::PB_KA2);
    model.SetParLimits(kPbKalpha2Mu, ReferenceLines::PB_KA2 - 0.3,
                       ReferenceLines::PB_KA2 + 0.3);
    model.SetParameter(kPbKalpha1Mu, ReferenceLines::PB_KA1);
    model.SetParLimits(kPbKalpha1Mu, ReferenceLines::PB_KA1 - 0.3,
                       ReferenceLines::PB_KA1 + 0.3);
  } else {
    model.FixParameter(kPbKalpha2Yield, 0.0);
    model.FixParameter(kPbKalpha1Yield, 0.0);
    model.FixParameter(kPbKalpha2Mu, ReferenceLines::PB_KA2);
    model.FixParameter(kPbKalpha1Mu, ReferenceLines::PB_KA1);
  }

  TFitResultPtr fr = h->Fit(&model, "S L Q R N");
  Bool_t width_railed = cfg.fixed_sigma <= 0 && ParameterAtLimit(model, kSigma);
  if (!fr.Get() || !fr->IsValid() || fr->Status() != 0 ||
      fr->CovMatrixStatus() != 3 || width_railed) {
    std::cerr << "ERROR: the escape window of " << tag
              << " did not converge cleanly" << std::endl;
    return out;
  }
  const Double_t *p = model.GetParameters();
  out.series_yield[0] = SeriesYield(p, 0);
  out.series_yield_err[0] =
      cfg.mode == kCdTied ? CD_KA_OVER_KB * model.GetParError(kCdKbetaYield)
                          : model.GetParError(kCdKalphaYield);
  out.series_yield[1] = model.GetParameter(kCdKbetaYield);
  out.series_yield_err[1] = model.GetParError(kCdKbetaYield);
  out.sigma = model.GetParameter(kSigma);
  WindowChi2(h, model, cfg.lo, cfg.hi, model.GetNumberFreeParameters(),
             out.chi2, out.ndf, out.baker_cousins);
  for (Int_t i = 0; i < kNParameters; i++)
    out.parameters.push_back(model.GetParameter(i));
  out.valid = kTRUE;
  return out;
}

// The escape window's model, Cd K-alpha escapes only, plus Pb K-alpha1 with
// its Gaussian yield in the parameter after the window's.
Double_t PbCheckModel(Double_t *x, Double_t *p) {
  return Continuum(x[0], p) + Escapes(x[0], p, 0) +
         LineWithTail(x[0], ReferenceLines::PB_KA1, p[kSigma], p[kNParameters],
                      p[kTailHeightRatio], p[kTauOverSigma], p[kBinWidth]);
}

// Pb K-alpha1 above the 01/15 window, on a linear continuum with the high
// side of the Cd K-alpha escapes fixed at the window's fit.
void CheckPbKalpha1(TH1D *h, const EscapeFit &fit) {
  const Int_t n = kNParameters + 1;
  TF1 model("pb_check", PbCheckModel, PB_CHECK_LO_KEV, PB_CHECK_HI_KEV, n);
  for (Int_t i = 0; i < kNParameters; i++)
    model.FixParameter(i, fit.parameters[i]);
  model.ReleaseParameter(kContinuum);
  model.ReleaseParameter(kSlope);
  model.SetParameter(kNParameters, 0.0);
  TFitResultPtr fr = h->Fit(&model, "S L Q R N");
  std::cout << std::fixed << std::setprecision(1) << "Pb K-alpha1 check ("
            << PB_CHECK_LO_KEV << "-" << PB_CHECK_HI_KEV
            << " keV): Gaussian yield " << std::setprecision(0)
            << model.GetParameter(kNParameters) << " +/- "
            << model.GetParError(kNParameters) << " against Cd K-alpha escapes "
            << fit.series_yield[0]
            << (fr.Get() && fr->IsValid() ? "" : " (fit not valid)")
            << std::endl;
}

// "Name value error fixed" lines of a saved RooFitUtils fit state.
Bool_t ReadSavedParameter(const TString &path, const TString &name,
                          Double_t &value) {
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return kFALSE;
  }
  std::string text;
  while (std::getline(in, text)) {
    std::istringstream fields(text);
    std::string key;
    Double_t number = 0;
    fields >> key >> number;
    if (TString(key.c_str()) == name) {
      value = number;
      return kTRUE;
    }
  }
  std::cerr << "ERROR: no " << name << " in " << path << std::endl;
  return kFALSE;
}

// The tail parameters of peak `index` of a saved state, as PeakFitResult
// carries them: the tail amplitude is the Gaussian yield times the
// full-range yield ratio the state stores.
Bool_t ReadSavedTail(const TString &path, Int_t index, PeakFitResult &peak) {
  Double_t gaus = 0, ratio = 0, rho = 0;
  if (!ReadSavedParameter(path, Form("GausAmplitude%d", index), gaus) ||
      !ReadSavedParameter(path, Form("LowExpTailAmplitude%d", index), ratio) ||
      !ReadSavedParameter(path, Form("LowExpTailRatio%d", index), rho))
    return kFALSE;
  peak.gaus_amplitude = gaus;
  peak.low_exp_tail_amplitude = gaus * ratio;
  peak.low_exp_tail_ratio = rho;
  return kTRUE;
}

TF1 *ReadCalibration(const TString &path) {
  TFile file(path, "READ");
  TF1 *cal =
      file.IsZombie() ? nullptr : dynamic_cast<TF1 *>(file.Get("calibration"));
  if (!cal) {
    std::cerr << "ERROR: no calibration in " << path << std::endl;
    return nullptr;
  }
  return (TF1 *)cal->Clone();
}

std::vector<TH1D *> BuildCountHists(const TString &path, const TString &name,
                                    Double_t lo, Double_t hi, Double_t bin) {
  std::vector<TH1D *> hists;
  TFile file(path, "READ");
  TTree *tree = file.IsZombie()
                    ? nullptr
                    : dynamic_cast<TTree *>(file.Get("calibrated_tree"));
  if (!tree) {
    std::cerr << "ERROR: no calibrated_tree in " << path << std::endl;
    return hists;
  }
  Int_t n_bins = (Int_t)TMath::Nint((hi - lo) / bin);
  for (Int_t k = 0; k < N_RATESUB_PHASES; k++) {
    Double_t start = lo + k * bin / N_RATESUB_PHASES;
    TH1D *h =
        new TH1D(Form("%s_phase%d", name.Data(), k), ";Energy [keV];Counts",
                 n_bins, start, start + n_bins * bin);
    h->SetDirectory(nullptr);
    hists.push_back(h);
  }
  Float_t energy = 0;
  tree->SetBranchAddress("depositedEnergykeV", &energy);
  Long64_t n_entries = tree->GetEntries();
  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    for (Int_t k = 0; k < N_RATESUB_PHASES; k++)
      hists[k]->Fill(energy);
  }
  tree->ResetBranchAddresses();
  return hists;
}

void DeleteHists(std::vector<TH1D *> &hists) {
  for (size_t k = 0; k < hists.size(); k++)
    delete hists[k];
  hists.clear();
}

std::vector<EscapeFit> FitPhases(const std::vector<TH1D *> &hists,
                                 const CaptureLine &line,
                                 const RateSubShape &shape,
                                 const WindowFitConfig &cfg,
                                 const TString &tag) {
  std::vector<EscapeFit> fits;
  for (size_t k = 0; k < hists.size(); k++)
    fits.push_back(FitEscapeWindow(hists[k], line, shape, cfg,
                                   Form("%s_phase%d", tag.Data(), (Int_t)k)));
  return fits;
}

// One series' fraction of the line's Gaussian yield over the phases: the
// mean, the mean fit error (with the line yield's) and the phases' rms.
struct Fraction {
  Bool_t valid = kFALSE;
  Double_t value = 0, stat = 0, phase_rms = 0;
};

Fraction SeriesFraction(const std::vector<EscapeFit> &fits,
                        const CaptureLine &line, Int_t series) {
  Fraction out;
  Double_t sum = 0, sum2 = 0, sum_err = 0;
  Int_t n = 0;
  for (size_t k = 0; k < fits.size(); k++) {
    if (!fits[k].valid)
      continue;
    Double_t f = fits[k].series_yield[series] / line.gaus_yield;
    Double_t f_err = TMath::Sqrt(
        TMath::Power(fits[k].series_yield_err[series] / line.gaus_yield, 2) +
        TMath::Power(f * line.gaus_yield_err / line.gaus_yield, 2));
    sum += f;
    sum2 += f * f;
    sum_err += f_err;
    n++;
  }
  if (n == 0)
    return out;
  out.valid = kTRUE;
  out.value = sum / n;
  out.stat = sum_err / n;
  out.phase_rms =
      TMath::Sqrt(TMath::Max(sum2 / n - out.value * out.value, 0.0));
  return out;
}

Double_t MeanBakerCousins(const std::vector<EscapeFit> &fits) {
  Double_t sum = 0;
  Int_t n = 0;
  for (size_t k = 0; k < fits.size(); k++) {
    if (!fits[k].valid)
      continue;
    sum += fits[k].baker_cousins;
    n++;
  }
  return n > 0 ? sum / n : 0.0;
}

void PrintFraction(const TString &label, const char *series,
                   const Fraction &f) {
  if (!f.valid)
    return;
  std::cout << label << " " << series << " escape fraction: " << std::fixed
            << std::setprecision(4) << f.value << " +/- " << f.stat
            << " (phase rms " << f.phase_rms << ")" << std::endl;
}

TGraph *CurveOf(const std::vector<Double_t> &p, Int_t series, Double_t lo,
                Double_t hi) {
  const Int_t n_points = 600;
  TGraph *g = new TGraph(n_points);
  for (Int_t i = 0; i < n_points; i++) {
    Double_t x = lo + (hi - lo) * i / (n_points - 1);
    g->SetPoint(i, x, Escapes(x, p.data(), series));
  }
  return g;
}

// A fit curve as the production fit figures draw theirs (FitFigure.cpp).
void StyleCurve(TGraph *g, Color_t color, Style_t style) {
  PlottingUtils::ConfigureGraph(g, color);
  g->SetLineStyle(style);
  g->SetLineWidth(2 * PlottingUtils::GetLineWidth());
}

// The fitted continuum alone (series = -2), or with one escape series (0-3)
// or all four (-1) on it.
TGraph *CurveOnContinuum(const std::vector<Double_t> &p, Int_t series,
                         Double_t lo, Double_t hi) {
  const Int_t n_points = 600;
  TGraph *g = new TGraph(n_points);
  for (Int_t i = 0; i < n_points; i++) {
    Double_t x = lo + (hi - lo) * i / (n_points - 1);
    Double_t lines = (series == -2) ? 0.0 : Escapes(x, p.data(), series);
    g->SetPoint(i, x, Continuum(x, p.data()) + lines);
  }
  return g;
}

// The window's counts as measured, the fitted line drawn on its own and each
// escape series drawn on it.
void DrawEscapeWindow(TH1D *h, const EscapeFit &fit, const TString &name) {
  const std::vector<Double_t> &p = fit.parameters;
  const Double_t lo = INSIDE_WINDOW_LO_KEV, hi = INSIDE_WINDOW_HI_KEV;
  Double_t bin_ev = 1000.0 * h->GetBinWidth(1);
  TGraphErrors *points = new TGraphErrors();
  TGraph *pulls = new TGraph();
  for (Int_t b = 1; b <= h->GetNbinsX(); b++) {
    Double_t x = h->GetBinCenter(b);
    if (x < lo || x > hi)
      continue;
    Double_t counts = h->GetBinContent(b);
    Double_t error = TMath::Sqrt(TMath::Max(counts, 1.0));
    Double_t expected = Continuum(x, p.data()) + Escapes(x, p.data(), -1);
    Int_t i = points->GetN();
    points->SetPoint(i, x, counts);
    points->SetPointError(i, 0.5 * h->GetBinWidth(1), error);
    pulls->SetPoint(pulls->GetN(), x, (counts - expected) / error);
  }

  TGraph *total = CurveOnContinuum(p, -1, lo, hi);
  TGraph *continuum = CurveOnContinuum(p, -2, lo, hi);
  TGraph *cd_kbeta = CurveOnContinuum(p, 1, lo, hi);
  TGraph *cd_kalpha = CurveOnContinuum(p, 0, lo, hi);
  TGraph *te_kalpha = CurveOnContinuum(p, 2, lo, hi);

  Double_t y_min = points->GetY()[0], y_max = points->GetY()[0];
  for (Int_t i = 0; i < points->GetN(); i++) {
    y_min = TMath::Min(y_min, points->GetY()[i] - points->GetEY()[i]);
    y_max = TMath::Max(y_max, points->GetY()[i] + points->GetEY()[i]);
  }

  std::cout << "FIGURE " << name << ": Pearson chi2 " << std::fixed
            << std::setprecision(1) << fit.chi2 << " / " << fit.ndf << " = "
            << std::setprecision(3) << fit.chi2 / fit.ndf << std::endl;
  std::vector<TString> labels;
  labels.push_back("Active background run");
  labels.push_back(Form("Fit, #chi^{2}/ndf = %.3f", fit.chi2 / fit.ndf));
  labels.push_back("Linear continuum");
  labels.push_back("Cd K#beta escapes");
  labels.push_back("Te K#alpha escapes");
  labels.push_back("Cd K#alpha escapes");

  const Int_t canvas_w = 1400, plot_h = 800, panel_h = 300;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(canvas_w, plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(canvas, panel_h, upper, lower);
  upper->cd();
  TH1D *frame = new TH1D("escape_frame", "", 100, lo, hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(
      frame, kBlack, Form(";Energy [keV];Counts / %.0f eV", bin_ev));
  // The legend's six rows above the spectrum.
  frame->SetMinimum(y_min - 0.1 * (y_max - y_min));
  frame->SetMaximum(y_max + 1.4 * (y_max - y_min));
  // Counts in thousands, as on the production fit figures.
  frame->GetYaxis()->SetMaxDigits(3);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw("AXIS");
  PlottingUtils::ConfigureGraph(points, kBlack, "");
  points->SetMarkerStyle(20);
  points->SetMarkerSize(0.9);
  StyleCurve(continuum, kGray + 2, 2);
  StyleCurve(cd_kbeta, kAzure + 2, 1);
  StyleCurve(te_kalpha, kGreen + 2, 1);
  StyleCurve(cd_kalpha, kMagenta + 1, 1);
  StyleCurve(total, kRed + 1, 1);
  continuum->Draw("L SAME");
  cd_kalpha->Draw("L SAME");
  te_kalpha->Draw("L SAME");
  cd_kbeta->Draw("L SAME");
  total->Draw("L SAME");
  points->Draw("P SAME");

  TLegend *leg = FigureStyle::AddFrameLegend(upper, labels, kTRUE, kTRUE, 70.0);
  leg->AddEntry(points, labels[0], "lep");
  leg->AddEntry(total, labels[1], "l");
  leg->AddEntry(continuum, labels[2], "l");
  leg->AddEntry(cd_kbeta, labels[3], "l");
  leg->AddEntry(te_kalpha, labels[4], "l");
  leg->AddEntry(cd_kalpha, labels[5], "l");

  lower->cd();
  TH1D *pull_frame = FigureStyle::DrawPullPanel(pulls, lo, hi, "Energy [keV]");
  pulls->SetMarkerSize(1.1);
  PlottingUtils::ScaleFigure(canvas, canvas_w, plot_h);
  // Each title beside its own labels.
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  std::vector<TString> y_labels;
  y_labels.push_back("000");
  FigureStyle::FitTitleToLabels(upper, frame, y_labels);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

// ------------------------------------- the escape model on the simulation

// The Geant4 run of captures in the crystal whose escape ratios the fits
// hold (geant4/Escape96Plots.cpp counts them), its deposits smeared event by
// event with the escape lines' measured width and fitted with production's
// escape model: the lines at the capture line's energy (95.902 keV, as
// simulated) minus the SRD 128 X-ray energies with Geant4's shares, Cd
// K-alpha held to Cd K-beta and the Te series to Cd K-alpha by Geant4's
// ratios, one free yield on a linear continuum. The smearing is Gaussian, so
// the lines are too (no tail). As a check, the two Cd series are fitted apart.
const char *const GEANT4_IN_CRYSTAL_DIR = "/root_files/geant4/czt96_surface";
const Double_t GEANT4_FIT_LO_KEV = 60.0, GEANT4_FIT_HI_KEV = 80.0;
const Double_t GEANT4_BIN_KEV = 0.25;
const UInt_t GEANT4_SMEAR_SEED = 20261126;
// The spectrum's top at this share of the log frame (the legend is beside it).
const Double_t GEANT4_DATA_SHARE = 0.93;
const Double_t GE_LINE_KEV = ReferenceLines::GE73_ENSDF;

struct Geant4Fit {
  Bool_t valid = kFALSE;
  Double_t chi2 = 0;
  Int_t ndf = 0;
  Double_t cd_ka_over_kb_apart = 0, cd_ka_over_kb_apart_err = 0;
  std::vector<Double_t> parameters;
};

// The run's deposits in every crystal, each moved by a Gaussian of width
// sigma (keV), in [lo, hi].
TH1D *SmearedGeant4Run(const TString &dir, Double_t sigma, Double_t lo,
                       Double_t hi, Double_t bin_kev) {
  const Int_t n_bins = (Int_t)((hi - lo) / bin_kev + 0.5);
  TH1D *h = new TH1D("geant4_in_crystal_smeared", "", n_bins, lo, hi);
  h->SetDirectory(nullptr);
  TRandom3 rng(GEANT4_SMEAR_SEED);
  for (Int_t crystal = 0; crystal < 4; crystal++) {
    TChain chain(Form("CZT_%d", crystal));
    // The sim merges its per-thread files into output0.root; older runs
    // kept one file per thread.
    if (gSystem->AccessPathName(dir + "/output0.root"))
      chain.Add(dir + "/output0_t*.root");
    else
      chain.Add(dir + "/output0.root");
    if (chain.GetEntries() == 0)
      continue;
    Double_t edep = 0;
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fEDep", 1);
    chain.SetBranchAddress("fEDep", &edep);
    const Long64_t n_entries = chain.GetEntries();
    for (Long64_t i = 0; i < n_entries; i++) {
      chain.GetEntry(i);
      if (edep < lo - 6.0 * sigma || edep > hi + 6.0 * sigma)
        continue;
      h->Fill(edep + rng.Gaus(0.0, sigma));
    }
  }
  return h;
}

// Production's escape model on the smeared run: Cd tied (cd_apart false) or
// the two Cd series free apart. Binned Poisson likelihood.
TF1 *FitGeant4Escapes(TH1D *h, Double_t sigma, Bool_t cd_apart) {
  TF1 *model = new TF1(cd_apart ? "geant4_escapes_apart" : "geant4_escapes",
                       EscapeWindowModel, GEANT4_FIT_LO_KEV, GEANT4_FIT_HI_KEV,
                       kNParameters);
  std::vector<Double_t> p(kNParameters, 0.0);
  p[kContinuum] = BinAt(h, CONTINUUM_PIVOT_KEV - 3.0);
  p[kSigma] = sigma;
  p[kTauOverSigma] = 1.0; // no tail: kTailHeightRatio stays zero
  p[kBinWidth] = h->GetBinWidth(1);
  p[kCaptureMu] = ReferenceLines::CD113_CAPTURE;
  p[kCdKbetaYield] = 0.2 * h->Integral(h->FindBin(68.5), h->FindBin(71.0));
  p[kCdKalphaYield] = CD_KA_OVER_KB * p[kCdKbetaYield];
  p[kCdKalphaOverCdKbeta] = cd_apart ? 0.0 : CD_KA_OVER_KB;
  p[kTeKalphaOverCdKalpha] = GEANT4_INSIDE_TE_KA_OVER_CD_KA;
  p[kTeKbetaOverTeKalpha] = GEANT4_INSIDE_TE_KB_OVER_TE_KA;
  p[kPbKalpha2Mu] = ReferenceLines::PB_KA2;
  p[kPbKalpha1Mu] = ReferenceLines::PB_KA1;
  model->SetParameters(p.data());
  for (Int_t k = 0; k < kNParameters; k++)
    model->FixParameter(k, p[k]);
  model->ReleaseParameter(kContinuum);
  model->ReleaseParameter(kSlope);
  model->ReleaseParameter(kCdKbetaYield);
  if (cd_apart)
    model->ReleaseParameter(kCdKalphaYield);
  h->Fit(model, "QLRN0");
  return model;
}

void DrawGeant4Fit(TH1D *h, const Geant4Fit &fit, Double_t sigma,
                   const TString &name) {
  const std::vector<Double_t> &p = fit.parameters;
  const Double_t lo = GEANT4_FIT_LO_KEV, hi = GEANT4_FIT_HI_KEV;
  TGraphErrors *points = new TGraphErrors();
  TGraph *pulls = new TGraph();
  Double_t data_min = 0, data_max = 0;
  for (Int_t b = 1; b <= h->GetNbinsX(); b++) {
    const Double_t x = h->GetBinCenter(b);
    const Double_t counts = h->GetBinContent(b);
    const Double_t error = TMath::Sqrt(TMath::Max(counts, 1.0));
    const Double_t expected = Continuum(x, p.data()) + Escapes(x, p.data(), -1);
    const Int_t i = points->GetN();
    points->SetPoint(i, x, counts);
    points->SetPointError(i, 0.5 * h->GetBinWidth(1), error);
    pulls->SetPoint(pulls->GetN(), x, (counts - expected) / error);
    data_min = (i == 0) ? counts : TMath::Min(data_min, counts);
    data_max = TMath::Max(data_max, counts);
  }

  const Int_t n_curve = 800;
  TGraph *total = new TGraph(n_curve);
  TGraph *continuum = new TGraph(n_curve);
  for (Int_t i = 0; i < n_curve; i++) {
    const Double_t x = lo + (hi - lo) * i / (n_curve - 1);
    total->SetPoint(i, x, Continuum(x, p.data()) + Escapes(x, p.data(), -1));
    continuum->SetPoint(i, x, Continuum(x, p.data()));
  }
  const Int_t n_series = 4;
  const char *series_labels[n_series] = {
      "Cd K#alpha escapes", "Cd K#beta escapes", "Te K#alpha escapes",
      "Te K#beta escapes"};
  const Int_t series_colors[n_series] = {kMagenta + 1, kAzure + 2, kGreen + 2,
                                         kOrange + 7};
  TGraph *series[n_series];
  for (Int_t s = 0; s < n_series; s++)
    series[s] = CurveOf(p, s, lo, hi);

  std::cout << "FIGURE " << name << ": Pearson chi2 " << std::fixed
            << std::setprecision(1) << fit.chi2 << " / " << fit.ndf << " = "
            << std::setprecision(3) << fit.chi2 / fit.ndf << std::endl;
  std::vector<TString> labels;
  labels.push_back(
      Form("Geant4, captures in the crystal, #sigma = %.0f eV", 1000 * sigma));
  labels.push_back(Form("Fit, #chi^{2}/ndf = %.3f", fit.chi2 / fit.ndf));
  labels.push_back("Continuum");
  for (Int_t s = 0; s < n_series; s++)
    labels.push_back(series_labels[s]);
  labels.push_back("^{73}Ge 68.752 keV #gamma ray");

  // Eight entries: the legend beside the plot and its pull panel.
  const Int_t plot_h = 800, panel_h = 300;
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels),
      plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(plot, panel_h, upper, lower);
  upper->cd();
  TH1D *frame = new TH1D("geant4_escape_frame", "", 100, lo, hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(
      frame, kBlack,
      Form(";Deposited Energy [keV];Counts / %.0f eV",
           1000.0 * h->GetBinWidth(1)));
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw("AXIS");
  // Log only: labels at 1 and 5 x 10^n chosen from the data.
  const std::vector<TString> y_labels = FigureStyle::ConfigureLogYAxis(
      frame, upper, 0.5 * data_min, data_max, GEANT4_DATA_SHARE);
  PlottingUtils::ConfigureGraph(points, kBlack, "");
  points->SetMarkerStyle(20);
  points->SetMarkerSize(0.9);
  StyleCurve(continuum, kGray + 2, 2);
  for (Int_t s = 0; s < n_series; s++)
    StyleCurve(series[s], series_colors[s], 1);
  StyleCurve(total, kRed + 1, 1);
  TLine *ge = new TLine(GE_LINE_KEV, frame->GetMinimum(), GE_LINE_KEV,
                        frame->GetMaximum());
  ge->SetLineColor(kViolet);
  ge->SetLineStyle(2);
  ge->SetLineWidth(PlottingUtils::GetLineWidth());
  ge->Draw();
  continuum->Draw("L SAME");
  for (Int_t s = 0; s < n_series; s++)
    series[s]->Draw("L SAME");
  total->Draw("L SAME");
  points->Draw("P SAME");

  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  leg->AddEntry(points, labels[0], "lep");
  leg->AddEntry(total, labels[1], "l");
  leg->AddEntry(continuum, labels[2], "l");
  for (Int_t s = 0; s < n_series; s++)
    leg->AddEntry(series[s], labels[3 + s], "l");
  leg->AddEntry(ge, labels[3 + n_series], "l");

  lower->cd();
  TH1D *pull_frame =
      FigureStyle::DrawPullPanel(pulls, lo, hi, "Deposited Energy [keV]");
  pulls->SetMarkerSize(1.1);
  PlottingUtils::ScaleFigure(plot);
  // Each title beside its own labels.
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  FigureStyle::FitTitleToLabels(upper, frame, y_labels);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLOG);
}

Geant4Fit FitGeant4InCrystal(const TString &project_root, Double_t sigma) {
  Geant4Fit out;
  TH1D *h =
      SmearedGeant4Run(project_root + GEANT4_IN_CRYSTAL_DIR, sigma,
                       GEANT4_FIT_LO_KEV, GEANT4_FIT_HI_KEV, GEANT4_BIN_KEV);
  if (h->GetEntries() == 0) {
    std::cout << "ERROR: no Geant4 output in " << project_root
              << GEANT4_IN_CRYSTAL_DIR << std::endl;
    delete h;
    return out;
  }
  TF1 *tied = FitGeant4Escapes(h, sigma, kFALSE);
  Double_t baker_cousins = 0;
  WindowChi2(h, *tied, GEANT4_FIT_LO_KEV, GEANT4_FIT_HI_KEV, 3, out.chi2,
             out.ndf, baker_cousins);
  out.parameters.assign(tied->GetParameters(),
                        tied->GetParameters() + kNParameters);
  TF1 *apart = FitGeant4Escapes(h, sigma, kTRUE);
  const Double_t ka = apart->GetParameter(kCdKalphaYield);
  const Double_t kb = apart->GetParameter(kCdKbetaYield);
  out.cd_ka_over_kb_apart = ka / kb;
  out.cd_ka_over_kb_apart_err =
      out.cd_ka_over_kb_apart *
      TMath::Sqrt(TMath::Power(apart->GetParError(kCdKalphaYield) / ka, 2) +
                  TMath::Power(apart->GetParError(kCdKbetaYield) / kb, 2));
  std::cout << std::fixed << "Geant4 in-crystal escapes, smeared with sigma "
            << std::setprecision(4) << sigma
            << " keV: " << (Long64_t)h->GetEntries() << " deposits in "
            << std::setprecision(1) << GEANT4_FIT_LO_KEV << "-"
            << GEANT4_FIT_HI_KEV << " keV; production's model chi2 " << out.chi2
            << " / " << out.ndf << "; Cd fitted apart: K-alpha / "
            << "K-beta " << std::setprecision(3) << out.cd_ka_over_kb_apart
            << " +/- " << out.cd_ka_over_kb_apart_err << " (held at "
            << CD_KA_OVER_KB << ")" << std::endl;
  out.valid = kTRUE;
  DrawGeant4Fit(h, out, sigma, "czt_escapes_geant4");
  delete tied;
  delete apart;
  delete h;
  return out;
}

// ------------------------------------------------ unbinned in-crystal fit

// The calibrated energies of a run's events inside [lo, hi].
std::vector<Double_t> ReadWindowEvents(const TString &path, Double_t lo,
                                       Double_t hi) {
  std::vector<Double_t> events;
  TFile file(path, "READ");
  TTree *tree = file.IsZombie()
                    ? nullptr
                    : dynamic_cast<TTree *>(file.Get("calibrated_tree"));
  if (!tree) {
    std::cerr << "ERROR: no calibrated_tree in " << path << std::endl;
    return events;
  }
  Float_t energy = 0;
  tree->SetBranchAddress("depositedEnergykeV", &energy);
  Long64_t n_entries = tree->GetEntries();
  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    if (energy >= lo && energy <= hi)
      events.push_back(energy);
  }
  tree->ResetBranchAddresses();
  return events;
}

// The density per keV of the escapes that come with one unit of Gaussian
// yield, series s weighted by series_weight[s].
struct EscapeTemplate {
  Double_t capture_mu = 0, sigma = 0;
  RateSubShape shape;
  Double_t series_weight[4] = {0, 0, 0, 0};

  Double_t operator()(Double_t x) const {
    Double_t sum = 0;
    for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
      Double_t weight = series_weight[ESCAPE_SERIES[i]];
      if (weight == 0)
        continue;
      sum += weight * GEANT4_INSIDE_SERIES_SHARE[i] *
             LineWithTail(x, capture_mu - ESCAPE_XRAY_KEV[i], sigma, 1.0,
                          shape.ge_tail_height_ratio, shape.ge_tau_over_sigma,
                          1.0);
    }
    return sum;
  }

  // Simpson's rule; the shape is smooth on this step.
  Double_t Integral(Double_t lo, Double_t hi) const {
    const Int_t n_steps = 4000;
    Double_t h = (hi - lo) / n_steps;
    Double_t sum = (*this)(lo) + (*this)(hi);
    for (Int_t k = 1; k < n_steps; k++)
      sum += (k % 2 == 1 ? 4.0 : 2.0) * (*this)(lo + k * h);
    return sum * h / 3.0;
  }
};

// The extended unbinned negative log-likelihood of a model that is a sum of
// fixed shapes with free coefficients: sum_j a_j I_j - sum_i ln(sum_j a_j
// phi_j(x_i)), phi_j the shapes at each event and I_j their integrals over
// the window. It is convex in the coefficients, with gradient
// I_j - sum_i phi_j / f and Hessian sum_i phi_j phi_k / f^2 (f the model at
// x_i), so Newton's method reaches its minimum exactly, and the inverse
// Hessian there is the covariance.
struct LinearUnbinnedNLL {
  std::vector<std::vector<Double_t>> values;
  std::vector<Double_t> integrals;

  // The value, or a huge one where the model is not positive at an event.
  Double_t Value(const std::vector<Double_t> &a) const {
    Double_t nll = 0;
    for (size_t j = 0; j < integrals.size(); j++)
      nll += a[j] * integrals[j];
    size_t n_events = values.empty() ? 0 : values[0].size();
    for (size_t i = 0; i < n_events; i++) {
      Double_t density = 0;
      for (size_t j = 0; j < values.size(); j++)
        density += a[j] * values[j][i];
      if (density <= 0)
        return 1e300;
      nll -= TMath::Log(density);
    }
    return nll;
  }

  void GradientAndHessian(const std::vector<Double_t> &a, TVectorD &gradient,
                          TMatrixDSym &hessian) const {
    const Int_t n = (Int_t)integrals.size();
    gradient.ResizeTo(n);
    hessian.ResizeTo(n, n);
    hessian.Zero();
    for (Int_t j = 0; j < n; j++)
      gradient[j] = integrals[j];
    size_t n_events = values[0].size();
    for (size_t i = 0; i < n_events; i++) {
      Double_t density = 0;
      for (Int_t j = 0; j < n; j++)
        density += a[j] * values[j][i];
      for (Int_t j = 0; j < n; j++) {
        gradient[j] -= values[j][i] / density;
        for (Int_t k = 0; k <= j; k++)
          hessian(j, k) += values[j][i] * values[k][i] / (density * density);
      }
    }
    for (Int_t j = 0; j < n; j++)
      for (Int_t k = 0; k < j; k++)
        hessian(k, j) = hessian(j, k);
  }
};

// One unbinned fit of the in-crystal window: the continuum's coefficients
// (per keV, about CONTINUUM_PIVOT_KEV) and the Gaussian yields of the escape
// templates.
struct UnbinnedFit {
  Bool_t valid = kFALSE;
  std::vector<Double_t> coefficients, errors;
  Int_t n_continuum = 0;
};

UnbinnedFit FitUnbinned(const std::vector<Double_t> &events, Double_t lo,
                        Double_t hi, Bool_t quadratic,
                        const std::vector<EscapeTemplate> &templates,
                        const std::vector<Double_t> &start_yields,
                        const TString &tag) {
  UnbinnedFit out;
  out.n_continuum = quadratic ? 3 : 2;
  LinearUnbinnedNLL nll;
  for (Int_t power = 0; power < out.n_continuum; power++) {
    std::vector<Double_t> column(events.size());
    for (size_t i = 0; i < events.size(); i++)
      column[i] = TMath::Power(events[i] - CONTINUUM_PIVOT_KEV, power);
    nll.values.push_back(column);
    nll.integrals.push_back(
        (TMath::Power(hi - CONTINUUM_PIVOT_KEV, power + 1) -
         TMath::Power(lo - CONTINUUM_PIVOT_KEV, power + 1)) /
        (power + 1));
  }
  for (size_t t = 0; t < templates.size(); t++) {
    std::vector<Double_t> column(events.size());
    for (size_t i = 0; i < events.size(); i++)
      column[i] = templates[t](events[i]);
    nll.values.push_back(column);
    nll.integrals.push_back(templates[t].Integral(lo, hi));
  }
  const Int_t n = (Int_t)nll.integrals.size();
  std::vector<Double_t> a(n, 0.0);
  a[0] = events.size() / (hi - lo);
  for (size_t t = 0; t < templates.size(); t++)
    a[out.n_continuum + t] = start_yields[t];

  Double_t current = nll.Value(a);
  TVectorD gradient;
  TMatrixDSym hessian;
  Bool_t converged = kFALSE;
  for (Int_t iteration = 0; iteration < 100 && !converged; iteration++) {
    nll.GradientAndHessian(a, gradient, hessian);
    TMatrixDSym covariance(hessian);
    Double_t determinant = 0;
    covariance.Invert(&determinant);
    if (!(determinant > 0))
      break;
    TVectorD step = covariance * gradient;
    // Halve the Newton step until the likelihood does not rise.
    Double_t scale = 1.0;
    std::vector<Double_t> trial(a);
    Double_t value = current;
    for (Int_t halving = 0; halving < 40; halving++) {
      for (Int_t j = 0; j < n; j++)
        trial[j] = a[j] - scale * step[j];
      value = nll.Value(trial);
      if (value <= current)
        break;
      scale *= 0.5;
    }
    if (value > current)
      break;
    // Converged once every coefficient moves by a negligible part of its
    // error.
    converged = kTRUE;
    for (Int_t j = 0; j < n; j++)
      if (TMath::Abs(scale * step[j]) > 1e-6 * TMath::Sqrt(covariance(j, j)))
        converged = kFALSE;
    a = trial;
    current = value;
  }
  if (!converged) {
    std::cerr << "ERROR: the unbinned fit " << tag << " did not converge"
              << std::endl;
    return out;
  }
  nll.GradientAndHessian(a, gradient, hessian);
  TMatrixDSym covariance(hessian);
  covariance.Invert();
  for (Int_t j = 0; j < n; j++) {
    out.coefficients.push_back(a[j]);
    out.errors.push_back(TMath::Sqrt(covariance(j, j)));
  }
  out.valid = kTRUE;
  return out;
}

// The tied template: per unit of Cd K-beta Gaussian yield, Cd K-alpha and the
// Te escapes by Geant4's ratios (cd_kalpha the Cd K-alpha one).
EscapeTemplate TiedTemplate(const CaptureLine &line, Double_t sigma,
                            const RateSubShape &shape, Double_t cd_kalpha) {
  EscapeTemplate t;
  t.capture_mu = line.mu;
  t.sigma = sigma;
  t.shape = shape;
  t.series_weight[0] = cd_kalpha;
  t.series_weight[1] = 1.0;
  t.series_weight[2] = cd_kalpha * GEANT4_INSIDE_TE_KA_OVER_CD_KA;
  t.series_weight[3] = t.series_weight[2] * GEANT4_INSIDE_TE_KB_OVER_TE_KA;
  return t;
}

// The unbinned fit as the binned model's parameters, for the figure and the
// Pb check: per-bin continuum on the display binning.
EscapeFit AsBinnedFit(const UnbinnedFit &fit, TH1D *display,
                      const CaptureLine &line, Double_t sigma,
                      const RateSubShape &shape) {
  EscapeFit out;
  const Double_t bin = display->GetBinWidth(1);
  std::vector<Double_t> p(kNParameters, 0.0);
  p[kContinuum] = fit.coefficients[0] * bin;
  p[kSlope] = fit.coefficients[1] * bin;
  p[kCurvature] = fit.n_continuum > 2 ? fit.coefficients[2] * bin : 0.0;
  p[kSigma] = sigma;
  p[kTailHeightRatio] = shape.ge_tail_height_ratio;
  p[kTauOverSigma] = shape.ge_tau_over_sigma;
  p[kBinWidth] = bin;
  p[kCaptureMu] = line.mu;
  p[kCdKbetaYield] = fit.coefficients[fit.n_continuum];
  p[kCdKalphaOverCdKbeta] = CD_KA_OVER_KB;
  p[kTeKalphaOverCdKalpha] = GEANT4_INSIDE_TE_KA_OVER_CD_KA;
  p[kTeKbetaOverTeKalpha] = GEANT4_INSIDE_TE_KB_OVER_TE_KA;
  p[kPbKalpha2Mu] = ReferenceLines::PB_KA2;
  p[kPbKalpha1Mu] = ReferenceLines::PB_KA1;
  TF1 model("unbinned_on_display", EscapeWindowModel, INSIDE_WINDOW_LO_KEV,
            INSIDE_WINDOW_HI_KEV, kNParameters);
  model.SetParameters(p.data());
  WindowChi2(display, model, INSIDE_WINDOW_LO_KEV, INSIDE_WINDOW_HI_KEV,
             fit.n_continuum + 1, out.chi2, out.ndf, out.baker_cousins);
  out.series_yield[0] = SeriesYield(p.data(), 0);
  out.series_yield[1] = SeriesYield(p.data(), 1);
  out.sigma = sigma;
  out.parameters = p;
  out.valid = kTRUE;
  return out;
}

// The in-crystal measurement and its uncertainty: the unbinned fit's
// statistics (with the capture line's) and the largest shift over the
// defensible alternatives (the photopeak width, a quadratic continuum).
struct InsideResult {
  Bool_t valid = kFALSE;
  Double_t cd_kbeta = 0, stat = 0, model = 0, total = 0;
  // The escape lines' width (the 01/15 capture line's) and the alternative,
  // the Am-241 and Ba-133 photopeak width near the Ge line; keV.
  Double_t capture_sigma = 0, photopeak_sigma = 0;
};

InsideResult MeasureInside(const TString &project_root,
                           const TString &calibrated_dir,
                           const TString &fits_dir) {
  InsideResult out;
  const TString day = "20260115";
  TF1 *cal = ReadCalibration(calibrated_dir + "/calibration_function_low_" +
                             day + ".root");
  if (!cal)
    return out;
  const Double_t slope = cal->GetParameter(1);
  const Double_t bin = COMB_PERIOD_RAW_KEV * slope;
  delete cal;

  const TString am_state =
      fits_dir + "/Am_59.5keV_PostReactor_Am241_" + day + ".roofits";
  const TString ba_state =
      fits_dir + "/Ba_80.98keV_PostReactor_Ba133_" + day + ".roofits";
  PeakFitResult am, ba;
  if (!ReadSavedTail(am_state, 1, am) || !ReadSavedTail(ba_state, 2, ba))
    return out;
  RateSubShape shape = RateSubShapeFromLines(am, ba);
  // The photopeak width near the escapes, interpolated like the tail and put
  // on the calibrated scale.
  Double_t am_sigma_raw = 0, ba_sigma_raw = 0;
  if (!ReadSavedParameter(am_state, "Sigma1", am_sigma_raw) ||
      !ReadSavedParameter(ba_state, "Sigma1", ba_sigma_raw))
    return out;
  const Double_t w = (E_NEAR_GE - E_AM241) / (E_BA133_81 - E_AM241);
  const Double_t photopeak_sigma =
      slope * (am_sigma_raw + w * (ba_sigma_raw - am_sigma_raw));
  std::cout << std::fixed << "01/15 tail from the Am-241 and Ba-133 fits: "
            << "height ratio " << std::setprecision(4)
            << shape.ge_tail_height_ratio << ", tau/sigma "
            << shape.ge_tau_over_sigma << "; photopeak sigma "
            << photopeak_sigma << " keV" << std::endl;

  // The denominator is production's own capture-line yield, the quantity it
  // multiplies the fractions by.
  const TString run_path = calibrated_dir + "/" +
                           Constants::NOSHIELDBACKGROUND_5PERCENT_20260115 +
                           ".root";
  std::vector<TH1D *> hists =
      BuildCountHists(run_path, "background_0115", INSIDE_WINDOW_LO_KEV - 2.0,
                      CAPTURE_FIT_HI + 2.0, bin);
  if ((Int_t)hists.size() != N_RATESUB_PHASES)
    return out;
  CaptureLine line =
      AverageCaptureLine(hists, shape.ge_tail_height_ratio,
                         shape.ge_tau_over_sigma, kTRUE, "background_0115");
  if (!line.valid) {
    DeleteHists(hists);
    return out;
  }
  std::cout << "01/15 95.902 keV line: " << std::setprecision(4) << line.mu
            << " keV, Gaussian yield " << std::setprecision(0)
            << line.gaus_yield << " +/- " << line.gaus_yield_err << ", sigma "
            << std::setprecision(4) << line.sigma << " keV" << std::endl;

  std::vector<Double_t> events =
      ReadWindowEvents(run_path, INSIDE_WINDOW_LO_KEV, INSIDE_WINDOW_HI_KEV);
  std::cout << "01/15 escape window " << std::setprecision(1)
            << INSIDE_WINDOW_LO_KEV << "-" << INSIDE_WINDOW_HI_KEV
            << " keV: " << events.size() << " events, fitted unbinned"
            << std::endl;

  // The measurement, its alternatives, and a check that is not one: the two
  // Cd series fitted apart.
  const Int_t n_variants = 4;
  const Int_t n_alternatives = 2;
  const char *names[n_variants] = {
      "Measurement (Cd series tied, capture-line width, linear continuum):",
      "Alternative, photopeak width:", "Alternative, quadratic continuum:",
      "Check, Cd series apart:"};
  const Double_t sigmas[n_variants] = {line.sigma, photopeak_sigma, line.sigma,
                                       line.sigma};
  const Bool_t quadratic[n_variants] = {kFALSE, kFALSE, kTRUE, kFALSE};
  UnbinnedFit measurement;
  for (Int_t v = 0; v < n_variants; v++) {
    std::vector<EscapeTemplate> templates;
    std::vector<Double_t> starts;
    if (v <= n_alternatives) {
      templates.push_back(TiedTemplate(line, sigmas[v], shape, CD_KA_OVER_KB));
      starts.push_back(line.gaus_yield * ESCAPES_INSIDE.series[1]);
    } else {
      // Cd K-alpha with its Te escapes, and Cd K-beta alone.
      EscapeTemplate kalpha = TiedTemplate(line, sigmas[v], shape, 1.0);
      kalpha.series_weight[1] = 0.0;
      EscapeTemplate kbeta = TiedTemplate(line, sigmas[v], shape, 1.0);
      for (Int_t s = 0; s < 4; s++)
        kbeta.series_weight[s] = (s == 1) ? 1.0 : 0.0;
      templates.push_back(kalpha);
      templates.push_back(kbeta);
      starts.push_back(line.gaus_yield * ESCAPES_INSIDE.series[0]);
      starts.push_back(line.gaus_yield * ESCAPES_INSIDE.series[1]);
    }
    UnbinnedFit fit =
        FitUnbinned(events, INSIDE_WINDOW_LO_KEV, INSIDE_WINDOW_HI_KEV,
                    quadratic[v], templates, starts, Form("inside%d", v));
    if (!fit.valid)
      continue;
    // Cd K-beta's yield is the last template's coefficient in either case.
    Double_t kbeta_yield = fit.coefficients.back();
    Double_t kbeta_err = fit.errors.back();
    Double_t fraction = kbeta_yield / line.gaus_yield;
    Double_t stat = TMath::Sqrt(
        TMath::Power(kbeta_err / line.gaus_yield, 2) +
        TMath::Power(fraction * line.gaus_yield_err / line.gaus_yield, 2));
    std::cout << names[v] << " Cd K-beta escape fraction " << std::fixed
              << std::setprecision(4) << fraction << " +/- " << stat;
    if (v > n_alternatives)
      std::cout << ", Cd K-alpha "
                << fit.coefficients[fit.n_continuum] / line.gaus_yield
                << " +/- " << fit.errors[fit.n_continuum] / line.gaus_yield;
    std::cout << std::endl;
    if (v == 0) {
      measurement = fit;
      out.cd_kbeta = fraction;
      out.stat = stat;
    } else if (v <= n_alternatives && measurement.valid) {
      out.model = TMath::Max(out.model, TMath::Abs(fraction - out.cd_kbeta));
    }
  }
  if (!measurement.valid) {
    DeleteHists(hists);
    return out;
  }
  out.total = TMath::Sqrt(out.stat * out.stat + out.model * out.model);
  out.capture_sigma = line.sigma;
  out.photopeak_sigma = photopeak_sigma;
  out.valid = kTRUE;

  EscapeFit shown = AsBinnedFit(measurement, hists[0], line, line.sigma, shape);
  CheckPbKalpha1(hists[0], shown);
  DrawEscapeWindow(hists[0], shown, "czt_escapes_0115");
  DeleteHists(hists);
  return out;
}

// A Cd-shield background run: the Cd K-beta escape beside free Pb K-alpha
// lines, and the in-crystal and shield fraction sets compared.
struct ShieldCheck {
  Bool_t valid = kFALSE;
  Fraction cd_kbeta_te_inside, cd_kbeta_te_shield;
  Double_t delta_baker_cousins = 0; // in-crystal set minus shield set
};

ShieldCheck CheckShieldRun(const TString &calibrated_dir,
                           const TString &pair_tag, const TString &run,
                           const RateSubShape &start_shape) {
  ShieldCheck out;
  TF1 *cal =
      ReadCalibration(calibrated_dir + "/calibration_function_low_20260113_" +
                      pair_tag + ".root");
  if (!cal)
    return out;
  const Double_t bin = COMB_PERIOD_RAW_KEV * cal->GetParameter(1);
  delete cal;
  std::vector<TH1D *> hists =
      BuildCountHists(calibrated_dir + "/" + run + ".root", pair_tag,
                      SHIELD_WINDOW_LO_KEV - 2.0, CAPTURE_FIT_HI + 2.0, bin);
  if ((Int_t)hists.size() != N_RATESUB_PHASES)
    return out;

  // Every line takes the 01/15 tail; a free tail rails on these spectra.
  const RateSubShape &shape = start_shape;
  CaptureLine line =
      AverageCaptureLine(hists, shape.ge_tail_height_ratio,
                         shape.ge_tau_over_sigma, kTRUE, pair_tag);
  if (!line.valid) {
    DeleteHists(hists);
    return out;
  }
  std::cout << pair_tag << " background run: 95.902 keV line Gaussian yield "
            << std::setprecision(0) << line.gaus_yield << " +/- "
            << line.gaus_yield_err << std::endl;
  WindowFitConfig tied;
  tied.lo = SHIELD_WINDOW_LO_KEV;
  tied.hi = SHIELD_WINDOW_HI_KEV;
  tied.with_pb = kTRUE;

  WindowFitConfig te_shield = tied;
  te_shield.te_kalpha_over_cd_kalpha =
      ESCAPES_OUTSIDE.series[2] / ESCAPES_OUTSIDE.series[0];
  te_shield.te_kbeta_over_te_kalpha =
      ESCAPES_OUTSIDE.series[3] / ESCAPES_OUTSIDE.series[2];
  WindowFitConfig inside_set = tied;
  inside_set.mode = kFixedSet;
  inside_set.fixed_set = &ESCAPES_INSIDE;
  WindowFitConfig shield_set = tied;
  shield_set.mode = kFixedSet;
  shield_set.fixed_set = &ESCAPES_OUTSIDE;

  std::vector<EscapeFit> fits_te_inside =
      FitPhases(hists, line, shape, tied, pair_tag + "_te_inside");
  std::vector<EscapeFit> fits_te_shield =
      FitPhases(hists, line, shape, te_shield, pair_tag + "_te_shield");
  std::vector<EscapeFit> fits_inside =
      FitPhases(hists, line, shape, inside_set, pair_tag + "_inside_set");
  std::vector<EscapeFit> fits_shield =
      FitPhases(hists, line, shape, shield_set, pair_tag + "_shield_set");
  out.cd_kbeta_te_inside = SeriesFraction(fits_te_inside, line, 1);
  out.cd_kbeta_te_shield = SeriesFraction(fits_te_shield, line, 1);
  out.delta_baker_cousins =
      MeanBakerCousins(fits_inside) - MeanBakerCousins(fits_shield);
  out.valid = out.cd_kbeta_te_inside.valid && out.cd_kbeta_te_shield.valid;
  PrintFraction(pair_tag + ", Te escapes the in-crystal way:", "Cd K-beta",
                out.cd_kbeta_te_inside);
  PrintFraction(pair_tag + ", Te escapes the shield way:", "Cd K-beta",
                out.cd_kbeta_te_shield);
  std::cout << pair_tag << ": Baker-Cousins chi2 of the in-crystal set minus "
            << "the shield set's " << std::setprecision(1)
            << out.delta_baker_cousins << " (positive favours the shield)"
            << std::endl;
  DeleteHists(hists);
  return out;
}

} // namespace

void CaptureEscapeFraction() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/czt/figures",
                                project_root + "/root_files");
  const TString calibrated_dir =
      project_root + "/root_files/calibrated_" + PRODUCTION_TAG;
  const TString fits_dir =
      project_root + "/plots/czt/chain/" + PRODUCTION_TAG + "/fits";

  InsideResult inside = MeasureInside(project_root, calibrated_dir, fits_dir);
  if (!inside.valid)
    return;
  const Geant4Fit geant4 =
      FitGeant4InCrystal(project_root, inside.capture_sigma);
  const Double_t ratio = CD_KA_OVER_KB;
  const Double_t kbeta = inside.cd_kbeta;
  const Double_t measured[4] = {ratio * kbeta, kbeta,
                                ratio * kbeta * GEANT4_INSIDE_TE_KA_OVER_CD_KA,
                                ratio * kbeta * GEANT4_INSIDE_TE_KA_OVER_CD_KA *
                                    GEANT4_INSIDE_TE_KB_OVER_TE_KA};
  std::cout << std::fixed << std::setprecision(4)
            << "In-crystal fractions: Cd K-alpha " << measured[0]
            << ", Cd K-beta " << measured[1] << ", Te K-alpha " << measured[2]
            << ", Te K-beta " << measured[3] << "; Cd K-beta +/- "
            << inside.stat << " (stat) +/- " << inside.model
            << " (alternatives) = +/- " << inside.total << " ("
            << std::setprecision(1) << 100.0 * inside.total / kbeta << "%)"
            << std::endl;
  // The constants production uses, at the precision they are written to.
  for (Int_t s = 0; s < 4; s++)
    if (TMath::Nint(1e5 * measured[s]) !=
        TMath::Nint(1e5 * ESCAPES_INSIDE.series[s]))
      std::cout << "WARNING: ESCAPES_INSIDE.series[" << s << "] is "
                << std::setprecision(5) << ESCAPES_INSIDE.series[s]
                << ", the measurement " << measured[s] << std::endl;
  if (TMath::Nint(1e3 * inside.total / kbeta) !=
      TMath::Nint(1e3 * ESCAPES_INSIDE_REL_ERR))
    std::cout << "WARNING: ESCAPES_INSIDE_REL_ERR is " << std::setprecision(3)
              << ESCAPES_INSIDE_REL_ERR << ", the measurement "
              << inside.total / kbeta << std::endl;

  // The Cd-shield fits take the 01/15 tail.
  PeakFitResult am, ba;
  RateSubShape start_shape;
  if (ReadSavedTail(fits_dir + "/Am_59.5keV_PostReactor_Am241_20260115.roofits",
                    1, am) &&
      ReadSavedTail(
          fits_dir + "/Ba_80.98keV_PostReactor_Ba133_20260115.roofits", 2, ba))
    start_shape = RateSubShapeFromLines(am, ba);
  ShieldCheck cd10 = CheckShieldRun(
      calibrated_dir, "CdShield10",
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113, start_shape);
  ShieldCheck cd25 = CheckShieldRun(
      calibrated_dir, "CdShield25",
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113, start_shape);

  const TString result_path = project_root + "/results/escape_fractions.result";
  std::ofstream out(result_path.Data());
  out << "# Written by czt/CaptureEscapeFraction.cpp. Fractions of the 95.902 "
         "keV line's Gaussian yield."
      << std::endl;
  out << "# In-crystal set, measured in the 01/15 beam background run "
         "(ESCAPES_INSIDE takes it)."
      << std::endl;
  out << std::fixed << std::setprecision(6);
  out << "inside_cd_kalpha " << measured[0] << std::endl;
  out << "inside_cd_kbeta " << measured[1] << std::endl;
  out << "inside_te_kalpha " << measured[2] << std::endl;
  out << "inside_te_kbeta " << measured[3] << std::endl;
  out << "inside_cd_kbeta_err_stat " << inside.stat << std::endl;
  out << "inside_cd_kbeta_err_alternatives " << inside.model << std::endl;
  out << "inside_cd_kbeta_err " << inside.total << std::endl;
  out << "inside_relative_err " << inside.total / kbeta << std::endl;
  out << "# The escape lines' Gaussian width, the 01/15 capture line's, and "
         "the alternative, the Am-241 and Ba-133 photopeak width near the Ge "
         "line (Escape96Plots.cpp smears the Geant4 spectrum with the first). "
         "keV."
      << std::endl;
  out << "capture_sigma " << inside.capture_sigma << std::endl;
  out << "photopeak_sigma " << inside.photopeak_sigma << std::endl;
  if (geant4.valid) {
    out << "# Production's escape model on the Geant4 in-crystal run smeared "
           "with capture_sigma, 60-80 keV: Pearson chi2 and ndf, and the Cd "
           "K-alpha / K-beta ratio with the two series fitted apart (the fits "
           "hold Geant4's counted ratio)."
        << std::endl;
    out << std::setprecision(2) << "geant4_fit_chi2 " << geant4.chi2
        << std::endl;
    out << "geant4_fit_ndf " << geant4.ndf << std::endl;
    out << std::setprecision(4) << "geant4_cd_ka_over_kb_apart "
        << geant4.cd_ka_over_kb_apart << std::endl;
    out << "geant4_cd_ka_over_kb_apart_err " << geant4.cd_ka_over_kb_apart_err
        << std::endl;
    out << "geant4_cd_ka_over_kb_held " << CD_KA_OVER_KB << std::endl;
    out << std::setprecision(6);
  }
  out << "# The escape lines, 95.902 keV minus the SRD 128 X-ray energies, two "
         "per series (Cd K-alpha, Cd K-beta, Te K-alpha, Te K-beta). keV."
      << std::endl;
  const char *series_keys[4] = {"cd_kalpha", "cd_kbeta", "te_kalpha",
                                "te_kbeta"};
  Int_t line_in_series[4] = {0, 0, 0, 0};
  out << std::setprecision(3);
  for (Int_t i = 0; i < N_ESCAPE_LINES; i++) {
    Int_t s = ESCAPE_SERIES[i];
    out << "escape_line_" << series_keys[s] << "_" << ++line_in_series[s] << " "
        << ReferenceLines::CD113_CAPTURE - ESCAPE_XRAY_KEV[i] << std::endl;
  }
  out << std::setprecision(6);
  out << "# Cd-shield background runs of 01/13: the Cd K-beta escape beside "
         "Pb K-alpha, and the Baker-Cousins chi2 of the in-crystal set minus "
         "the shield set's."
      << std::endl;
  const ShieldCheck *checks[2] = {&cd10, &cd25};
  const char *keys[2] = {"cd10", "cd25"};
  for (Int_t i = 0; i < 2; i++) {
    if (!checks[i]->valid)
      continue;
    const Fraction &f = checks[i]->cd_kbeta_te_inside;
    out << keys[i] << "_cd_kbeta " << f.value << std::endl;
    out << keys[i] << "_cd_kbeta_err "
        << TMath::Sqrt(f.stat * f.stat + f.phase_rms * f.phase_rms)
        << std::endl;
    out << keys[i] << "_cd_kbeta_te_shield "
        << checks[i]->cd_kbeta_te_shield.value << std::endl;
    out << keys[i] << "_delta_chi2_inside_minus_shield " << std::setprecision(2)
        << checks[i]->delta_baker_cousins << std::setprecision(6) << std::endl;
  }
  std::cout << "Wrote " << result_path << std::endl;
}
