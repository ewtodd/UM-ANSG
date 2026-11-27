// Where the CdTe testbed's closure bias comes from (CdTeTestbed.cpp: the
// nominal fit reads Ge about 5 eV low on pseudo-spectra drawn from itself).
// In one cell (a day's Ge runs in one rise slice) the nominal fit to the data
// is the truth; each Poisson pseudo-spectrum is fitted five ways, each staged
// as FitInRun:
//   staged           the production fit;
//   unbounded escape the Cd K-beta escape's yield allowed below zero: does
//                    its bound at zero push Ge?
//   escape held      the escape's yield held at the truth's;
//   tails held       the shared tails held at the truth's: do the free tails
//                    trade with Ge?
//   both held        the tails and the escape's yield held at the truth's.
// (The Ge, Pb K-alpha1 and K-beta1 yields' [0, 1e9] limits leave Minuit's
// errors on them meaningless, but removing them changes no fit.)
// One TRIAL line per pseudo-spectrum and fit: Ge (in-run and on the source
// scale) less the truth's, the shape and ratio parameters, the fit's status
// and the parameters that end at a limit.
//
// Usage: root -b -q 'TestbedBias.cpp+(<cell 0-7>, <trials>)'; cell = day
// index (GE_DAYS) x N_RISE_SLICES + slice.
#include "CdTeAnalysis.hpp"

#include <TRandom3.h>

namespace CdTeTestbedBias {

using namespace CdTe;

enum {
  F_STAGED,
  F_UNBOUNDED_ESCAPE,
  F_ESCAPE_HELD,
  F_TAILS_HELD,
  F_BOTH_HELD,
  N_FITS
};
const char *const FIT_NAMES[N_FITS] = {
    "staged", "unbounded_escape", "escape_held", "tails_held", "both_held"};

const Int_t N_WATCHED = 10;
const Int_t WATCHED[N_WATCHED] = {
    Q_SIG0,           Q_SIG1,  Q_LOW_EXP, Q_LOW_EXP_DECAY, Q_HIGH_EXP,
    Q_HIGH_EXP_DECAY, Q_ESC_A, Q_KA2_R,   Q_KB2_R,         Q_BKG};
const char *const WATCHED_NAMES[N_WATCHED] = {
    "sig0",       "sig1", "low_exp", "low_decay", "high_exp",
    "high_decay", "esc",  "ka2_r",   "kb2_r",     "bkg"};

void NameParameters(TF1 *f) {
  for (Int_t k = 0; k < N_Q; k++)
    f->SetParName(k, Form("p%d", k));
  for (Int_t k = 0; k < N_WATCHED; k++)
    f->SetParName(WATCHED[k], WATCHED_NAMES[k]);
  f->SetParName(Q_GE_E, "ge_e");
  f->SetParName(Q_GE_A, "ge_a");
  f->SetParName(Q_KA1_A, "ka1_a");
  f->SetParName(Q_KB1_A, "kb1_a");
  f->SetParName(Q_A, "a");
  f->SetParName(Q_B, "b");
}

TH1D *Pseudo(TH1D *data, TF1 *truth, TRandom3 &random, const TString &name) {
  TH1D *h = (TH1D *)data->Clone(name);
  h->SetDirectory(nullptr);
  for (Int_t b = h->FindBin(g_lo); b <= h->FindBin(g_hi); b++) {
    Double_t expected = TMath::Max(truth->Eval(h->GetBinCenter(b)), 0.0);
    Double_t counts = random.Poisson(expected);
    h->SetBinContent(b, counts);
    h->SetBinError(b, TMath::Sqrt(TMath::Max(counts, 1.0)));
  }
  return h;
}

// FitInRun's stages with the study's changes; status is the last fit's.
TF1 *FitVariant(TH1D *h, const TString &tag, Int_t fit, TF1 *truth,
                Int_t &status) {
  TF1 *f = new TF1("variant_" + tag, InRunModel, g_lo, g_hi, N_Q);
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
  if (fit != F_UNBOUNDED_ESCAPE)
    f->SetParLimits(Q_ESC_A, 0, f->GetParameter(Q_KA1_A));
  ReleaseInRunShape(f);
  if (fit == F_TAILS_HELD || fit == F_BOTH_HELD)
    for (Int_t k = Q_LOW_EXP; k <= Q_HIGH_EXP_DECAY; k++)
      f->FixParameter(k, truth->GetParameter(k));
  if (fit == F_ESCAPE_HELD || fit == F_BOTH_HELD)
    f->FixParameter(Q_ESC_A, truth->GetParameter(Q_ESC_A));
  h->Fit(f, "QLRN0");
  TFitResultPtr result = h->Fit(f, "QLRN0MS");
  status = result->Status();
  g_inrun_cov.Zero();
  if (result.Get() != nullptr && result->CovMatrixStatus() > 0) {
    g_inrun_cov = result->GetCovarianceMatrix();
  } else {
    for (Int_t k = 0; k < N_Q; k++)
      g_inrun_cov(k, k) = f->GetParError(k) * f->GetParError(k);
  }
  NameParameters(f);
  return f;
}

// The parameters that end at a limit, to 1e-4 of their range (at most 0.01).
TString AtLimit(TF1 *f) {
  TString out;
  for (Int_t k = 0; k < N_Q; k++) {
    Double_t lo = 0, hi = 0;
    f->GetParLimits(k, lo, hi);
    if (!(hi > lo))
      continue;
    Double_t p = f->GetParameter(k);
    Double_t tolerance = 1e-4 * TMath::Min(hi - lo, 100.0);
    if (p - lo < tolerance || hi - p < tolerance)
      out += Form("%s%s", out.Length() ? "," : "", f->GetParName(k));
  }
  return out.Length() ? out : TString("-");
}

void PrintTrial(Int_t cell, Int_t trial, Int_t fit, TF1 *f, Int_t status,
                Double_t ge_truth, Double_t ge_truth_in_run,
                const SourceLine &line) {
  Transferred t = GainTransfer(line, f);
  std::cout << "TRIAL " << cell << " " << trial << " " << FIT_NAMES[fit]
            << std::fixed << std::setprecision(2) << " d_ge "
            << 1000 * (t.e - ge_truth) << " d_ge_in_run "
            << 1000 * (f->GetParameter(Q_GE_E) - ge_truth_in_run) << " err "
            << 1000 * t.err;
  for (Int_t k = 0; k < N_WATCHED; k++)
    std::cout << " " << WATCHED_NAMES[k] << " " << std::setprecision(5)
              << f->GetParameter(WATCHED[k]);
  std::cout << " ge_a " << std::setprecision(1) << f->GetParameter(Q_GE_A)
            << " ka1_a " << f->GetParameter(Q_KA1_A) << " deviance "
            << std::setprecision(3) << f->GetChisquare() << " status " << status
            << " at_limit " << AtLimit(f) << std::endl;
}

} // namespace CdTeTestbedBias

void TestbedBias(Int_t cell = 0, Int_t trials = 100) {
  using namespace CdTeTestbedBias;
  g_lo = 62.0;
  g_hi = 89.0;
  ResetKbGroups();
  SetEscapeTie(KB_KA_NOMINAL, TE_CDKB);
  g_tails = TAILS_FREE;
  g_kappa = KAPPA_SOURCE;
  const Int_t d = cell / N_RISE_SLICES, s = cell % N_RISE_SLICES;
  TRandom3 random(9100 + cell);

  SourceLine line = FitSourceLine(GE_DAYS[d], s);
  TH1D *h = FillRise(OfDay(GeRuns(), GE_DAYS[d]), Form("data_%d", cell),
                     RISE_SLICE_LO[s], RISE_SLICE_HI[s]);
  PeakFit capture = FitCaptureLine(h, Form("capture_%d", cell));
  g_capture_meas = capture.mu;
  TF1 *truth = FitInRun(h, Form("truth_%d", cell));
  NameParameters(truth);
  const Double_t ge_truth = GainTransfer(line, truth).e;
  const Double_t ge_truth_in_run = truth->GetParameter(Q_GE_E);
  std::cout << "TRUTH cell " << cell << " 10/" << GE_DAYS[d] << " slice "
            << std::fixed << std::setprecision(2) << RISE_SLICE_LO[s] << "-"
            << RISE_SLICE_HI[s] << " ns: Ge " << std::setprecision(4)
            << ge_truth << " (in-run " << ge_truth_in_run << ")" << std::endl;
  PrintTrial(cell, -1, F_STAGED, truth, 0, ge_truth, ge_truth_in_run, line);

  for (Int_t t = 0; t < trials; t++) {
    TH1D *pseudo = Pseudo(h, truth, random, Form("pseudo_%d_%d", cell, t));
    for (Int_t fit = 0; fit < N_FITS; fit++) {
      Int_t status = 0;
      TF1 *f = FitVariant(pseudo, Form("%s_%d_%d", FIT_NAMES[fit], cell, t),
                          fit, truth, status);
      PrintTrial(cell, t, fit, f, status, ge_truth, ge_truth_in_run, line);
      delete f;
    }
    delete pseudo;
  }
}
