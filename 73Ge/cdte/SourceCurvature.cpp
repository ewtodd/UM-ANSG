// The relative curvature kappa of the CdTe's energy scale, measured in the
// Ba-133 source runs per rise-time slice with the beam fit's scale and
// lineshape (CdTeAnalysis.hpp): every line at its true energy on
// x = a + b (E - 70 + kappa (E - 70)^2), every line in PeakShape, the shape
// shared and its width linear in energy. kappa carries over to the beam fits,
// whose gain differs. The lines, over 55.5-60.5 and 70-88.5 keV (the
// charge-loss shelf and the Pb K-beta escapes between them are left out):
//   - Ba-133 81.00 keV and 79.61 keV, the yields at ENSDF's 32.9 : 2.65;
//   - their Cd K-alpha and K-beta escapes, at the line's energy minus the
//     X-ray's, with one free escape fraction and K-beta tied to K-alpha at
//     the beam fit's nominal ratio (the K-beta escapes sit just below the
//     low window);
//   - Pb K-alpha2/K-alpha1 and the Pb K-beta groups from the shielding.
// Ge (68.75 keV) lies inside that span, between the 57.8 keV escape and
// Pb K-alpha2.
//
// The slices are fitted without and with the rise-gain correction: without
// it, kappa changes from slice to slice beyond its errors, because the gain
// changes inside a slice and the lines do not share one rise distribution
// (RiseTime.cpp shows both per line). Am-241 59.54 keV from the Am source
// runs is placed on each slice's scale as a check; its runs and its gain
// against rise differ from the Ba runs', so it does not test the scale to
// tens of eV.
//
// Usage: root -b -q 'SourceCurvature.cpp+'
#include "CdTeAnalysis.hpp"

namespace CdTeCurvature {

using namespace CdTe;

const Double_t LOW_LO = 55.5, LOW_HI = 60.5, HIGH_LO = 70.0, HIGH_HI = 88.5;

enum {
  S_A,
  S_B,
  S_KAPPA, // x = a + b (E - 70 + kappa (E - 70)^2)
  S_SIG0,
  S_SIG1, // sigma(E) = sig0 + sig1 (E - 70)
  S_STEP, // the shape every line shares, in PeakShape's order
  S_LOW_EXP,
  S_LOW_EXP_DECAY,
  S_HIGH_EXP,
  S_HIGH_EXP_DECAY,
  S_LOW_C0,
  S_LOW_C1, // continuum in the low window, linear in (x - 58)
  S_HIGH_C0,
  S_HIGH_C1, // continuum in the high window, linear in (x - 79)
  S_BA81_A,
  S_BA80_R, // 79.61 keV yield over 81.00 keV's
  S_ESC_KA, // Cd K-alpha escape fraction of each Ba line
  S_KA1_A,
  S_KA2_R,
  S_KB1_A,
  S_KB2_R,
  N_S
};

Double_t Position(Double_t e, const Double_t *p) {
  return CurvedScale(e, p[S_A], p[S_B], p[S_KAPPA]);
}

Double_t Line(Double_t x, Double_t e, const Double_t *p) {
  Double_t sigma = p[S_SIG0] + p[S_SIG1] * (e - 70.0);
  return PeakShape(x, Position(e, p), sigma, &p[S_STEP]);
}

Double_t BaLines(Double_t x, const Double_t *p) {
  const Double_t parent_e[2] = {ReferenceLines::BA133_81,
                                ReferenceLines::BA133_79};
  const Double_t parent_a[2] = {p[S_BA81_A], p[S_BA81_A] * p[S_BA80_R]};
  Double_t v = 0;
  for (Int_t k = 0; k < 2; k++) {
    Double_t lines = Line(x, parent_e[k], p);
    for (Int_t i = 0; i < ReferenceLines::N_ESCAPE_LINES; i++) {
      Int_t series = ReferenceLines::ESCAPE_SERIES[i];
      if (series > 1)
        continue;
      Double_t fraction =
          series == 0 ? p[S_ESC_KA] : KB_KA_NOMINAL * p[S_ESC_KA];
      lines += fraction * GEANT4_SERIES_SHARE[i] *
               Line(x, parent_e[k] - ReferenceLines::ESCAPE_XRAY_KEV[i], p);
    }
    v += parent_a[k] * lines;
  }
  return v;
}

Double_t PbLines(Double_t x, const Double_t *p) {
  Double_t v = p[S_KA1_A] * (Line(x, ReferenceLines::PB_KA1, p) +
                             p[S_KA2_R] * Line(x, ReferenceLines::PB_KA2, p));
  Double_t kb1 = 0, kb2 = 0;
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB1; k++)
    kb1 += ReferenceLines::PB_KB1[k].weight *
           Line(x, ReferenceLines::PB_KB1[k].energy, p);
  for (Int_t k = 0; k < ReferenceLines::N_PB_KB2; k++)
    kb2 += ReferenceLines::PB_KB2[k].weight *
           Line(x, ReferenceLines::PB_KB2[k].energy, p);
  return v + p[S_KB1_A] * (kb1 + p[S_KB2_R] * kb2);
}

Double_t Continuum(Double_t x, const Double_t *p) {
  if (x < 0.5 * (LOW_HI + HIGH_LO))
    return p[S_LOW_C0] + p[S_LOW_C1] * (x - 58.0);
  return p[S_HIGH_C0] + p[S_HIGH_C1] * (x - 79.0);
}

Double_t Model(Double_t *x, Double_t *p) {
  if (x[0] > LOW_HI && x[0] < HIGH_LO) {
    TF1::RejectPoint();
    return 0;
  }
  return Continuum(x[0], p) + BaLines(x[0], p) + PbLines(x[0], p);
}

Double_t DrawSourceContinuum(Double_t x, const Double_t *p) {
  return Continuum(x, p);
}
Double_t DrawSourceBa(Double_t x, const Double_t *p) {
  return Continuum(x, p) + BaLines(x, p);
}
Double_t DrawSourcePb(Double_t x, const Double_t *p) {
  return Continuum(x, p) + PbLines(x, p);
}

struct SliceCurvature {
  Double_t kappa = 0, kappa_err = 0, b = 0, chi2_ndf = 0;
  Double_t am_resid = 0, am_resid_err = 0;
  Int_t status = -1, cov_status = -1;
};

// Staged as the beam fit: the scale fixed at the 81.00 keV peak and the high
// tail off, then the scale released, then the curvature, the Pb ratios, the
// width's slope and the high tail.
TF1 *FitSource(TH1D *h, const TString &tag, Int_t &status, Int_t &cov_status) {
  TF1 *f = new TF1("source_" + tag, Model, LOW_LO, HIGH_HI, N_S);
  f->SetNpx(3000);
  h->GetXaxis()->SetRangeUser(80.0, 82.0);
  Double_t ba81_mode = h->GetBinCenter(h->GetMaximumBin());
  Double_t ba81_height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  Double_t a0 = ba81_mode - (ReferenceLines::BA133_81 - 70.0);
  const Double_t sigma0 = 0.35;
  Double_t ba81 = ba81_height * TMath::Sqrt(2.0 * TMath::Pi()) * sigma0 / BIN;
  f->SetParameter(S_SIG0, sigma0);
  f->SetParLimits(S_SIG0, 0.1, 1.5);
  StartShape(f, S_STEP);
  f->SetParameter(S_LOW_C0, h->GetBinContent(h->FindBin(60.3)));
  f->SetParameter(S_LOW_C1, 0);
  f->SetParameter(S_HIGH_C0, h->GetBinContent(h->FindBin(71.0)));
  f->SetParameter(S_HIGH_C1, 0);
  f->SetParameter(S_BA81_A, ba81);
  f->SetParLimits(S_BA81_A, 0, 1e10);
  f->SetParameter(S_KA1_A, 0.02 * ba81);
  f->SetParLimits(S_KA1_A, 0, 1e9);
  f->SetParameter(S_KB1_A, 0.005 * ba81);
  f->SetParLimits(S_KB1_A, 0, 1e9);
  f->FixParameter(S_A, a0);
  f->FixParameter(S_B, 1.0);
  f->FixParameter(S_KAPPA, 0.0);
  f->FixParameter(S_SIG1, 0.0);
  f->FixParameter(S_BA80_R, 2.65 / 32.9);
  f->FixParameter(S_ESC_KA, 0.1);
  f->FixParameter(S_KA2_R, 0.6);
  f->FixParameter(S_KB2_R, 0.3);
  h->Fit(f, "QLRN0");

  f->ReleaseParameter(S_A);
  f->SetParLimits(S_A, a0 - 1.0, a0 + 1.0);
  f->ReleaseParameter(S_B);
  f->SetParLimits(S_B, 0.95, 1.05);
  f->ReleaseParameter(S_ESC_KA);
  f->SetParLimits(S_ESC_KA, 0.0, 0.5);
  h->Fit(f, "QLRN0");

  f->ReleaseParameter(S_KA2_R);
  f->SetParLimits(S_KA2_R, 0.3, 1.2);
  f->ReleaseParameter(S_KB2_R);
  f->SetParLimits(S_KB2_R, 0.1, 0.6);
  f->ReleaseParameter(S_SIG1);
  f->SetParLimits(S_SIG1, -0.05, 0.05);
  ReleaseHighTail(f, S_STEP);
  f->ReleaseParameter(S_KAPPA);
  f->SetParLimits(S_KAPPA, -0.01, 0.01);
  h->Fit(f, "QLRN0");
  TFitResultPtr result = h->Fit(f, "QLRN0MS");
  status = result->Status();
  cov_status = result->CovMatrixStatus();
  return f;
}

// kappa's error from the likelihood profile: the fit redone with kappa held
// at its value +- step, the deviance's rise on each side read as a parabola.
// The fit's covariance is forced positive-definite in some slices, so its
// error is not used.
Double_t ProfileError(TH1D *h, TF1 *f, Double_t step) {
  const Int_t n_par = f->GetNpar();
  std::vector<Double_t> best(f->GetParameters(), f->GetParameters() + n_par);
  const Double_t deviance = f->GetChisquare();
  Double_t limit_lo = 0, limit_hi = 0;
  f->GetParLimits(S_KAPPA, limit_lo, limit_hi);
  Double_t side_err[2] = {0, 0};
  for (Int_t side = 0; side < 2; side++) {
    f->SetParameters(&best[0]);
    f->FixParameter(S_KAPPA, best[S_KAPPA] + (side == 0 ? -step : step));
    h->Fit(f, "QLRN0");
    Double_t rise = f->GetChisquare() - deviance;
    side_err[side] = rise > 0 ? step / TMath::Sqrt(rise) : -1;
  }
  f->SetParameters(&best[0]);
  f->ReleaseParameter(S_KAPPA);
  f->SetParLimits(S_KAPPA, limit_lo, limit_hi);
  if (side_err[0] < 0 || side_err[1] < 0)
    return -1;
  return 0.5 * (side_err[0] + side_err[1]);
}

} // namespace CdTeCurvature

void SourceCurvature() {
  using namespace CdTeCurvature;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/source_curvature",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);

  const Int_t n = N_RISE_SLICES;
  const Double_t *lo_ns = RISE_SLICE_LO;
  const Double_t *hi_ns = RISE_SLICE_HI;
  for (Int_t pass = 0; pass < 2; pass++) {
    g_rise_gain_correction = (pass == 1);
    TString pass_name = g_rise_gain_correction ? "corrected" : "as measured";
    std::vector<SliceCurvature> out(n);
    for (Int_t s = 0; s < n; s++) {
      TString tag = Form("%.0f_%.0f", lo_ns[s], hi_ns[s]);
      TH1D *ba = FillRise(BaSourceRuns(), "ba_" + tag, lo_ns[s], hi_ns[s]);
      TH1D *am = FillRise(AmSourceRuns(), "am_" + tag, lo_ns[s], hi_ns[s]);
      SliceCurvature &r = out[s];
      TF1 *f = FitSource(ba, tag, r.status, r.cov_status);
      const Double_t *p = f->GetParameters();
      r.b = p[S_B];
      r.kappa = p[S_KAPPA];
      r.chi2_ndf = f->GetChisquare() / f->GetNDF();
      Double_t hesse_err = f->GetParError(S_KAPPA);
      r.kappa_err = ProfileError(ba, f, TMath::Max(hesse_err, 5e-5));
      if (r.kappa_err < 0)
        r.kappa_err = hesse_err;
      PeakFit am_peak =
          FitPeak(am, ReferenceLines::AM241, 56.0, 61.5, "am_" + tag, kFALSE);
      r.am_resid = CurvedScaleInverse(am_peak.mu, p[S_A], p[S_B], p[S_KAPPA]) -
                   ReferenceLines::AM241;
      r.am_resid_err = am_peak.mu_err / p[S_B];
      std::cout << "SOURCE-CURV [" << pass_name << "] rise " << tag
                << " ns: kappa " << std::scientific << std::setprecision(2)
                << r.kappa << " +- " << r.kappa_err << " (Hesse " << hesse_err
                << ")" << std::fixed << " /keV  b " << std::setprecision(4)
                << r.b << "  chi2/ndf " << std::setprecision(2) << r.chi2_ndf
                << "  status " << r.status << " cov " << r.cov_status
                << "  esc Ka " << std::setprecision(3) << p[S_ESC_KA]
                << "  | Am 59.54 on this scale " << std::showpos
                << std::setprecision(1) << 1000 * r.am_resid << std::noshowpos
                << " +- " << 1000 * r.am_resid_err << " eV" << std::endl;
      std::cout << "SOURCE-CURV [" << pass_name << "]   shape: sigma "
                << std::setprecision(3) << p[S_SIG0] << " + " << p[S_SIG1]
                << " (E - 70)  step " << p[S_STEP] << "  low exp "
                << p[S_LOW_EXP] << " x " << p[S_LOW_EXP_DECAY]
                << " sigma  high exp " << p[S_HIGH_EXP] << " x "
                << p[S_HIGH_EXP_DECAY] << " sigma" << std::endl;
      if (g_rise_gain_correction) {
        const Double_t window_lo[2] = {LOW_LO, HIGH_LO};
        const Double_t window_hi[2] = {LOW_HI, HIGH_HI};
        const char *window_name[2] = {"low", "high"};
        for (Int_t w = 0; w < 2; w++) {
          std::vector<TGraph *> components;
          components.push_back(Component(window_lo[w], window_hi[w], kGray + 2,
                                         2, DrawSourceContinuum, p));
          components.push_back(Component(window_lo[w], window_hi[w], kBlue + 1,
                                         1, DrawSourceBa, p));
          components.push_back(Component(window_lo[w], window_hi[w],
                                         kOrange + 7, 1, DrawSourcePb, p));
          Plot(ba, f, components, window_lo[w], window_hi[w],
               Form("psi_cdte_source_curvature_%s_%s", tag.Data(),
                    window_name[w]),
               Form("Ba-133 source runs, rise %.1f-%.1f ns", lo_ns[s],
                    hi_ns[s]));
        }
      }
      delete f;
      delete ba;
      delete am;
    }
    Double_t sw = 0, swx = 0;
    for (Int_t s = 0; s < n; s++) {
      sw += 1.0 / (out[s].kappa_err * out[s].kappa_err);
      swx += out[s].kappa / (out[s].kappa_err * out[s].kappa_err);
    }
    Double_t mean = swx / sw, mean_err = 1.0 / TMath::Sqrt(sw), chi2 = 0;
    for (Int_t s = 0; s < n; s++)
      chi2 += TMath::Power((out[s].kappa - mean) / out[s].kappa_err, 2);
    std::cout << "SOURCE-CURV [" << pass_name << "] weighted mean kappa "
              << std::scientific << std::setprecision(2) << mean << " +- "
              << mean_err << std::fixed << " /keV, chi2/ndf "
              << std::setprecision(1) << chi2 << "/" << n - 1 << std::endl;
  }
}
