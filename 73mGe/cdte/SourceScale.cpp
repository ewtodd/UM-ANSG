// Ge on the source runs' scale, as the paper draft's CdTe calibration does it:
// in each rise slice, a straight line through the Am-241 59.5409 keV and
// Ba-133 80.9979 keV lines of the source runs (heights over the rise gain,
// like the beam runs'), applied to the positions the nominal in-run fit
// (GeBudget.cpp) finds for Ge and the Pb K-alpha lines in the Ge beam runs.
// Am lies below Ge, so on this scale Ge is interpolated where the in-run
// scale, fixed by the Pb lines above it, extrapolates. Pb K-alpha on the
// source scale shows how well the source runs' scale carries over to the
// beam runs.
//
// Usage: root -b -q 'SourceScale.cpp+'
#include "CdTeAnalysis.hpp"

namespace CdTeSourceScale {

using namespace CdTe;

// GeBudget.cpp's nominal escape tie.
const Double_t KB_KA = 0.275, TE_CDKB = 0.352;
const Double_t BA133_81 = ReferenceLines::BA133_81;

// keV of spectrum position x on the line through (x_am, Am) and (x_ba, Ba),
// with its error from the two positions' errors.
Double_t SourceEnergy(Double_t x, Double_t x_am, Double_t x_am_err,
                      Double_t x_ba, Double_t x_ba_err, Double_t &err) {
  const Double_t e_am = ReferenceLines::AM241;
  Double_t slope = (BA133_81 - e_am) / (x_ba - x_am);
  Double_t e = e_am + slope * (x - x_am);
  // d e / d x_am and d e / d x_ba
  Double_t u = (x - x_am) / (x_ba - x_am);
  Double_t d_am = slope * (u - 1.0), d_ba = -slope * u;
  err = TMath::Sqrt(d_am * d_am * x_am_err * x_am_err +
                    d_ba * d_ba * x_ba_err * x_ba_err);
  return e;
}

} // namespace CdTeSourceScale

void SourceScale() {
  using namespace CdTeSourceScale;
  g_lo = 62.0;
  g_hi = 89.0;
  g_kappa = KAPPA_SOURCE;
  ResetKbGroups();
  SetEscapeTie(KB_KA, TE_CDKB);
  Double_t sw = 0, swx = 0, sw_inrun = 0, swx_inrun = 0;
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    TString tag = Form("%d", s);
    TH1D *am = FillRise(AmSourceRuns(), "am_" + tag, RISE_SLICE_LO[s],
                        RISE_SLICE_HI[s]);
    TH1D *ba = FillRise(BaSourceRuns(), "ba_" + tag, RISE_SLICE_LO[s],
                        RISE_SLICE_HI[s]);
    PeakFit am_peak =
        FitPeak(am, ReferenceLines::AM241, 56.0, 61.5, "am_" + tag, kFALSE);
    PeakFit ba_peak = FitPeak(ba, BA133_81, 80.2, 82.5, "ba_" + tag, kFALSE);

    TH1D *h =
        FillRise(GeRuns(), "ge_" + tag, RISE_SLICE_LO[s], RISE_SLICE_HI[s]);
    g_tails = TAILS_FREE;
    PeakFit capture = FitCaptureLine(h, "capture_" + tag);
    g_capture_meas = capture.mu;
    TF1 *f = FitInRun(h, "inrun_" + tag);
    const Double_t *p = f->GetParameters();
    Double_t x_ge = ScalePosition(p[Q_GE_E], p);
    // Position error: Ge's energy error carried onto the spectrum scale.
    Double_t x_ge_err = p[Q_B] * f->GetParError(Q_GE_E);
    Double_t x_ka1 = ScalePosition(ReferenceLines::PB_KA1, p);
    Double_t x_ka2 = ScalePosition(ReferenceLines::PB_KA2, p);
    Double_t cal_err = 0, unused = 0;
    Double_t e_ge = SourceEnergy(x_ge, am_peak.mu, am_peak.mu_err, ba_peak.mu,
                                 ba_peak.mu_err, cal_err);
    Double_t slope =
        (BA133_81 - ReferenceLines::AM241) / (ba_peak.mu - am_peak.mu);
    Double_t e_ge_err =
        TMath::Sqrt(TMath::Power(slope * x_ge_err, 2) + cal_err * cal_err);
    Double_t e_ka1 = SourceEnergy(x_ka1, am_peak.mu, am_peak.mu_err, ba_peak.mu,
                                  ba_peak.mu_err, unused);
    Double_t e_ka2 = SourceEnergy(x_ka2, am_peak.mu, am_peak.mu_err, ba_peak.mu,
                                  ba_peak.mu_err, unused);
    std::cout << "SOURCE-SCALE slice " << std::fixed << std::setprecision(1)
              << RISE_SLICE_LO[s] << "-" << RISE_SLICE_HI[s] << " ns: Am at "
              << std::setprecision(4) << am_peak.mu << " +- "
              << std::setprecision(1) << 1000 * am_peak.mu_err
              << " eV, Ba 81 at " << std::setprecision(4) << ba_peak.mu
              << " +- " << std::setprecision(1) << 1000 * ba_peak.mu_err
              << " eV (spectrum keV) | Ge: in-run " << std::setprecision(4)
              << p[Q_GE_E] << ", source scale " << e_ge << " +- "
              << std::setprecision(1) << 1000 * e_ge_err << " eV (cal "
              << 1000 * cal_err << ") | Pb Ka1 on the source scale "
              << std::showpos << 1000 * (e_ka1 - ReferenceLines::PB_KA1)
              << " eV, Ka2 " << 1000 * (e_ka2 - ReferenceLines::PB_KA2)
              << std::noshowpos << " eV" << std::endl;
    sw += 1.0 / (e_ge_err * e_ge_err);
    swx += e_ge / (e_ge_err * e_ge_err);
    Double_t w_inrun = 1.0 / TMath::Power(f->GetParError(Q_GE_E), 2);
    sw_inrun += w_inrun;
    swx_inrun += w_inrun * p[Q_GE_E];
    delete f;
    delete h;
    delete am;
    delete ba;
  }
  std::cout << "SOURCE-SCALE weighted mean: source scale " << std::fixed
            << std::setprecision(4) << swx / sw << " +- "
            << std::setprecision(1) << 1000 / TMath::Sqrt(sw)
            << " eV; in-run scale " << std::setprecision(4)
            << swx_inrun / sw_inrun << " +- " << std::setprecision(1)
            << 1000 / TMath::Sqrt(sw_inrun) << " eV" << std::endl;
}
