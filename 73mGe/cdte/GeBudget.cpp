// The CdTe Ge energy and its budget, from the in-run calibrated fit
// (CdTeAnalysis.hpp) in each rise-time slice of the Ge beam runs, the slices
// up to BUDGET_MAX_RISE_NS combined. Each systematic is the shift of the
// slices' weighted mean, with the nominal weights, when one choice changes:
//   statistics  the weighted mean's error, scaled by sqrt(chi2/ndf) of the
//               slices up to 300 ns (the slices scatter more than their
//               errors);
//   outlier     the largest shift without the highest, the lowest, or both
//               slices (CombineGeResult's outlier removal);
//   lineshape   Ge's tail fitted on its own rather than shared;
//   escapes     the Te Ka/Cd Kb tie varied by its error and taken per slice,
//               and the Cd Kb/Ka ratio from the no-Ge runs instead of
//               Geant4's (NoGeEscapes.cpp);
//   Pb K-beta   the K-beta groups' sub-line intensities: Kb3/Kb1 (and the
//               K-N2/K-N3 split taken as it) +-10%, the outer sub-lines
//               (Kb5, K-N4,5, K-O2,3) +-20%;
//   extrapolation  Ge lies below the lines that fix the scale; with a free
//               curvature of the scale Ge moves by about twice the error the
//               curvature adds, so the whole shift is charged (the larger of
//               the two);
//   slice range the slices up to 300 ns instead (printed against its
//               statistical spread, the quadrature difference of the errors);
//   slice combination  the slices scatter beyond their errors; a
//               random-effects mean instead of the scaled weighted mean.
//
// Usage: root -b -q 'GeBudget.cpp+'. Writes results/cdte_ge.result.
#include "CdTeAnalysis.hpp"

#include <TBox.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>

#include <fstream>

// The slices the result combines; above, the model fails (the no-Ge runs'
// escape fit turns unphysical and Ge rails at its limits).
const Double_t BUDGET_MAX_RISE_NS = 320.0;
const Double_t CHECK_MAX_RISE_NS = 300.0;

namespace CdTeBudget {

using namespace CdTe;

// The escape tie from NoGeEscapes.cpp: the no-Ge Pb runs, rise < 300 ns, Pb
// Ka2/Ka1 free; per slice; and the Cd Kb/Ka ratio with Pb Ka2/Ka1 fixed at
// DDEP's 0.595, against Geant4's 0.275 for captures in CdTe. Above 320 ns
// the per-slice fit is unphysical (Te Kb at 3-9% of the capture line), so
// those slices keep the rise < 300 ns value.
const Double_t TE_CDKB = 0.353, TE_CDKB_ERR = 0.043;
const Double_t TE_CDKB_SLICE[N_RISE_SLICES] = {0.084, 0.319,   0.338,  0.657,
                                               0.376, TE_CDKB, TE_CDKB};
const Double_t KB_KA_NOMINAL = 0.275, KB_KA_ALT = 0.324;

const Double_t KB3_VARIATION = 0.10;
const Double_t KB_OUTER_VARIATION = 0.20;

enum {
  V_NOMINAL,
  V_GE_OWN_TAIL,
  V_TE_UP,
  V_TE_DOWN,
  V_TE_SLICE,
  V_KBKA_ALT,
  V_KB3_UP,
  V_KB3_DOWN,
  V_KB_OUTER_UP,
  V_KB_OUTER_DOWN,
  V_SCALE2,
  N_VARIANTS
};
const char *const VARIANT_NAMES[N_VARIANTS] = {"nominal",
                                               "Ge tail on its own",
                                               "Te Ka/Cd Kb +1 sigma",
                                               "Te Ka/Cd Kb -1 sigma",
                                               "Te Ka/Cd Kb per rise slice",
                                               "Cd Kb/Ka 0.324",
                                               "Kb3/Kb1 +10%",
                                               "Kb3/Kb1 -10%",
                                               "outer Kb sub-lines +20%",
                                               "outer Kb sub-lines -20%",
                                               "scale curvature free"};

void ApplyVariant(Int_t variant, Int_t slice) {
  g_ge_own_tail = (variant == V_GE_OWN_TAIL);
  g_free_scale2 = (variant == V_SCALE2);
  Double_t te = TE_CDKB;
  if (variant == V_TE_UP)
    te = TE_CDKB + TE_CDKB_ERR;
  if (variant == V_TE_DOWN)
    te = TE_CDKB - TE_CDKB_ERR;
  if (variant == V_TE_SLICE)
    te = TE_CDKB_SLICE[slice];
  SetEscapeTie((variant == V_KBKA_ALT) ? KB_KA_ALT : KB_KA_NOMINAL, te);
  ResetKbGroups();
  Double_t kb3_factor = 1.0, outer_factor = 1.0;
  if (variant == V_KB3_UP)
    kb3_factor = 1.0 + KB3_VARIATION;
  if (variant == V_KB3_DOWN)
    kb3_factor = 1.0 - KB3_VARIATION;
  if (variant == V_KB_OUTER_UP)
    outer_factor = 1.0 + KB_OUTER_VARIATION;
  if (variant == V_KB_OUTER_DOWN)
    outer_factor = 1.0 - KB_OUTER_VARIATION;
  g_pb_kb1[1].weight *= kb3_factor;   // K-M2, Kb3
  g_pb_kb2[1].weight *= kb3_factor;   // K-N2
  g_pb_kb1[2].weight *= outer_factor; // Kb5
  g_pb_kb2[2].weight *= outer_factor; // K-N4,5
  g_pb_kb2[3].weight *= outer_factor; // K-O2,3
}

// Each slice's Ge per variant, and its fit error.
struct SliceResults {
  std::vector<Double_t> value[N_VARIANTS], error[N_VARIANTS];
};

Double_t WeightedMean(const std::vector<Double_t> &v,
                      const std::vector<Double_t> &e, Int_t skip1, Int_t skip2,
                      Double_t &err) {
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < v.size(); i++) {
    if ((Int_t)i == skip1 || (Int_t)i == skip2)
      continue;
    sw += 1.0 / (e[i] * e[i]);
    swx += v[i] / (e[i] * e[i]);
  }
  err = 1.0 / TMath::Sqrt(sw);
  return swx / sw;
}

// The DerSimonian-Laird random-effects mean: the slices' excess scatter as a
// common variance tau^2 added to each slice's own.
Double_t RandomEffectsMean(const std::vector<Double_t> &v,
                           const std::vector<Double_t> &e, Double_t &err,
                           Double_t &tau) {
  Double_t fixed_err = 0;
  Double_t fixed = WeightedMean(v, e, -1, -1, fixed_err);
  Double_t q = 0, sw = 0, sw2 = 0;
  for (size_t i = 0; i < v.size(); i++) {
    Double_t w = 1.0 / (e[i] * e[i]);
    q += w * (v[i] - fixed) * (v[i] - fixed);
    sw += w;
    sw2 += w * w;
  }
  Double_t tau2 = TMath::Max(0.0, (q - (v.size() - 1.0)) / (sw - sw2 / sw));
  tau = TMath::Sqrt(tau2);
  std::vector<Double_t> inflated;
  for (size_t i = 0; i < e.size(); i++)
    inflated.push_back(TMath::Sqrt(e[i] * e[i] + tau2));
  return WeightedMean(v, inflated, -1, -1, err);
}

// The slices up to max_hi_ns.
void Select(const SliceResults &all, Double_t max_hi_ns, SliceResults &used,
            std::vector<Int_t> &slice_index) {
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    if (RISE_SLICE_HI[s] > max_hi_ns)
      continue;
    slice_index.push_back(s);
    for (Int_t k = 0; k < N_VARIANTS; k++) {
      used.value[k].push_back(all.value[k][s]);
      used.error[k].push_back(all.error[k][s]);
    }
  }
}

struct Budget {
  Double_t mean = 0, weighted_err = 0, scale = 1, stat = 0;
  Double_t outlier = 0, lineshape = 0, escapes = 0, kb = 0;
  Double_t extrapolation = 0, total = 0;
  Double_t shift[N_VARIANTS];
  Double_t scale2_err = 0;
};

Budget Summarize(const TString &label, const SliceResults &all,
                 Double_t max_hi_ns) {
  Budget b;
  SliceResults used;
  std::vector<Int_t> slice_index;
  Select(all, max_hi_ns, used, slice_index);
  const std::vector<Double_t> &v = used.value[V_NOMINAL];
  const std::vector<Double_t> &e = used.error[V_NOMINAL];
  const Int_t n = (Int_t)v.size();
  b.mean = WeightedMean(v, e, -1, -1, b.weighted_err);

  std::vector<Double_t> good_v, good_e;
  for (Int_t i = 0; i < n; i++)
    if (RISE_SLICE_HI[slice_index[i]] <= 300.0) {
      good_v.push_back(v[i]);
      good_e.push_back(e[i]);
    }
  Double_t good_err = 0;
  Double_t good_mean = WeightedMean(good_v, good_e, -1, -1, good_err);
  Double_t chi2_good = 0;
  for (size_t i = 0; i < good_v.size(); i++)
    chi2_good += TMath::Power((good_v[i] - good_mean) / good_e[i], 2);
  b.scale = TMath::Sqrt(TMath::Max(chi2_good / (good_v.size() - 1), 1.0));
  b.stat = b.weighted_err * b.scale;

  Int_t i_hi = 0, i_lo = 0;
  for (Int_t i = 1; i < n; i++) {
    if (v[i] > v[i_hi])
      i_hi = i;
    if (v[i] < v[i_lo])
      i_lo = i;
  }
  Double_t unused = 0;
  b.outlier =
      TMath::Max(TMath::Abs(WeightedMean(v, e, i_hi, -1, unused) - b.mean),
                 TMath::Abs(WeightedMean(v, e, i_lo, -1, unused) - b.mean));
  b.outlier = TMath::Max(
      b.outlier, TMath::Abs(WeightedMean(v, e, i_hi, i_lo, unused) - b.mean));

  for (Int_t k = 1; k < N_VARIANTS; k++)
    b.shift[k] = WeightedMean(used.value[k], e, -1, -1, unused) - b.mean;
  b.lineshape = TMath::Abs(b.shift[V_GE_OWN_TAIL]);
  Double_t escape_ratio =
      TMath::Max(TMath::Abs(b.shift[V_TE_UP]), TMath::Abs(b.shift[V_TE_DOWN]));
  Double_t escape_slices = TMath::Abs(b.shift[V_TE_SLICE]);
  Double_t escape_kbka = TMath::Abs(b.shift[V_KBKA_ALT]);
  b.escapes =
      TMath::Sqrt(escape_ratio * escape_ratio + escape_slices * escape_slices +
                  escape_kbka * escape_kbka);
  Double_t kb3 = TMath::Max(TMath::Abs(b.shift[V_KB3_UP]),
                            TMath::Abs(b.shift[V_KB3_DOWN]));
  Double_t kb_outer = TMath::Max(TMath::Abs(b.shift[V_KB_OUTER_UP]),
                                 TMath::Abs(b.shift[V_KB_OUTER_DOWN]));
  b.kb = TMath::Sqrt(kb3 * kb3 + kb_outer * kb_outer);
  WeightedMean(used.value[V_SCALE2], used.error[V_SCALE2], -1, -1,
               b.scale2_err);
  Double_t curvature_cost = TMath::Sqrt(TMath::Max(
      b.scale2_err * b.scale2_err - b.weighted_err * b.weighted_err, 0.0));
  b.extrapolation = TMath::Max(TMath::Abs(b.shift[V_SCALE2]), curvature_cost);
  b.total = TMath::Sqrt(b.stat * b.stat + b.outlier * b.outlier +
                        b.lineshape * b.lineshape + b.escapes * b.escapes +
                        b.kb * b.kb + b.extrapolation * b.extrapolation);

  TString tag = "BUDGET [" + label + "]";
  std::cout << tag << " weighted mean of " << n << " slices: " << std::fixed
            << std::setprecision(4) << b.mean << " keV" << std::endl;
  std::cout << tag << "   statistics " << std::setprecision(1)
            << 1000 * b.weighted_err << " eV x sqrt(chi2/ndf) "
            << std::setprecision(2) << b.scale
            << " (rise <= 300 ns: " << chi2_good << "/" << good_v.size() - 1
            << ") = " << std::setprecision(1) << 1000 * b.stat << " eV"
            << std::endl;
  std::cout << tag << "   outlier removal " << 1000 * b.outlier
            << " eV (highest " << RISE_SLICE_LO[slice_index[i_hi]] << "-"
            << RISE_SLICE_HI[slice_index[i_hi]] << " ns, lowest "
            << RISE_SLICE_LO[slice_index[i_lo]] << "-"
            << RISE_SLICE_HI[slice_index[i_lo]] << " ns)" << std::endl;
  for (Int_t k = 1; k < N_VARIANTS; k++)
    std::cout << tag << "   " << VARIANT_NAMES[k] << ": shift " << std::showpos
              << 1000 * b.shift[k] << std::noshowpos << " eV" << std::endl;
  std::cout << tag << "   escapes: ratio " << 1000 * escape_ratio
            << " (+) per slice " << 1000 * escape_slices << " (+) Kb/Ka "
            << 1000 * escape_kbka << " = " << 1000 * b.escapes << " eV"
            << std::endl;
  std::cout << tag << "   Pb K-beta: Kb3/Kb1 " << 1000 * kb3
            << " (+) outer sub-lines " << 1000 * kb_outer << " = "
            << 1000 * b.kb << " eV" << std::endl;
  std::cout << tag << "   extrapolation: shift " << std::showpos
            << 1000 * b.shift[V_SCALE2] << std::noshowpos
            << " eV with the curvature free; the weighted mean's error "
            << 1000 * b.weighted_err << " -> " << 1000 * b.scale2_err << " eV, "
            << 1000 * curvature_cost
            << " eV in quadrature: " << 1000 * b.extrapolation << " eV"
            << std::endl;
  std::cout << tag << " RESULT: " << std::setprecision(4) << b.mean
            << " keV +/- " << std::setprecision(1) << 1000 * b.stat
            << " (stat) +/- " << 1000 * b.outlier << " (outlier) +/- "
            << 1000 * b.lineshape << " (lineshape) +/- " << 1000 * b.escapes
            << " (escapes) +/- " << 1000 * b.kb << " (Pb K-beta) +/- "
            << 1000 * b.extrapolation << " (extrapolation) = +/- "
            << 1000 * b.total << " eV   (" << std::showpos
            << 1000 * (b.mean - ReferenceLines::GE73M_ENSDF) << std::noshowpos
            << " eV vs ENSDF)" << std::endl;
  return b;
}

// The result: the budget of the slices up to BUDGET_MAX_RISE_NS, plus how the
// slices are chosen and combined.
struct FinalResult {
  Budget budget;
  Double_t slice_range = 0, slice_combination = 0, total = 0;
};

FinalResult Finalize(const SliceResults &all, const Budget &b,
                     const Budget &narrow) {
  FinalResult result;
  result.budget = b;
  Double_t range_shift = narrow.mean - b.mean;
  Double_t range_spread =
      TMath::Sqrt(TMath::Max(narrow.weighted_err * narrow.weighted_err -
                                 b.weighted_err * b.weighted_err,
                             0.0));
  result.slice_range = TMath::Abs(range_shift);
  std::cout << "SLICES up to " << std::fixed << std::setprecision(0)
            << CHECK_MAX_RISE_NS << " instead of " << BUDGET_MAX_RISE_NS
            << " ns: shift " << std::showpos << std::setprecision(1)
            << 1000 * range_shift << std::noshowpos
            << " eV, statistical spread " << 1000 * range_spread << " eV"
            << std::endl;

  SliceResults used;
  std::vector<Int_t> slice_index;
  Select(all, BUDGET_MAX_RISE_NS, used, slice_index);
  Double_t re_err = 0, tau = 0;
  Double_t re_mean = RandomEffectsMean(used.value[V_NOMINAL],
                                       used.error[V_NOMINAL], re_err, tau);
  result.slice_combination = TMath::Abs(re_mean - b.mean);
  std::cout << "SLICES random-effects mean: " << std::fixed
            << std::setprecision(4) << re_mean << " +/- "
            << std::setprecision(1) << 1000 * re_err
            << " eV (between-slice spread " << 1000 * tau
            << " eV), against the scaled weighted mean " << std::setprecision(4)
            << b.mean << " +/- " << std::setprecision(1) << 1000 * b.stat
            << " eV: shift " << std::showpos << 1000 * (re_mean - b.mean)
            << std::noshowpos << " eV" << std::endl;

  result.total =
      TMath::Sqrt(b.total * b.total + result.slice_range * result.slice_range +
                  result.slice_combination * result.slice_combination);
  std::cout << "FINAL CdTe: " << std::setprecision(4) << b.mean << " keV +/- "
            << std::setprecision(1) << 1000 * b.stat << " (stat) +/- "
            << 1000 * b.outlier << " (outlier) +/- " << 1000 * b.lineshape
            << " (lineshape) +/- " << 1000 * b.escapes << " (escapes) +/- "
            << 1000 * b.kb << " (Pb K-beta) +/- " << 1000 * b.extrapolation
            << " (extrapolation) +/- " << 1000 * result.slice_range
            << " (slice range) +/- " << 1000 * result.slice_combination
            << " (slice combination) = +/- " << 1000 * result.total << " eV   ("
            << std::showpos << 1000 * (b.mean - ReferenceLines::GE73M_ENSDF)
            << std::noshowpos << " eV vs ENSDF)" << std::endl;
  return result;
}

// Ge per rise slice (nominal fit, its fit error), the slices used filled and
// the rest open, with the result as a band and ENSDF for reference.
void PlotGeVsRise(const SliceResults &all, Double_t mean, Double_t total) {
  TGraphErrors *used = new TGraphErrors();
  TGraphErrors *excluded = new TGraphErrors();
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    Double_t slice_lo = (RISE_SLICE_LO[s] > 0) ? RISE_SLICE_LO[s] : 250.0;
    TGraphErrors *g = (RISE_SLICE_HI[s] > BUDGET_MAX_RISE_NS) ? excluded : used;
    Int_t i = g->GetN();
    g->SetPoint(i, 0.5 * (slice_lo + RISE_SLICE_HI[s]),
                all.value[V_NOMINAL][s]);
    g->SetPointError(i, 0.5 * (RISE_SLICE_HI[s] - slice_lo),
                     all.error[V_NOMINAL][s]);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(used, kBlue + 1,
                                ";Rise time 10-90% [ns];E_{#gamma} [keV]");
  PlottingUtils::ConfigureGraph(excluded, kGray + 2, "");
  used->SetMarkerStyle(20);
  excluded->SetMarkerStyle(24);
  const Double_t x_lo = 240.0, x_hi = 460.0;
  used->GetXaxis()->SetLimits(x_lo, x_hi);
  used->GetYaxis()->SetRangeUser(67.9, 69.5);
  used->Draw("AP");
  TBox *band = new TBox(x_lo, mean - total, x_hi, mean + total);
  band->SetFillColorAlpha(kBlue + 1, 0.2);
  band->Draw();
  TLine *line = new TLine(x_lo, mean, x_hi, mean);
  line->SetLineColor(kBlue + 1);
  line->SetLineStyle(2);
  line->Draw();
  TLine *ensdf = new TLine(x_lo, ReferenceLines::GE73M_ENSDF, x_hi,
                           ReferenceLines::GE73M_ENSDF);
  ensdf->SetLineColor(kRed + 1);
  ensdf->SetLineStyle(3);
  ensdf->Draw();
  used->Draw("P SAME");
  excluded->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.88, 0.66, 0.9);
  leg->AddEntry(used, "Slice, in-run calibration (fit error)", "pe");
  leg->AddEntry(
      excluded,
      Form("Above %.0f ns: model fails, not used", BUDGET_MAX_RISE_NS), "pe");
  leg->AddEntry(band,
                Form("Budget up to %.0f ns: %.4f #pm %.4f keV",
                     BUDGET_MAX_RISE_NS, mean, total),
                "f");
  leg->AddEntry(ensdf, "ENSDF 68.752 keV", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_cdte_ge_vs_rise_slice", "",
                            PlotSaveOptions::kLINEAR);
}

void WriteResult(const TString &path, const FinalResult &result) {
  const Budget &b = result.budget;
  std::ofstream out(path.Data());
  out << std::fixed << std::setprecision(6);
  out << "# CdTe (PSI) Ge energy, GeBudget.cpp; keV" << std::endl;
  out << "ge " << b.mean << std::endl;
  out << "total " << result.total << std::endl;
  out << "stat " << b.stat << std::endl;
  out << "outlier " << b.outlier << std::endl;
  out << "lineshape " << b.lineshape << std::endl;
  out << "escapes " << b.escapes << std::endl;
  out << "pb_kbeta " << b.kb << std::endl;
  out << "extrapolation " << b.extrapolation << std::endl;
  out << "slice_range " << result.slice_range << std::endl;
  out << "slice_combination " << result.slice_combination << std::endl;
  std::cout << "Wrote " << path << std::endl;
}

} // namespace CdTeBudget

void GeBudget() {
  using namespace CdTeBudget;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root +
                                    "/plots/studies/psi_cdte/ge_rise_slices",
                                project_root + "/root_files");
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;
  SliceResults all;
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    TH1D *h = FillRise(GeRuns(), Form("budget_h%d", s), RISE_SLICE_LO[s],
                       RISE_SLICE_HI[s]);
    g_capture_meas = FitCaptureLine(h, Form("budget%d", s)).mu;
    for (Int_t v = 0; v < N_VARIANTS; v++) {
      ApplyVariant(v, s);
      TF1 *f = FitInRun(h, Form("budget%d_v%d", s, v));
      all.value[v].push_back(f->GetParameter(Q_GE_E));
      all.error[v].push_back(f->GetParError(Q_GE_E));
      std::cout << "BUDGET slice " << std::fixed << std::setprecision(0)
                << RISE_SLICE_LO[s] << "-" << RISE_SLICE_HI[s] << " ns, "
                << VARIANT_NAMES[v] << ": Ge " << std::fixed
                << std::setprecision(4) << f->GetParameter(Q_GE_E) << " +/- "
                << std::setprecision(1) << 1000 * f->GetParError(Q_GE_E)
                << " eV, chi2/ndf " << std::setprecision(2)
                << f->GetChisquare() / f->GetNDF() << ", Cd Kb escape / Ge "
                << std::setprecision(3)
                << f->GetParameter(Q_ESC_A) / f->GetParameter(Q_GE_A)
                << std::endl;
      if (v == V_NOMINAL)
        DrawInRun(h, f,
                  Form("psi_cdte_ge_slice_%.0f_%.0fns", RISE_SLICE_LO[s],
                       RISE_SLICE_HI[s]),
                  Form("Rise time %.0f-%.0f ns, Ge runs", RISE_SLICE_LO[s],
                       RISE_SLICE_HI[s]));
      delete f;
    }
    ApplyVariant(V_NOMINAL, s);
    delete h;
  }

  Summarize("all slices", all, 450.0);
  Budget narrow = Summarize(Form("rise <= %.0f ns", CHECK_MAX_RISE_NS), all,
                            CHECK_MAX_RISE_NS);
  Budget b = Summarize(Form("rise <= %.0f ns", BUDGET_MAX_RISE_NS), all,
                       BUDGET_MAX_RISE_NS);
  FinalResult result = Finalize(all, b, narrow);
  PlotGeVsRise(all, b.mean, result.total);
  WriteResult(project_root + "/results/cdte_ge.result", result);
}
