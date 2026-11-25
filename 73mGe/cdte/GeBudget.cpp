// The CdTe Ge energy and its budget. The Ge beam runs are split into cells,
// one day in one rise-time slice (CdTeAnalysis.hpp). In each cell the in-run
// fit gives the positions of Ge and the Pb K-alpha lines, and Ge is read on
// that day's source-run scale (Am-241 and Ba-133), carried to the beam runs
// by the Pb K-alpha gain (GainTransfer). The cells up to BUDGET_MAX_RISE_NS
// are combined in a weighted mean, their statistical errors as weights. Each
// systematic is the shift of that mean, with the nominal weights, when one
// choice changes:
//   statistics  the weighted mean's error, scaled by sqrt(chi2/ndf) of the
//               cells (when they scatter more than their errors);
//   calibration the source runs' Am and Ba positions by their fit errors
//               (each cell has its own source runs: its day, its slice);
//   outlier     the largest shift without the highest, the lowest, or both
//               cells (CombineGeResult's outlier removal);
//   lineshape   Ge's width fitted on its own rather than shared (Ge's rise
//               distribution, so its mix of widths, is its own). With its
//               tails free as well Ge takes over the escape lines beside it,
//               so that is printed as a check;
//   escapes     the Te Ka/Cd Kb tie varied by its error and taken per slice,
//               and the Cd Kb/Ka ratio from the no-Ge runs instead of
//               Geant4's (NoGeEscapes.cpp);
//   background  the continuum's slope, held at the side bands' value
//               (CdTeAnalysis.hpp): the largest shift with that slope +-1
//               sigma or a flat continuum;
//   Pb K-beta   the K-beta groups' sub-line intensities: Kb3/Kb1 (and the
//               K-N2/K-N3 split taken as it) +-10%, the outer sub-lines
//               (Kb5, K-N4,5, K-O2,3) +-20%;
//   curvature   the scale's relative curvature kappa, in the source scale and
//               the in-run fit alike, at the source runs' value +- its error
//               (SourceCurvature.cpp). The beam fit with kappa free is
//               printed as a check of that value, not charged;
//   slice range without the slowest slice's cells (printed against its
//               statistical spread, the quadrature difference of the errors);
//   slice combination  a random-effects mean of the cells instead of the
//               scaled weighted mean.
// Checks, not charged: the step (held at the source runs' value) free, and
// the tails (free) taken from the capture line, or the low tail from the
// source runs. The former in-run scale (Ge's energy straight from the in-run
// fit) is printed for comparison.
//
// Usage: root -b -q 'GeBudget.cpp+'. Writes results/cdte_ge.result.
#include "CdTeAnalysis.hpp"

#include <TBox.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>

#include <fstream>

// The result combines every slice (below 281 ns, CdTeAnalysis.hpp); the
// slice-range check leaves out the slowest.
const Double_t BUDGET_MAX_RISE_NS = CdTe::RESULT_MAX_RISE_NS;
const Double_t CHECK_MAX_RISE_NS = CdTe::RISE_SLICE_HI[CdTe::N_RISE_SLICES - 2];

namespace CdTeBudget {

using namespace CdTe;

// The escape tie from NoGeEscapes.cpp, the no-Ge Pb runs with Pb Ka2/Ka1 at
// DDEP's 0.595 (free, it goes to 0.30-0.54, trading with the Cd Ka escapes
// under Pb Ka2): rise < 281 ns; per slice; and the Cd Kb/Ka ratio, against
// Geant4's 0.275 for captures in CdTe; the continuum's slope held at the side
// bands', the step at the source runs' and the tails free.
// A slice whose own fit is unphysical keeps the rise < 281 ns value.
const Double_t TE_CDKB = 0.352, TE_CDKB_ERR = 0.040;
const Double_t TE_CDKB_SLICE[N_RISE_SLICES] = {0.245, 0.353, 0.395, 0.341};
const Double_t KB_KA_NOMINAL = 0.275, KB_KA_ALT = 0.298;

const Double_t KB3_VARIATION = 0.10;
const Double_t KB_OUTER_VARIATION = 0.20;

// A cell: day index d (GE_DAYS) and rise slice s.
const Int_t N_CELLS = N_GE_DAYS * N_RISE_SLICES;
Int_t CellDayIndex(Int_t cell) { return cell / N_RISE_SLICES; }
Int_t CellSlice(Int_t cell) { return cell % N_RISE_SLICES; }
Int_t CellDay(Int_t cell) { return GE_DAYS[CellDayIndex(cell)]; }

enum {
  V_NOMINAL,
  V_GE_OWN_WIDTH,
  V_GE_OWN_SHAPE,
  V_TE_UP,
  V_TE_DOWN,
  V_TE_SLICE,
  V_KBKA_ALT,
  V_KB3_UP,
  V_KB3_DOWN,
  V_KB_OUTER_UP,
  V_KB_OUTER_DOWN,
  V_KAPPA_UP,
  V_KAPPA_DOWN,
  V_KAPPA_FREE,
  V_SLOPE_UP,
  V_SLOPE_DOWN,
  V_FLAT,
  V_STEP_FREE,
  V_TAILS_CAPTURE,
  V_TAILS_SOURCE,
  N_VARIANTS
};
const char *const VARIANT_NAMES[N_VARIANTS] = {
    "nominal",
    "Ge width on its own",
    "Ge width and tails on their own (check)",
    "Te Ka/Cd Kb +1 sigma",
    "Te Ka/Cd Kb -1 sigma",
    "Te Ka/Cd Kb per rise slice",
    "Cd Kb/Ka from the no-Ge runs",
    "Kb3/Kb1 +10%",
    "Kb3/Kb1 -10%",
    "outer Kb sub-lines +20%",
    "outer Kb sub-lines -20%",
    "kappa +1 sigma",
    "kappa -1 sigma",
    "kappa free (check)",
    "continuum slope +1 sigma",
    "continuum slope -1 sigma",
    "flat continuum",
    "step free (check)",
    "capture line's tails (check)",
    "source low tail (check)"};

void ApplyVariant(Int_t variant, Int_t slice) {
  g_tails = TAILS_FREE;
  if (variant == V_TAILS_CAPTURE)
    g_tails = TAILS_CAPTURE;
  if (variant == V_TAILS_SOURCE)
    g_tails = TAILS_SOURCE;
  g_ge_own_width = (variant == V_GE_OWN_WIDTH);
  g_ge_own_shape = (variant == V_GE_OWN_SHAPE);
  g_free_kappa = (variant == V_KAPPA_FREE);
  g_free_step = (variant == V_STEP_FREE);
  g_flat_continuum = (variant == V_FLAT);
  g_slope_shift = 0;
  if (variant == V_SLOPE_UP)
    g_slope_shift = 1.0;
  if (variant == V_SLOPE_DOWN)
    g_slope_shift = -1.0;
  g_kappa = KAPPA_SOURCE;
  if (variant == V_KAPPA_UP)
    g_kappa = KAPPA_SOURCE + KAPPA_SOURCE_ERR;
  if (variant == V_KAPPA_DOWN)
    g_kappa = KAPPA_SOURCE - KAPPA_SOURCE_ERR;
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

// Per cell: Ge per variant and its statistical error; for the nominal, Ge's
// derivatives with respect to its source line's Am and Ba positions and
// their errors, and Ge on the in-run scale; and the beam's own kappa from the
// kappa-free check.
struct CellResults {
  std::vector<Double_t> value[N_VARIANTS], error[N_VARIANTS];
  std::vector<Double_t> d_am, d_ba, x_am_err, x_ba_err;
  std::vector<Int_t> cell;
  std::vector<Double_t> in_run, in_run_err;
  std::vector<Double_t> kappa, kappa_err;
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

Double_t Chi2About(const std::vector<Double_t> &v,
                   const std::vector<Double_t> &e, Double_t mean) {
  Double_t chi2 = 0;
  for (size_t i = 0; i < v.size(); i++)
    chi2 += TMath::Power((v[i] - mean) / e[i], 2);
  return chi2;
}

// The DerSimonian-Laird random-effects mean: the cells' excess scatter as a
// common variance tau^2 added to each cell's own.
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

// The cells of the slices up to max_hi_ns.
CellResults Select(const CellResults &all, Double_t max_hi_ns) {
  CellResults used;
  for (size_t i = 0; i < all.cell.size(); i++) {
    if (RISE_SLICE_HI[CellSlice(all.cell[i])] > max_hi_ns)
      continue;
    for (Int_t k = 0; k < N_VARIANTS; k++) {
      used.value[k].push_back(all.value[k][i]);
      used.error[k].push_back(all.error[k][i]);
    }
    used.d_am.push_back(all.d_am[i]);
    used.d_ba.push_back(all.d_ba[i]);
    used.x_am_err.push_back(all.x_am_err[i]);
    used.x_ba_err.push_back(all.x_ba_err[i]);
    used.cell.push_back(all.cell[i]);
    used.in_run.push_back(all.in_run[i]);
    used.in_run_err.push_back(all.in_run_err[i]);
    used.kappa.push_back(all.kappa[i]);
    used.kappa_err.push_back(all.kappa_err[i]);
  }
  return used;
}

// The weighted mean's error from the source lines' Am and Ba positions.
Double_t CalibrationError(const CellResults &used) {
  const std::vector<Double_t> &e = used.error[V_NOMINAL];
  Double_t sw = 0, variance = 0;
  for (size_t i = 0; i < e.size(); i++) {
    Double_t w = 1.0 / (e[i] * e[i]);
    sw += w;
    variance += w * w *
                (TMath::Power(used.d_am[i] * used.x_am_err[i], 2) +
                 TMath::Power(used.d_ba[i] * used.x_ba_err[i], 2));
  }
  return TMath::Sqrt(variance) / sw;
}

TString CellName(Int_t cell) {
  Int_t s = CellSlice(cell);
  return Form("10/%02d %.1f-%.1f ns", CellDay(cell), RISE_SLICE_LO[s],
              RISE_SLICE_HI[s]);
}

struct Budget {
  Double_t mean = 0, weighted_err = 0, chi2 = 0, scale = 1, stat = 0;
  Double_t calibration = 0, outlier = 0, lineshape = 0, escapes = 0;
  Double_t background = 0, kb = 0, curvature = 0, total = 0;
  Double_t shift[N_VARIANTS];
};

Budget Summarize(const TString &label, const CellResults &all,
                 Double_t max_hi_ns) {
  Budget b;
  CellResults used = Select(all, max_hi_ns);
  const std::vector<Double_t> &v = used.value[V_NOMINAL];
  const std::vector<Double_t> &e = used.error[V_NOMINAL];
  const Int_t n = (Int_t)v.size();
  b.mean = WeightedMean(v, e, -1, -1, b.weighted_err);
  b.chi2 = Chi2About(v, e, b.mean);
  b.scale = TMath::Sqrt(TMath::Max(b.chi2 / (n - 1), 1.0));
  b.stat = b.weighted_err * b.scale;
  b.calibration = CalibrationError(used);

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
  b.lineshape = TMath::Abs(b.shift[V_GE_OWN_WIDTH]);
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
  b.background = TMath::Max(TMath::Abs(b.shift[V_SLOPE_UP]),
                            TMath::Abs(b.shift[V_SLOPE_DOWN]));
  b.background = TMath::Max(b.background, TMath::Abs(b.shift[V_FLAT]));
  b.curvature = TMath::Max(TMath::Abs(b.shift[V_KAPPA_UP]),
                           TMath::Abs(b.shift[V_KAPPA_DOWN]));
  Double_t beam_kappa_err = 0;
  Double_t beam_kappa =
      WeightedMean(used.kappa, used.kappa_err, -1, -1, beam_kappa_err);
  Double_t in_run_err = 0;
  Double_t in_run =
      WeightedMean(used.in_run, used.in_run_err, -1, -1, in_run_err);
  Double_t in_run_chi2 = Chi2About(used.in_run, used.in_run_err, in_run);
  b.total = TMath::Sqrt(b.stat * b.stat + b.calibration * b.calibration +
                        b.outlier * b.outlier + b.lineshape * b.lineshape +
                        b.escapes * b.escapes + b.background * b.background +
                        b.kb * b.kb + b.curvature * b.curvature);

  TString tag = "BUDGET [" + label + "]";
  std::cout << tag << " weighted mean of " << n << " cells: " << std::fixed
            << std::setprecision(4) << b.mean << " keV" << std::endl;
  std::cout << tag << "   statistics " << std::setprecision(1)
            << 1000 * b.weighted_err << " eV x sqrt(chi2/ndf) "
            << std::setprecision(2) << b.scale << " (" << b.chi2 << "/" << n - 1
            << ") = " << std::setprecision(1) << 1000 * b.stat << " eV"
            << std::endl;
  std::cout << tag << "   calibration (source runs' Am and Ba positions) "
            << 1000 * b.calibration << " eV" << std::endl;
  std::cout << tag << "   outlier removal " << 1000 * b.outlier
            << " eV (highest " << CellName(used.cell[i_hi]) << ", lowest "
            << CellName(used.cell[i_lo]) << ")" << std::endl;
  for (Int_t k = 1; k < N_VARIANTS; k++)
    std::cout << tag << "   " << VARIANT_NAMES[k] << ": shift " << std::showpos
              << 1000 * b.shift[k] << std::noshowpos << " eV" << std::endl;
  std::cout << tag << "   escapes: ratio " << 1000 * escape_ratio
            << " (+) per slice " << 1000 * escape_slices << " (+) Kb/Ka "
            << 1000 * escape_kbka << " = " << 1000 * b.escapes << " eV"
            << std::endl;
  std::cout << tag << "   background: the largest of the slope +-1 sigma and "
            << "flat shifts, " << 1000 * b.background << " eV. Checks: step "
            << "free, Ge " << std::showpos << 1000 * b.shift[V_STEP_FREE]
            << " eV; capture line's tails, Ge "
            << 1000 * b.shift[V_TAILS_CAPTURE] << " eV; source low tail, Ge "
            << 1000 * b.shift[V_TAILS_SOURCE] << std::noshowpos << " eV"
            << std::endl;
  std::cout << tag << "   Pb K-beta: Kb3/Kb1 " << 1000 * kb3
            << " (+) outer sub-lines " << 1000 * kb_outer << " = "
            << 1000 * b.kb << " eV" << std::endl;
  std::cout << tag << "   curvature: kappa " << std::scientific
            << std::setprecision(2) << KAPPA_SOURCE << " +- "
            << KAPPA_SOURCE_ERR << std::fixed
            << " /keV from the source runs: " << std::setprecision(1)
            << 1000 * b.curvature
            << " eV. Check, kappa free in the beam fit: " << std::scientific
            << std::setprecision(2) << beam_kappa << " +- " << beam_kappa_err
            << std::fixed << " /keV, Ge " << std::showpos
            << std::setprecision(1) << 1000 * b.shift[V_KAPPA_FREE]
            << std::noshowpos << " eV" << std::endl;
  std::cout << tag
            << "   for comparison, the in-run scale: " << std::setprecision(4)
            << in_run << " +/- " << std::setprecision(1) << 1000 * in_run_err
            << " eV (chi2 " << in_run_chi2 << "/" << n - 1 << ")" << std::endl;
  std::cout << tag << " RESULT: " << std::setprecision(4) << b.mean
            << " keV +/- " << std::setprecision(1) << 1000 * b.stat
            << " (stat) +/- " << 1000 * b.calibration << " (calibration) +/- "
            << 1000 * b.outlier << " (outlier) +/- " << 1000 * b.lineshape
            << " (lineshape) +/- " << 1000 * b.escapes << " (escapes) +/- "
            << 1000 * b.background << " (background) +/- " << 1000 * b.kb
            << " (Pb K-beta) +/- " << 1000 * b.curvature
            << " (curvature) = +/- " << 1000 * b.total << " eV   ("
            << std::showpos << 1000 * (b.mean - ReferenceLines::GE73M_ENSDF)
            << std::noshowpos << " eV vs ENSDF)" << std::endl;
  return b;
}

// The nominal cells' weighted mean per day (over its slices) and per slice
// (over the days), with chi2 within each and between them.
void PrintGroups(const CellResults &used) {
  const std::vector<Double_t> &v = used.value[V_NOMINAL];
  const std::vector<Double_t> &e = used.error[V_NOMINAL];
  for (Int_t by_day = 1; by_day >= 0; by_day--) {
    const Int_t n_groups = by_day ? N_GE_DAYS : N_RISE_SLICES;
    std::vector<Double_t> means, errs;
    for (Int_t g = 0; g < n_groups; g++) {
      std::vector<Double_t> gv, ge;
      for (size_t i = 0; i < v.size(); i++) {
        Int_t key =
            by_day ? CellDayIndex(used.cell[i]) : CellSlice(used.cell[i]);
        if (key != g)
          continue;
        gv.push_back(v[i]);
        ge.push_back(e[i]);
      }
      if (gv.empty())
        continue;
      Double_t err = 0;
      Double_t mean = WeightedMean(gv, ge, -1, -1, err);
      means.push_back(mean);
      errs.push_back(err);
      std::cout << "GROUP " << std::fixed << std::setprecision(1);
      if (by_day)
        std::cout << "10/" << GE_DAYS[g];
      else
        std::cout << "slice " << RISE_SLICE_LO[g] << "-" << RISE_SLICE_HI[g]
                  << " ns";
      std::cout << ": " << std::setprecision(4) << mean << " +/- "
                << std::setprecision(1) << 1000 * err << " eV (chi2 "
                << Chi2About(gv, ge, mean) << "/" << gv.size() - 1 << ")"
                << std::endl;
    }
    Double_t err = 0;
    Double_t mean = WeightedMean(means, errs, -1, -1, err);
    std::cout << "GROUP between " << (by_day ? "days" : "slices") << ": chi2 "
              << std::setprecision(1) << Chi2About(means, errs, mean) << "/"
              << means.size() - 1 << std::endl;
  }
}

// The result: the budget of the cells up to BUDGET_MAX_RISE_NS, plus how the
// slices are chosen and the cells combined.
struct FinalResult {
  Budget budget;
  Double_t slice_range = 0, slice_combination = 0, total = 0;
};

FinalResult Finalize(const CellResults &all, const Budget &b,
                     const Budget &narrow) {
  FinalResult result;
  result.budget = b;
  Double_t range_shift = narrow.mean - b.mean;
  Double_t range_spread =
      TMath::Sqrt(TMath::Max(narrow.weighted_err * narrow.weighted_err -
                                 b.weighted_err * b.weighted_err,
                             0.0));
  result.slice_range = TMath::Abs(range_shift);
  std::cout << "SLICES up to " << std::fixed << std::setprecision(1)
            << CHECK_MAX_RISE_NS << " instead of " << BUDGET_MAX_RISE_NS
            << " ns: shift " << std::showpos << std::setprecision(1)
            << 1000 * range_shift << std::noshowpos
            << " eV, statistical spread " << 1000 * range_spread << " eV"
            << std::endl;

  CellResults used = Select(all, BUDGET_MAX_RISE_NS);
  PrintGroups(used);
  Double_t re_err = 0, tau = 0;
  Double_t re_mean = RandomEffectsMean(used.value[V_NOMINAL],
                                       used.error[V_NOMINAL], re_err, tau);
  result.slice_combination = TMath::Abs(re_mean - b.mean);
  std::cout << "CELLS random-effects mean: " << std::fixed
            << std::setprecision(4) << re_mean << " +/- "
            << std::setprecision(1) << 1000 * re_err
            << " eV (between-cell spread " << 1000 * tau
            << " eV), against the scaled weighted mean " << std::setprecision(4)
            << b.mean << " +/- " << std::setprecision(1) << 1000 * b.stat
            << " eV: shift " << std::showpos << 1000 * (re_mean - b.mean)
            << std::noshowpos << " eV" << std::endl;

  result.total =
      TMath::Sqrt(b.total * b.total + result.slice_range * result.slice_range +
                  result.slice_combination * result.slice_combination);
  std::cout << "FINAL CdTe: " << std::setprecision(4) << b.mean << " keV +/- "
            << std::setprecision(1) << 1000 * b.stat << " (stat) +/- "
            << 1000 * b.calibration << " (calibration) +/- " << 1000 * b.outlier
            << " (outlier) +/- " << 1000 * b.lineshape << " (lineshape) +/- "
            << 1000 * b.escapes << " (escapes) +/- " << 1000 * b.background
            << " (background) +/- " << 1000 * b.kb << " (Pb K-beta) +/- "
            << 1000 * b.curvature << " (curvature) +/- "
            << 1000 * result.slice_range << " (slice range) +/- "
            << 1000 * result.slice_combination << " (slice combination) = +/- "
            << 1000 * result.total << " eV   (" << std::showpos
            << 1000 * (b.mean - ReferenceLines::GE73M_ENSDF) << std::noshowpos
            << " eV vs ENSDF)" << std::endl;
  return result;
}

// Ge per cell (nominal, its statistical error), each day in its own colour,
// with the result as a band and ENSDF for reference.
void PlotGeVsRise(const CellResults &all, Double_t mean, Double_t total) {
  const Int_t colors[N_GE_DAYS] = {kBlue + 1, kRed + 1};
  const Double_t x_lo = 245.0, x_hi = 285.0;
  TGraphErrors *days[N_GE_DAYS];
  Double_t y_lo = mean - total, y_hi = mean + total;
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    days[d] = new TGraphErrors();
  for (size_t i = 0; i < all.cell.size(); i++) {
    Int_t d = CellDayIndex(all.cell[i]), s = CellSlice(all.cell[i]);
    Double_t slice_lo = (RISE_SLICE_LO[s] > 0) ? RISE_SLICE_LO[s] : 250.0;
    Double_t x = 0.5 * (slice_lo + RISE_SLICE_HI[s]) + 0.4 * (d - 0.5);
    Double_t y = all.value[V_NOMINAL][i], err = all.error[V_NOMINAL][i];
    Int_t point = days[d]->GetN();
    days[d]->SetPoint(point, x, y);
    days[d]->SetPointError(point, 0, err);
    y_lo = TMath::Min(y_lo, y - err);
    y_hi = TMath::Max(y_hi, y + err);
  }
  std::vector<TString> labels;
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    labels.push_back(Form("10/%02d", GE_DAYS[d]));
  labels.push_back(Form("Result: %.4f #pm %.4f keV", mean, total));
  labels.push_back("ENSDF 68.752 keV");
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  const Double_t margin = 0.1 * (y_hi - y_lo);
  // Room for the five-digit keV labels beside the axis title.
  plot->SetLeftMargin(0.17);
  for (Int_t d = 0; d < N_GE_DAYS; d++) {
    PlottingUtils::ConfigureGraph(days[d], colors[d],
                                  ";Rise time 10-90% [ns];E_{#gamma} [keV]");
    days[d]->GetYaxis()->SetTitleOffset(1.5);
    days[d]->SetMarkerStyle(20 + d);
    days[d]->GetXaxis()->SetLimits(x_lo, x_hi);
    days[d]->GetYaxis()->SetRangeUser(y_lo - margin, y_hi + margin);
    days[d]->Draw(d == 0 ? "AP" : "P SAME");
    leg->AddEntry(days[d], labels[d], "pe");
  }
  TBox *band = new TBox(x_lo, mean - total, x_hi, mean + total);
  band->SetFillColorAlpha(kGray + 2, 0.2);
  band->Draw();
  TLine *line = new TLine(x_lo, mean, x_hi, mean);
  line->SetLineColor(kGray + 2);
  line->SetLineStyle(2);
  line->Draw();
  TLine *ensdf = new TLine(x_lo, ReferenceLines::GE73M_ENSDF, x_hi,
                           ReferenceLines::GE73M_ENSDF);
  ensdf->SetLineColor(kMagenta + 2);
  ensdf->SetLineStyle(3);
  ensdf->Draw();
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    days[d]->Draw("P SAME");
  leg->AddEntry(band, labels[N_GE_DAYS], "f");
  leg->AddEntry(ensdf, labels[N_GE_DAYS + 1], "l");
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
  out << "calibration " << b.calibration << std::endl;
  out << "outlier " << b.outlier << std::endl;
  out << "lineshape " << b.lineshape << std::endl;
  out << "escapes " << b.escapes << std::endl;
  out << "background " << b.background << std::endl;
  out << "pb_kbeta " << b.kb << std::endl;
  out << "curvature " << b.curvature << std::endl;
  out << "slice_range " << result.slice_range << std::endl;
  out << "slice_combination " << result.slice_combination << std::endl;
  std::cout << "Wrote " << path << std::endl;
}

} // namespace CdTeBudget

void GeBudget() {
  using namespace CdTeBudget;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/ge_budget",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;

  SourceLine lines[N_GE_DAYS][N_RISE_SLICES];
  FitSourceLines(lines);

  CellResults all;
  for (Int_t cell = 0; cell < N_CELLS; cell++) {
    const Int_t d = CellDayIndex(cell), s = CellSlice(cell);
    const TString name = CellName(cell);
    TH1D *h = FillRise(OfDay(GeRuns(), GE_DAYS[d]), Form("budget_h%d", cell),
                       RISE_SLICE_LO[s], RISE_SLICE_HI[s]);
    SetSourceLowTail(s, s);
    g_tails = TAILS_FREE;
    PeakFit capture = FitCaptureLine(h, Form("budget%d", cell));
    g_capture_meas = capture.mu;
    std::cout << "BUDGET " << name << ", capture line: sigma " << std::fixed
              << std::setprecision(3) << capture.sigma << ", low exp "
              << capture.tail[0] << " +- " << capture.tail_err[0] << " x "
              << capture.tail[1] << " +- " << capture.tail_err[1]
              << " sigma, high exp " << capture.tail[2] << " x "
              << capture.tail[3] << " sigma, chi2/ndf " << std::setprecision(2)
              << capture.chi2 << std::endl;
    all.cell.push_back(cell);
    all.x_am_err.push_back(lines[d][s].x_am_err);
    all.x_ba_err.push_back(lines[d][s].x_ba_err);
    for (Int_t v = 0; v < N_VARIANTS; v++) {
      ApplyVariant(v, s);
      TF1 *f = FitInRun(h, Form("budget%d_v%d", cell, v));
      Transferred t = GainTransfer(lines[d][s], f);
      all.value[v].push_back(t.e);
      all.error[v].push_back(t.err);
      if (v == V_NOMINAL) {
        all.d_am.push_back(t.d_am);
        all.d_ba.push_back(t.d_ba);
        all.in_run.push_back(f->GetParameter(Q_GE_E));
        all.in_run_err.push_back(f->GetParError(Q_GE_E));
      }
      if (v == V_KAPPA_FREE) {
        all.kappa.push_back(f->GetParameter(Q_KAPPA));
        all.kappa_err.push_back(f->GetParError(Q_KAPPA));
      }
      std::cout << "BUDGET " << name << ", " << VARIANT_NAMES[v] << ": Ge "
                << std::fixed << std::setprecision(4) << t.e << " +/- "
                << std::setprecision(1) << 1000 * t.err << " eV (in-run "
                << std::setprecision(4) << f->GetParameter(Q_GE_E) << " +/- "
                << std::setprecision(1) << 1000 * f->GetParError(Q_GE_E)
                << " eV), gain " << std::setprecision(5) << t.gain
                << ", Pb Ka2 " << std::showpos << std::setprecision(1)
                << 1000 * t.ka2_residual << " eV, Kb1 " << 1000 * t.kb1_residual
                << std::noshowpos << " eV, chi2/ndf " << std::setprecision(2)
                << f->GetChisquare() / f->GetNDF() << ", Cd Kb escape / Ge "
                << std::setprecision(3)
                << f->GetParameter(Q_ESC_A) / f->GetParameter(Q_GE_A)
                << ", Pb Ka2/Ka1 " << f->GetParameter(Q_KA2_R) << ", step "
                << f->GetParameter(Q_STEP) << ", low exp "
                << f->GetParameter(Q_LOW_EXP) << " x "
                << f->GetParameter(Q_LOW_EXP_DECAY) << " sigma, continuum "
                << std::setprecision(1) << Bkg(g_lo, f->GetParameters())
                << " -> " << Bkg(g_hi, f->GetParameters()) << std::endl;
      if (v == V_NOMINAL)
        DrawInRun(h, f,
                  Form("psi_cdte_ge_10%02d_slice_%.0f_%.0fns", GE_DAYS[d],
                       RISE_SLICE_LO[s], RISE_SLICE_HI[s]),
                  Form("10/%02d, rise time %.1f-%.1f ns, Ge runs", GE_DAYS[d],
                       RISE_SLICE_LO[s], RISE_SLICE_HI[s]));
      delete f;
    }
    ApplyVariant(V_NOMINAL, s);
    delete h;
  }

  Budget narrow = Summarize(Form("rise <= %.1f ns", CHECK_MAX_RISE_NS), all,
                            CHECK_MAX_RISE_NS);
  Budget b = Summarize(Form("rise < %.0f ns", BUDGET_MAX_RISE_NS), all,
                       BUDGET_MAX_RISE_NS);
  FinalResult result = Finalize(all, b, narrow);
  PlotGeVsRise(all, b.mean, result.total);
  WriteResult(project_root + "/results/cdte_ge.result", result);
}
