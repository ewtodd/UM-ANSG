// Final-budget combiner for the Ge-73m gamma energy.
//
// Reads the per-run rows written by CalibrationLow.cpp:
//   method(insitu|ratesub) ge_mu fit_err cal_err bkg_err chi2 gain_err label
//
// Method 1 (in-situ) sets the central value; Method 2 (rate-sub) enters only
// as the method systematic. Runs are combined by BLUE with the calibration
// split into its correlated (Am reference) and independent (gain transfer)
// parts.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include <RtypesCore.h>
#include <TString.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

// Which run to combine. Must match the RESULT_TAG the rows were written with.
const TString RESULT_TAG = "amxfer";

// Lineshape variants to compare. Both .result files must exist. The nominal
// is the no-high-tail fit: that component rails at both bounds when enabled,
// a flat pedestal degenerate with the background rather than a model the
// data supports.
const TString PRIMARY_VARIANT = "_nohitail";
const TString ALTERNATE_VARIANT = "_hitail";

// Alternate lineshape variant (high-side exponential tail). OFF: on Am-241
// that component fits to 6.8e-09 +/- 1.1e-02, an unconstrainable null
// parameter, and enabling it degrades the fits (chi2/ndf 2.21 vs 0.57,
// errors swinging 2-3x on identical data). Its separation from the nominal
// is not a systematic.
const Bool_t USE_LINESHAPE_VARIANT = kFALSE;

// Fit-quality cut on the in-situ fits. chi2/ndf carries no information about
// where mu landed, so it cannot bias the result. Taken on the UNION of
// failures across both lineshape variants so both are combined over the same
// runs. A rejected run keeps its rate-sub row: this removes a failed fit,
// not a dataset.
const Double_t INSITU_CHI2_MAX = 3.0;

struct Row {
  TString method;
  TString label;
  Double_t mu = 0, fit_err = 0, cal_err = 0, bkg_err = 0, chi2 = 0;
  Double_t gain_err = 0;
};

std::vector<Row> ReadResult(const TString &path) {
  std::vector<Row> rows;
  std::ifstream in(path.Data());
  if (!in.is_open()) {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return rows;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream ss(line);
    Row r;
    std::string method, label;
    // Columns: method mu fit_err cal_err bkg_err chi2 gain_err label
    Double_t mu, fe, ce, be, c2, ge;
    if (!(ss >> method >> mu >> fe >> ce >> be >> c2 >> ge >> label))
      continue;
    r.method = method.c_str();
    r.mu = mu;
    r.fit_err = fe;
    r.cal_err = ce;
    r.bkg_err = be;
    r.chi2 = c2;
    r.gain_err = ge;
    r.label = label.c_str();
    rows.push_back(r);
  }
  return rows;
}

Double_t WeightedMean(const std::vector<Double_t> &v,
                      const std::vector<Double_t> &e) {
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < v.size(); i++) {
    if (e[i] <= 0)
      continue;
    Double_t w = 1.0 / (e[i] * e[i]);
    sw += w;
    swx += w * v[i];
  }
  return (sw > 0) ? swx / sw : 0;
}

// Combined estimate; err already CONTAINS the calibration, never add it again.
struct Blue {
  Double_t mean = 0, err = 0, err_nocal = 0, cal_part = 0, chi2ndf = 0;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

Bool_t InvertSym(std::vector<std::vector<Double_t>> &m) {
  Int_t n = (Int_t)m.size();
  std::vector<std::vector<Double_t>> a(n, std::vector<Double_t>(2 * n, 0.0));
  for (Int_t i = 0; i < n; i++) {
    for (Int_t j = 0; j < n; j++)
      a[i][j] = m[i][j];
    a[i][n + i] = 1.0;
  }
  for (Int_t c = 0; c < n; c++) {
    Int_t piv = c;
    for (Int_t r = c + 1; r < n; r++)
      if (std::fabs(a[r][c]) > std::fabs(a[piv][c]))
        piv = r;
    if (std::fabs(a[piv][c]) < 1e-300)
      return kFALSE;
    std::swap(a[c], a[piv]);
    Double_t d = a[c][c];
    for (Int_t j = 0; j < 2 * n; j++)
      a[c][j] /= d;
    for (Int_t r = 0; r < n; r++) {
      if (r == c)
        continue;
      Double_t f = a[r][c];
      if (f == 0.0)
        continue;
      for (Int_t j = 0; j < 2 * n; j++)
        a[r][j] -= f * a[c][j];
    }
  }
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = 0; j < n; j++)
      m[i][j] = a[i][n + j];
  return kTRUE;
}

// BLUE with an explicit covariance V_ij = diag_i^2 d_ij + sum_k c_k,i c_k,j,
// each c_k a fully correlated component. Two here: the Am reference every
// run inherits, and the gain transfer, correlated through that same
// reference. Only fit errors stay on the diagonal.
Blue CombineBlueGen(const std::vector<Double_t> &x,
                    const std::vector<Double_t> &diag,
                    const std::vector<std::vector<Double_t>> &corr) {
  Blue b;
  Int_t n = (Int_t)x.size();
  if (n == 0)
    return b;
  std::vector<std::vector<Double_t>> V(n, std::vector<Double_t>(n, 0.0));
  for (Int_t i = 0; i < n; i++) {
    V[i][i] += diag[i] * diag[i];
    for (Int_t j = 0; j < n; j++)
      for (size_t k = 0; k < corr.size(); k++)
        V[i][j] += corr[k][i] * corr[k][j];
  }
  if (!InvertSym(V))
    return b;
  Double_t A = 0, B = 0;
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = 0; j < n; j++) {
      A += V[i][j];
      B += V[i][j] * x[j];
    }
  if (!(A > 0))
    return b;
  b.n = n;
  b.mean = B / A;
  Double_t chi2 = 0;
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = 0; j < n; j++)
      chi2 += (x[i] - b.mean) * V[i][j] * (x[j] - b.mean);
  b.chi2ndf = (n > 1) ? chi2 / (n - 1) : 0;
  Double_t scale = (b.chi2ndf > 1.0) ? std::sqrt(b.chi2ndf) : 1.0;
  b.err = std::sqrt(1.0 / A) * scale;
  // Diagnostic: the same estimate with every correlated block switched off, so
  // the correlated share of the total is visible rather than asserted.
  Double_t S0 = 0;
  for (Int_t i = 0; i < n; i++)
    if (diag[i] > 0)
      S0 += 1.0 / (diag[i] * diag[i]);
  b.err_nocal = (S0 > 0) ? std::sqrt(1.0 / S0) * scale : 0.0;
  b.cal_part = (b.err > b.err_nocal)
                   ? std::sqrt(b.err * b.err - b.err_nocal * b.err_nocal)
                   : 0.0;
  b.valid = kTRUE;
  return b;
}

// In-situ labels that pass the fit-quality cut in BOTH variants. A run cut in
// one variant and kept in the other makes the lineshape systematic a comparison
// of different DATASETS rather than of two models of the same data, so the cut
// must be taken on the union of failures.
std::set<TString> KeptInsituLabels(const std::vector<Row> &pri,
                                   const std::vector<Row> &alt) {
  std::set<TString> failed, all_insitu;
  const std::vector<Row> *srcs[2] = {&pri, &alt};
  for (Int_t s = 0; s < 2; s++)
    for (size_t i = 0; i < srcs[s]->size(); i++) {
      const Row &r = (*srcs[s])[i];
      if (r.method != "insitu")
        continue;
      all_insitu.insert(r.label);
      if (r.chi2 > INSITU_CHI2_MAX)
        failed.insert(r.label);
    }
  std::set<TString> kept;
  for (std::set<TString>::const_iterator it = all_insitu.begin();
       it != all_insitu.end(); ++it)
    if (failed.find(*it) == failed.end())
      kept.insert(*it);
  for (std::set<TString>::const_iterator it = failed.begin();
       it != failed.end(); ++it)
    std::cout << "  [cut] in-situ row " << *it << " fails chi2/ndf > "
              << INSITU_CHI2_MAX
              << " in at least one lineshape variant; dropped from BOTH"
              << std::endl;
  return kept;
}

// One lineshape variant, combined over the rows the mode selects.
// mode 0 = in-situ only, 1 = in-situ + rate-sub, 2 = rate-sub only.
// err ALREADY CONTAINS the calibration; never add a calibration term again.
struct Variant {
  Double_t mean = 0, err = 0, err_nocal = 0, cal_part = 0, chi2ndf = 0;
  Int_t n = 0;
};

Variant BuildVariant(const std::vector<Row> &rows, Int_t mode,
                     const std::set<TString> &kept_insitu) {
  // Three error sources per run: fit (independent, diagonal), gain transfer
  // (all measured against the same reference centroids -> correlated), and
  // the Am reference itself, sqrt(cal^2 - gain^2) (shared -> correlated).
  // The 15th/16th self-calibrate and inherit no reference, so their
  // calibration goes on the DIAGONAL rather than in a shared block.
  std::vector<Double_t> mu, diag, c_ref, c_xfer;
  for (size_t i = 0; i < rows.size(); i++) {
    const Row &r = rows[i];
    if (r.method == "insitu") {
      if (mode == 2 || kept_insitu.find(r.label) == kept_insitu.end())
        continue;
    } else if (r.method == "ratesub") {
      if (mode == 0)
        continue;
    } else {
      continue;
    }
    Bool_t am_ref_group =
        r.label.Contains("(01/13)") || r.label.Contains("(01/14)");
    Double_t g = r.gain_err;
    Double_t ref2 = r.cal_err * r.cal_err - g * g;
    Double_t ref = (ref2 > 0) ? std::sqrt(ref2) : 0.0;
    mu.push_back(r.mu);
    if (am_ref_group) {
      diag.push_back(r.fit_err);
      c_ref.push_back(ref);
      c_xfer.push_back(g);
    } else {
      // Own calibration: fold it into the diagonal, contribute to neither
      // shared block.
      diag.push_back(std::sqrt(r.fit_err * r.fit_err + r.cal_err * r.cal_err));
      c_ref.push_back(0.0);
      c_xfer.push_back(0.0);
    }
  }
  Variant v;
  if (mu.empty())
    return v;
  std::vector<std::vector<Double_t>> corr;
  corr.push_back(c_ref);
  corr.push_back(c_xfer);
  Blue b = CombineBlueGen(mu, diag, corr);
  v.mean = b.mean;
  v.err = b.err;
  v.err_nocal = b.err_nocal;
  v.cal_part = b.cal_part;
  v.chi2ndf = b.chi2ndf;
  v.n = b.n;
  return v;
}

// The two lineshape variants averaged. They are the SAME DATA refit with a
// different model, so they are ~100% correlated: the central value is their
// uncertainty-weighted average, but the base error is the weighted MEAN of the
// two errors, not 1/sqrt(sum 1/sigma^2), which would return less than either
// input. Their separation is the lineshape systematic.
struct Averaged {
  Variant p, a;
  Double_t central = 0, base = 0, lineshape = 0;
  Bool_t have_alt = kFALSE, valid = kFALSE;
};

Averaged AverageVariants(const std::vector<Row> &pri,
                         const std::vector<Row> &alt, Int_t mode,
                         const std::set<TString> &kept) {
  Averaged o;
  o.p = BuildVariant(pri, mode, kept);
  if (o.p.n == 0)
    return o;
  o.a = BuildVariant(alt, mode, kept);
  o.valid = kTRUE;
  if (o.a.n > 0 && o.a.err > 0 && o.p.err > 0) {
    Double_t wp = 1.0 / (o.p.err * o.p.err);
    Double_t wa = 1.0 / (o.a.err * o.a.err);
    o.central = (o.p.mean * wp + o.a.mean * wa) / (wp + wa);
    o.base = (o.p.err * wp + o.a.err * wa) / (wp + wa);
    o.lineshape = std::fabs(o.p.mean - o.a.mean);
    o.have_alt = kTRUE;
  } else {
    o.central = o.p.mean;
    o.base = o.p.err;
  }
  return o;
}

struct Scheme {
  Double_t central = 0, base = 0, lineshape = 0, method = 0, total = 0;
  Bool_t valid = kFALSE;
};

// Print one scheme's FULL derivation, so the quoted number and the printed
// budget come from the same arithmetic and the column can be added up.
Scheme ReportScheme(const TString &title, const std::vector<Row> &pri,
                    const std::vector<Row> &alt, Int_t mode,
                    const std::set<TString> &kept, Double_t ratesub_central,
                    Bool_t have_ratesub) {
  Scheme s;
  Averaged o = AverageVariants(pri, alt, mode, kept);
  if (!o.valid)
    return s;

  std::cout << std::endl;
  std::cout << "========== " << title << " ==========" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  const Variant *vs[2] = {&o.p, &o.a};
  const TString names[2] = {PRIMARY_VARIANT, ALTERNATE_VARIANT};
  for (Int_t k = 0; k < 2; k++) {
    if (vs[k]->n == 0)
      continue;
    std::cout << "  " << std::left << std::setw(11) << names[k] << std::right
              << " n=" << vs[k]->n << "  mean " << vs[k]->mean << "  err "
              << vs[k]->err << "  (indep " << vs[k]->err_nocal
              << " (+) common-cal " << vs[k]->cal_part << ")  chi2/ndf "
              << std::setprecision(2) << vs[k]->chi2ndf << std::setprecision(5)
              << std::endl;
  }

  s.central = o.central;
  s.base = o.base;
  s.lineshape = o.lineshape;
  s.method = have_ratesub ? std::fabs(o.central - ratesub_central) : 0;
  s.total = std::sqrt(s.base * s.base + s.lineshape * s.lineshape +
                      s.method * s.method);
  s.valid = kTRUE;

  std::cout << "  central   = " << s.central
            << "   (variant means, weighted by 1/err^2)" << std::endl;
  std::cout << "  base      = " << s.base
            << "   (weighted MEAN of the two variant errors -- same data)"
            << std::endl;
  std::cout << "  lineshape = " << s.lineshape << "   (variant separation)"
            << std::endl;
  if (have_ratesub)
    std::cout << "  method    = " << s.method
              << "   (MEASURED gap |M1 - M2|, M2 = " << ratesub_central << ")"
              << std::endl;
  std::cout << "  TOTAL     = sqrt(" << s.base << "^2 + " << s.lineshape
            << "^2 + " << s.method << "^2) = " << s.total << std::endl;
  std::cout << std::setprecision(4) << "  E = " << s.central << " +/- "
            << s.total << " keV" << std::endl;
  return s;
}

void CombineGeResult() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  // PRIMARY is the nominal lineshape; ALTERNATE is the same analysis with the
  // high-side exponential tail toggled. Their difference is the lineshape
  // systematic -- the two-file mechanism this combiner was written around,
  // restored with a better-motivated pair than the old shared/ka1only split.
  std::vector<Row> shared = ReadResult(
      project_root + "/results/ge_" + RESULT_TAG + PRIMARY_VARIANT + ".result");
  if (shared.empty()) {
    std::cerr << "No results for variant '" << PRIMARY_VARIANT
              << "'; run CalibrationLow.cpp first." << std::endl;
    return;
  }
  std::vector<Row> alt;
  if (USE_LINESHAPE_VARIANT)
    alt = ReadResult(project_root + "/results/ge_" + RESULT_TAG +
                     ALTERNATE_VARIANT + ".result");
  if (USE_LINESHAPE_VARIANT && alt.empty())
    std::cerr << "WARNING: no results for alternate variant '"
              << ALTERNATE_VARIANT
              << "'; the lineshape systematic will be ABSENT from the budget, "
                 "not zero. Run CalibrationLow.cpp with USE_HIGH_EXP_TAIL "
                 "flipped."
              << std::endl;

  // bkg_err is |cal(precal Ge) - postcal Ge|, a diagnostic only: with the
  // unit bug fixed the post-cal fit returns its own seed, so it is a
  // tautology rather than a systematic. Printed, never in the budget.
  std::vector<Double_t> all_be, all_fe;
  for (size_t i = 0; i < shared.size(); i++) {
    all_be.push_back(shared[i].bkg_err);
    all_fe.push_back(shared[i].fit_err);
  }
  Double_t background_sys = WeightedMean(all_be, all_fe);
  std::cout << std::fixed << std::setprecision(4);

  // Fit-quality cut, taken on the UNION of failures across the two lineshape
  // variants so both are combined over the same runs.
  std::set<TString> kept = KeptInsituLabels(shared, alt);

  // Rate-sub central value, averaged over the two variants exactly as the
  // in-situ value is, so the method systematic compares like with like. Taking
  // it from the primary variant alone would fold the lineshape shift into the
  // method term and count the same effect twice.
  Averaged rs = AverageVariants(shared, alt, 2, kept);

  // TWO INDEPENDENT MEASUREMENTS OF THE SAME QUANTITY, on disjoint datasets and
  // with entirely different response functions, backgrounds and treatments of
  // the Pb lines. Their agreement is the analysis's main internal check, and
  // their DIFFERENCE is exactly what is quoted as the method systematic -- the
  // budget line is not an assumed model spread, it is this measured gap.
  std::cout << std::endl;
  std::cout << "===== Two independent measurements of E(Ge-73m) ====="
            << std::endl;
  Averaged in_only = AverageVariants(shared, alt, 0, kept);
  if (in_only.valid)
    std::cout << "  Method 1  in-situ simultaneous fit    (" << in_only.p.n
              << " datasets, 01/13)   : " << in_only.central << " +/- "
              << in_only.base << " keV" << std::endl;
  if (rs.valid)
    std::cout << "  Method 2  live-time rate subtraction  (" << rs.p.n
              << " datasets, 01/14-16): " << rs.central << " +/- " << rs.base
              << " keV" << std::endl;
  if (in_only.valid && rs.valid) {
    Double_t diff = std::fabs(in_only.central - rs.central);
    Double_t ddiff = std::sqrt(in_only.base * in_only.base + rs.base * rs.base);
    std::cout << "  ---> they agree to " << diff * 1000.0 << " eV, against "
              << ddiff * 1000.0 << " eV on the difference ("
              << std::setprecision(2) << (ddiff > 0 ? diff / ddiff : 0)
              << " sigma)" << std::setprecision(4) << std::endl;
    std::cout << "       The QUOTED value is Method 1. That " << diff * 1000.0
              << " eV gap IS the method systematic below -- measured, not"
              << std::endl;
    std::cout << "       assumed. It is a null result: the comparison resolves "
              << "method dependence" << std::endl;
    std::cout << "       only to ~" << ddiff * 1000.0
              << " eV, so it shows no bias rather than bounding one at "
              << diff * 1000.0 << " eV." << std::endl;
  }
  std::cout << "  [diag] |cal(pre)-post| = " << background_sys
            << "  (post-cal drift from seed; NOT in any budget)" << std::endl;

  // Two schemes, each printed with its full derivation so the column adds up.
  //
  // A is what this analysis quotes: the central value comes from the in-situ
  // runs and the rate-sub enters only as the method systematic. B additionally
  // gives the rate-sub weight in the central value -- arguable, but it
  // double-uses those runs, since their difference from the in-situ value is
  // still charged as a systematic. Both are printed; A is quoted.
  Scheme sA = ReportScheme("Scheme A: in-situ only (rate-sub as systematic)",
                           shared, alt, 0, kept, rs.central, rs.valid);
  Scheme sB = ReportScheme("Scheme B: in-situ + rate-sub, both weighted",
                           shared, alt, 1, kept, rs.central, rs.valid);
  Double_t sys = std::sqrt(sA.method * sA.method + sA.lineshape * sA.lineshape);
  Double_t total = sA.valid ? sA.total : 0.0;
  Double_t combined_mean = sA.valid ? sA.central : 0.0;

  std::cout << std::endl;
  std::cout << "========== FINAL: Ge-73m gamma energy ==========" << std::endl;
  std::cout << "  QUOTED = Method 1 (in-situ); Method 2 enters ONLY as the"
            << " method systematic" << std::endl;
  std::cout << "  E = " << combined_mean << " +/- " << sA.base
            << " (stat+cal) +/- " << sys << " (sys) keV" << std::endl;
  std::cout << "  E = " << combined_mean << " +/- " << total
            << " keV  (combined)" << std::endl;
  if (sB.valid)
    std::cout << "  [alt] Scheme B = " << sB.central << " +/- " << sB.total
              << " keV" << std::endl;
  std::cout << "  literature 68.752(7), collaborator 68.755" << std::endl;

  // Comparison table in the collaborator's format: Eg in keV, uncertainties in
  // eV. Columns Fit | Slope | Offset | Bkg | Comb. -- the calibration term is
  // split into the orthogonal pivot slope/offset (slope^2+offset^2 = cal^2).
  std::cout << std::endl;
  std::cout << "----- Comparison table (Eg keV, uncertainties eV) -----"
            << std::endl;
  std::cout << "  (cal = full calibration-curve covariance propagated to "
               "E(mu_Ge); bkg = |cal(precal Ge) - postcal Ge|)"
            << std::endl;
  std::cout << std::left << std::setw(34) << "Run" << std::right
            << std::setw(10) << "Eg[keV]" << std::setw(7) << "Fit"
            << std::setw(8) << "Cal" << std::setw(7) << "Bkg" << std::setw(9)
            << "Comb." << std::endl;
  for (size_t i = 0; i < shared.size(); i++) {
    if (shared[i].method != "insitu" && shared[i].method != "ratesub")
      continue;
    Double_t fit_ev = shared[i].fit_err * 1000.0;
    Double_t cal_ev = shared[i].cal_err * 1000.0;
    // Per-run background = |cal(precal Ge) - postcal Ge| for this run.
    Double_t bkg_ev = shared[i].bkg_err * 1000.0;
    Double_t comb_ev =
        std::sqrt(fit_ev * fit_ev + cal_ev * cal_ev + bkg_ev * bkg_ev);
    std::cout << std::left << std::setw(34) << shared[i].label << std::right
              << std::fixed << std::setprecision(4) << std::setw(10)
              << shared[i].mu << std::setprecision(1) << std::setw(7) << fit_ev
              << std::setw(8) << cal_ev << std::setw(7) << bkg_ev
              << std::setw(9) << comb_ev << std::endl;
  }
  // Fit / Cal here are the independent and correlated halves of the quoted
  // base error; Comb. is the quoted total. bkg is shown but not summed.
  Double_t stat_ev = in_only.valid ? in_only.p.err_nocal * 1000.0 : 0;
  Double_t cal_ev = in_only.valid ? in_only.p.cal_part * 1000.0 : 0;
  Double_t bkg_ev = (background_sys > 0) ? background_sys * 1000.0 : 0;
  Double_t comb_ev = total * 1000.0;
  std::cout << std::left << std::setw(34) << "COMBINED (this work, CZT)"
            << std::right << std::fixed << std::setprecision(4) << std::setw(10)
            << combined_mean << std::setprecision(1) << std::setw(7) << stat_ev
            << std::setw(8) << cal_ev << std::setw(7) << bkg_ev << std::setw(9)
            << comb_ev << std::endl;
}
