// Final-budget combiner for the Ge-73m gamma energy.
//
// Reads the per-run rows written by CalibrationLow.cpp:
//   method(insitu|ratesub) ge_mu fit_err cal_err bkg_err chi2 gain_err
//   cal_off label
//
// Method 1 (in-situ) sets the central value; Method 2 (rate-sub) enters only
// as the method systematic. Runs are combined by BLUE with the calibration
// split into its correlated (Am reference) and independent (gain transfer)
// parts.
//
// The calibration term is reported as offset (+) multiplicative. Offset is the
// calibration error at its best-constrained energy (cal_off, the intercept
// error at the pivot for a straight line); multiplicative is the rest: the
// growth of the calibration error from that energy to the Ge line, plus the
// gain transfer, which is a pure scale factor. offset^2 + mult^2 = cal^2 per
// run, and the combined columns are the same sums taken with the BLUE
// weights.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include <RtypesCore.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

// Which run to combine. Must match the RESULT_TAG the rows were written with.
const TString RESULT_TAG = "amxfer";

// Fit-quality cut on the in-situ fits. chi2/ndf carries no information about
// where mu landed, so it cannot bias the result. A rejected run keeps its
// rate-sub row: this removes a failed fit, not a dataset.
const Double_t INSITU_CHI2_MAX = 3.0;

struct Row {
  TString method;
  TString label;
  Double_t mu = 0, fit_err = 0, cal_err = 0, bkg_err = 0, chi2 = 0;
  Double_t gain_err = 0, cal_off = 0;
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
    std::vector<std::string> tok;
    std::string t;
    while (ss >> t)
      tok.push_back(t);
    // method, then the numeric columns, then the label; cal_off was added
    // last and is read as 0 from a file that predates it.
    if (tok.size() < 8)
      continue;
    Row r;
    r.method = tok[0].c_str();
    r.label = tok.back().c_str();
    r.mu = std::atof(tok[1].c_str());
    r.fit_err = std::atof(tok[2].c_str());
    r.cal_err = std::atof(tok[3].c_str());
    r.bkg_err = std::atof(tok[4].c_str());
    r.chi2 = std::atof(tok[5].c_str());
    r.gain_err = std::atof(tok[6].c_str());
    r.cal_off = (tok.size() >= 9) ? std::atof(tok[7].c_str()) : 0.0;
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
  Double_t scale = 1;
  std::vector<Double_t> w;
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
  b.w.assign(n, 0.0);
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = 0; j < n; j++) {
      A += V[i][j];
      B += V[i][j] * x[j];
      b.w[i] += V[i][j];
    }
  if (!(A > 0))
    return b;
  for (Int_t i = 0; i < n; i++)
    b.w[i] /= A;
  b.n = n;
  b.mean = B / A;
  Double_t chi2 = 0;
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = 0; j < n; j++)
      chi2 += (x[i] - b.mean) * V[i][j] * (x[j] - b.mean);
  b.chi2ndf = (n > 1) ? chi2 / (n - 1) : 0;
  Double_t scale = (b.chi2ndf > 1.0) ? std::sqrt(b.chi2ndf) : 1.0;
  b.scale = scale;
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

std::set<TString> KeptInsituLabels(const std::vector<Row> &rows) {
  std::set<TString> kept;
  for (size_t i = 0; i < rows.size(); i++) {
    const Row &r = rows[i];
    if (r.method != "insitu")
      continue;
    if (r.chi2 > INSITU_CHI2_MAX) {
      std::cout << "  [cut] in-situ row " << r.label << " fails chi2/ndf > "
                << INSITU_CHI2_MAX << "; dropped" << std::endl;
      continue;
    }
    kept.insert(r.label);
  }
  return kept;
}

// One combination over the rows the mode selects.
// mode 0 = in-situ only, 1 = in-situ + rate-sub, 2 = rate-sub only.
// err ALREADY CONTAINS the calibration; never add a calibration term again.
struct Combination {
  Double_t mean = 0, err = 0, err_nocal = 0, cal_part = 0, chi2ndf = 0;
  // Exact decomposition of err^2 with the BLUE weights:
  // stat^2 + offset^2 + mult^2 = err^2.
  Double_t stat_part = 0, off_part = 0, mult_part = 0;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

Combination Combine(const std::vector<Row> &rows, Int_t mode,
                    const std::set<TString> &kept_insitu) {
  // Three error sources per run: fit (independent, diagonal), gain transfer
  // (all measured against the same reference centroids -> correlated), and
  // the Am reference itself, sqrt(cal^2 - gain^2) (shared -> correlated).
  // The 15th/16th self-calibrate and inherit no reference, so their
  // calibration goes on the DIAGONAL rather than in a shared block.
  std::vector<Double_t> mu, diag, c_ref, c_xfer;
  // Per-row offset / multiplicative pieces of the calibration, kept apart from
  // the BLUE inputs so err can be decomposed with the weights afterwards. Rows
  // in the shared block sum coherently; own-calibration rows in quadrature.
  std::vector<Double_t> fit, off, lev, own_off, own_lev;
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
    Double_t o = std::min(r.cal_off, ref);
    Double_t l = std::sqrt(std::max(0.0, ref * ref - o * o));
    mu.push_back(r.mu);
    fit.push_back(r.fit_err);
    if (am_ref_group) {
      diag.push_back(r.fit_err);
      c_ref.push_back(ref);
      c_xfer.push_back(g);
      off.push_back(o);
      lev.push_back(l);
      own_off.push_back(0.0);
      own_lev.push_back(0.0);
    } else {
      // Own calibration: fold it into the diagonal, contribute to neither
      // shared block.
      diag.push_back(std::sqrt(r.fit_err * r.fit_err + r.cal_err * r.cal_err));
      c_ref.push_back(0.0);
      c_xfer.push_back(0.0);
      off.push_back(0.0);
      lev.push_back(0.0);
      own_off.push_back(o);
      own_lev.push_back(l);
    }
  }
  Combination v;
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
  v.valid = b.valid;
  if (b.valid) {
    Double_t stat2 = 0, so = 0, sl = 0, sx = 0, own_o2 = 0, own_l2 = 0;
    for (Int_t i = 0; i < b.n; i++) {
      Double_t wi = b.w[i];
      stat2 += wi * wi * fit[i] * fit[i];
      so += wi * off[i];
      sl += wi * lev[i];
      sx += wi * c_xfer[i];
      own_o2 += wi * wi * own_off[i] * own_off[i];
      own_l2 += wi * wi * own_lev[i] * own_lev[i];
    }
    v.stat_part = std::sqrt(stat2) * b.scale;
    v.off_part = std::sqrt(so * so + own_o2) * b.scale;
    v.mult_part = std::sqrt(sl * sl + sx * sx + own_l2) * b.scale;
  }
  return v;
}

struct Scheme {
  Double_t central = 0, base = 0, method = 0, total = 0;
  Bool_t valid = kFALSE;
};

// Print one scheme's FULL derivation, so the quoted number and the printed
// budget come from the same arithmetic and the column can be added up.
Scheme ReportScheme(const TString &title, const std::vector<Row> &rows,
                    Int_t mode, const std::set<TString> &kept,
                    Double_t ratesub_central, Bool_t have_ratesub) {
  Scheme s;
  Combination c = Combine(rows, mode, kept);
  if (!c.valid)
    return s;

  std::cout << std::endl;
  std::cout << "========== " << title << " ==========" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  std::cout << "  n=" << c.n << "  mean " << c.mean << "  err " << c.err
            << "  (indep " << c.err_nocal << " (+) common-cal " << c.cal_part
            << ")  chi2/ndf " << std::setprecision(2) << c.chi2ndf
            << std::setprecision(5) << std::endl;
  std::cout << "  split of err with the BLUE weights: stat " << c.stat_part
            << " (+) cal offset " << c.off_part << " (+) cal multiplicative "
            << c.mult_part << std::endl;

  s.central = c.mean;
  s.base = c.err;
  s.method = have_ratesub ? std::fabs(c.mean - ratesub_central) : 0;
  s.total = std::sqrt(s.base * s.base + s.method * s.method);
  s.valid = kTRUE;

  std::cout << "  central   = " << s.central << std::endl;
  std::cout << "  base      = " << s.base << "   (BLUE, fit (+) calibration)"
            << std::endl;
  if (have_ratesub)
    std::cout << "  method    = " << s.method
              << "   (MEASURED gap |M1 - M2|, M2 = " << ratesub_central << ")"
              << std::endl;
  std::cout << "  TOTAL     = sqrt(" << s.base << "^2 + " << s.method
            << "^2) = " << s.total << std::endl;
  std::cout << std::setprecision(4) << "  E = " << s.central << " +/- "
            << s.total << " keV" << std::endl;
  return s;
}

void CombineGeResult() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  std::vector<Row> rows =
      ReadResult(project_root + "/results/ge_" + RESULT_TAG + ".result");
  if (rows.empty()) {
    std::cerr << "No results for tag '" << RESULT_TAG
              << "'; run CalibrationLow.cpp first." << std::endl;
    return;
  }

  // bkg_err is |cal(precal Ge) - postcal Ge|, a diagnostic only: with the
  // unit bug fixed the post-cal fit returns its own seed, so it is a
  // tautology rather than a systematic. Printed, never in the budget.
  std::vector<Double_t> all_be, all_fe;
  for (size_t i = 0; i < rows.size(); i++) {
    all_be.push_back(rows[i].bkg_err);
    all_fe.push_back(rows[i].fit_err);
  }
  Double_t background_sys = WeightedMean(all_be, all_fe);
  std::cout << std::fixed << std::setprecision(4);

  std::set<TString> kept = KeptInsituLabels(rows);
  Combination rs = Combine(rows, 2, kept);

  // TWO INDEPENDENT MEASUREMENTS OF THE SAME QUANTITY, on disjoint datasets and
  // with entirely different response functions, backgrounds and treatments of
  // the Pb lines. Their agreement is the analysis's main internal check, and
  // their DIFFERENCE is exactly what is quoted as the method systematic -- the
  // budget line is not an assumed model spread, it is this measured gap.
  std::cout << std::endl;
  std::cout << "===== Two independent measurements of E(Ge-73m) ====="
            << std::endl;
  Combination in_only = Combine(rows, 0, kept);
  if (in_only.valid)
    std::cout << "  Method 1  in-situ simultaneous fit    (" << in_only.n
              << " datasets, 01/13)   : " << in_only.mean << " +/- "
              << in_only.err << " keV" << std::endl;
  if (rs.valid)
    std::cout << "  Method 2  live-time rate subtraction  (" << rs.n
              << " datasets, 01/14-16): " << rs.mean << " +/- " << rs.err
              << " keV" << std::endl;
  if (in_only.valid && rs.valid) {
    Double_t diff = std::fabs(in_only.mean - rs.mean);
    Double_t ddiff = std::sqrt(in_only.err * in_only.err + rs.err * rs.err);
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
                           rows, 0, kept, rs.mean, rs.valid);
  Scheme sB = ReportScheme("Scheme B: in-situ + rate-sub, both weighted", rows,
                           1, kept, rs.mean, rs.valid);
  Double_t total = sA.valid ? sA.total : 0.0;
  Double_t combined_mean = sA.valid ? sA.central : 0.0;

  std::cout << std::endl;
  std::cout << "========== FINAL: Ge-73m gamma energy ==========" << std::endl;
  std::cout << "  QUOTED = Method 1 (in-situ); Method 2 enters ONLY as the"
            << " method systematic" << std::endl;
  std::cout << "  E = " << combined_mean << " +/- " << sA.base
            << " (stat+cal) +/- " << sA.method << " (sys) keV" << std::endl;
  std::cout << "  E = " << combined_mean << " +/- " << total
            << " keV  (combined)" << std::endl;
  if (sB.valid)
    std::cout << "  [alt] Scheme B = " << sB.central << " +/- " << sB.total
              << " keV" << std::endl;
  std::cout << "  literature 68.752(7), collaborator 68.755" << std::endl;

  // Comparison table in the collaborator's format: Eg in keV, uncertainties in
  // eV. Columns Fit | Offset | Mult. | Total -- offset^2 + mult^2 = cal^2, see
  // the file-top note for the definition. Per run, Total = fit (+) cal; the
  // combined row's Total is the quoted number and includes the method term.
  std::cout << std::endl;
  std::cout << "----- Comparison table (Eg keV, uncertainties eV) -----"
            << std::endl;
  std::cout << "  (offset = cal error at its best-constrained energy; mult. = "
               "its growth to the Ge line + gain transfer)"
            << std::endl;
  std::cout << std::left << std::setw(46) << "Run" << std::right
            << std::setw(10) << "Eg[keV]" << std::setw(8) << "Fit"
            << std::setw(8) << "Offset" << std::setw(8) << "Mult."
            << std::setw(8) << "Total" << std::endl;
  for (size_t i = 0; i < rows.size(); i++) {
    if (rows[i].method != "insitu" && rows[i].method != "ratesub")
      continue;
    Double_t fit_ev = rows[i].fit_err * 1000.0;
    Double_t cal_ev = rows[i].cal_err * 1000.0;
    Double_t off_ev = std::min(rows[i].cal_off * 1000.0, cal_ev);
    Double_t mult_ev =
        std::sqrt(std::max(0.0, cal_ev * cal_ev - off_ev * off_ev));
    Double_t tot_ev = std::sqrt(fit_ev * fit_ev + cal_ev * cal_ev);
    std::cout << std::left << std::setw(46) << rows[i].label << std::right
              << std::fixed << std::setprecision(5) << std::setw(10)
              << rows[i].mu << std::setprecision(2) << std::setw(8) << fit_ev
              << std::setw(8) << off_ev << std::setw(8) << mult_ev
              << std::setprecision(1) << std::setw(8) << tot_ev << std::endl;
  }
  Double_t stat_ev = in_only.valid ? in_only.stat_part * 1000.0 : 0;
  Double_t off_ev = in_only.valid ? in_only.off_part * 1000.0 : 0;
  Double_t mult_ev = in_only.valid ? in_only.mult_part * 1000.0 : 0;
  std::cout << std::left << std::setw(46) << "COMBINED (01/13 in-situ, BLUE)"
            << std::right << std::fixed << std::setprecision(5) << std::setw(10)
            << combined_mean << std::setprecision(2) << std::setw(8) << stat_ev
            << std::setw(8) << off_ev << std::setw(8) << mult_ev
            << std::setprecision(1) << std::setw(8) << total * 1000.0
            << std::endl;
  std::cout << "  combined Total = sqrt(" << std::setprecision(2) << stat_ev
            << "^2 + " << off_ev << "^2 + " << mult_ev << "^2 + method "
            << sA.method * 1000.0 << "^2)" << std::endl;
}
