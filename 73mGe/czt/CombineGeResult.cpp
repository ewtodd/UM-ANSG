// Final-budget combiner for the Ge-73m gamma energy.
//
// Reads the per-run rows written by CalibrateAndFit.cpp:
//   method(insitu|ratesub) ge_mu fit_err cal_err bkg_err chi2 gain_err
//   cal_off sys_err anchor_err label
//
// The quoted value is the BLUE of every run of both methods: Method 1
// (in-situ, 01/13) and Method 2 (rate subtraction, 01/15-16, with its own
// systematics in sys_err: bin phase, residual-background shape, calibration
// transfer). The calibration is split into its correlated (Am reference) and
// independent (gain transfer) parts, and the transfer model (gain against
// offset, from a second result set) is a third correlated part. The 95.902 keV
// escape correction (from a result set without it) is charged in full, in two
// more correlated parts, one per set of escape fractions. On top comes the
// published budget's outlier removal, in quadrature: the largest shift of the
// BLUE without its highest run, its lowest, or both. Method 1 alone and the
// uncorrelated weighted mean are printed for comparison.
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

// Fit-quality cut on the in-situ fits. chi2/ndf carries no information about
// where mu landed, so it cannot bias the result. A rejected run drops out of
// Method 1; no in-situ run also has a rate-sub row.
const Double_t INSITU_CHI2_MAX = 3.0;

// Combine the quoted CZT value with the PSI CdTe one (results/cdte_ge.result,
// from cdte/GeBudget.cpp) in an error-weighted mean. Off until the CdTe
// analysis is final.
const Bool_t COMBINE_WITH_CDTE = kFALSE;

struct Row {
  TString method;
  TString label;
  Double_t mu = 0, fit_err = 0, cal_err = 0, bkg_err = 0, chi2 = 0;
  Double_t gain_err = 0, cal_off = 0, sys_err = 0;
  // The signal-run Pb anchor's offset error, part of cal_err but independent
  // between runs.
  Double_t anchor_err = 0;
  // The row's Ge under the alternative transfer model minus this one; one
  // modelling choice, so fully correlated between the rows it moves.
  Double_t model_shift = 0;
  // The 95.902 keV escapes' correction, this row's Ge with them minus
  // without, charged in full as a systematic.
  Double_t escape_shift = 0;
};

// The Cd-shield runs' escapes come from photons out of the shield, with
// Geant4's fractions; every other run's from captures inside the crystal,
// with the fractions measured in the 01/15-16 background runs. Each set's
// correction is fully correlated among its runs and independent of the
// other's.
Bool_t EscapesFromShield(const Row &r) { return r.label.Contains("Cd_Shield"); }

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
    // method, then the numeric columns, then the label; cal_off, sys_err and
    // anchor_err were added later and are read as 0 from a file that predates
    // them.
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
    r.sys_err = (tok.size() >= 10) ? std::atof(tok[8].c_str()) : 0.0;
    r.anchor_err = (tok.size() >= 11) ? std::atof(tok[9].c_str()) : 0.0;
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
  Double_t mean = 0, err = 0, chi2ndf = 0;
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
  Double_t mean = 0, err = 0, chi2ndf = 0;
  // Exact decomposition of err^2 with the BLUE weights:
  // stat^2 + offset^2 + mult^2 + sys^2 + model^2 + escape^2 = err^2.
  Double_t stat_part = 0, off_part = 0, mult_part = 0, sys_part = 0;
  Double_t model_part = 0, escape_part = 0;
  // Each combined run: its label, BLUE weight, Ge and total error (the
  // square root of its covariance diagonal).
  std::vector<TString> labels;
  std::vector<Double_t> weights, values, total_errs;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

Combination Combine(const std::vector<Row> &rows, Int_t mode,
                    const std::set<TString> &kept_insitu) {
  // Error sources per run: fit (independent, diagonal), gain transfer (all
  // measured against the same reference centroids -> correlated), the
  // signal-run Pb anchor (its own run's Pb fit -> independent, diagonal), and
  // the Am reference itself, sqrt(cal^2 - gain^2 - anchor^2) (shared ->
  // correlated).
  // The 15th/16th self-calibrate and inherit no reference, so their
  // calibration goes on the DIAGONAL rather than in a shared block, as do
  // their own systematics (bin phase, residual background, transfer), which
  // are specific to each day.
  // The transfer-model shift is a third correlated block, and the escape
  // correction two more, one per set of escape fractions.
  std::vector<Double_t> mu, diag, c_ref, c_xfer, c_model, c_esc_inside,
      c_esc_shield;
  // Per-row offset / multiplicative pieces of the calibration, kept apart from
  // the BLUE inputs so err can be decomposed with the weights afterwards. Rows
  // in the shared block sum coherently; own-calibration rows in quadrature.
  std::vector<Double_t> fit, off, lev, own_off, own_lev, sys;
  std::vector<TString> labels;
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
    Bool_t am_ref_group = r.label.Contains("(01/13)");
    Double_t g = r.gain_err;
    Double_t a = r.anchor_err;
    Double_t ref2 = r.cal_err * r.cal_err - g * g - a * a;
    Double_t ref = (ref2 > 0) ? std::sqrt(ref2) : 0.0;
    Double_t o = std::min(r.cal_off, ref);
    Double_t l = std::sqrt(std::max(0.0, ref * ref - o * o));
    mu.push_back(r.mu);
    labels.push_back(r.label);
    fit.push_back(r.fit_err);
    sys.push_back(r.sys_err);
    c_model.push_back(r.model_shift);
    c_esc_inside.push_back(EscapesFromShield(r) ? 0.0 : r.escape_shift);
    c_esc_shield.push_back(EscapesFromShield(r) ? r.escape_shift : 0.0);
    if (am_ref_group) {
      diag.push_back(
          std::sqrt(r.fit_err * r.fit_err + r.sys_err * r.sys_err + a * a));
      c_ref.push_back(ref);
      c_xfer.push_back(g);
      off.push_back(o);
      lev.push_back(l);
      own_off.push_back(a);
      own_lev.push_back(0.0);
    } else {
      // Own calibration: fold it into the diagonal, contribute to neither
      // shared block.
      diag.push_back(std::sqrt(r.fit_err * r.fit_err + r.cal_err * r.cal_err +
                               r.sys_err * r.sys_err));
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
  corr.push_back(c_model);
  corr.push_back(c_esc_inside);
  corr.push_back(c_esc_shield);
  Blue b = CombineBlueGen(mu, diag, corr);
  v.mean = b.mean;
  v.err = b.err;
  v.chi2ndf = b.chi2ndf;
  v.n = b.n;
  v.valid = b.valid;
  for (Int_t i = 0; i < b.n && b.valid; i++) {
    Double_t var = diag[i] * diag[i];
    for (size_t k = 0; k < corr.size(); k++)
      var += corr[k][i] * corr[k][i];
    v.labels.push_back(labels[i]);
    v.weights.push_back(b.w[i]);
    v.values.push_back(mu[i]);
    v.total_errs.push_back(std::sqrt(var));
  }
  if (b.valid) {
    Double_t stat2 = 0, so = 0, sl = 0, sx = 0, own_o2 = 0, own_l2 = 0;
    Double_t sys2 = 0, sm = 0, se_inside = 0, se_shield = 0;
    for (Int_t i = 0; i < b.n; i++) {
      Double_t wi = b.w[i];
      stat2 += wi * wi * fit[i] * fit[i];
      sys2 += wi * wi * sys[i] * sys[i];
      sm += wi * c_model[i];
      se_inside += wi * c_esc_inside[i];
      se_shield += wi * c_esc_shield[i];
      so += wi * off[i];
      sl += wi * lev[i];
      sx += wi * c_xfer[i];
      own_o2 += wi * wi * own_off[i] * own_off[i];
      own_l2 += wi * wi * own_lev[i] * own_lev[i];
    }
    v.stat_part = std::sqrt(stat2) * b.scale;
    v.off_part = std::sqrt(so * so + own_o2) * b.scale;
    v.mult_part = std::sqrt(sl * sl + sx * sx + own_l2) * b.scale;
    v.sys_part = std::sqrt(sys2) * b.scale;
    v.model_part = std::fabs(sm) * b.scale;
    v.escape_part =
        std::sqrt(se_inside * se_inside + se_shield * se_shield) * b.scale;
  }
  return v;
}

Bool_t UsedInBoth(const Row &r, const std::set<TString> &kept) {
  if (r.method == "ratesub")
    return kTRUE;
  return r.method == "insitu" && kept.find(r.label) != kept.end();
}

// The published budget's robustness term: the both-methods BLUE recomputed
// without its highest run, without its lowest, and without both; the largest
// shift from the full combination is the nuisance.
Double_t OutlierNuisance(const std::vector<Row> &rows,
                         const std::set<TString> &kept, Double_t central) {
  Int_t i_hi = -1, i_lo = -1;
  for (size_t i = 0; i < rows.size(); i++) {
    if (!UsedInBoth(rows[i], kept))
      continue;
    if (i_hi < 0 || rows[i].mu > rows[i_hi].mu)
      i_hi = (Int_t)i;
    if (i_lo < 0 || rows[i].mu < rows[i_lo].mu)
      i_lo = (Int_t)i;
  }
  if (i_hi < 0 || i_lo < 0)
    return 0;
  const char *names[3] = {"without the highest", "without the lowest",
                          "without both"};
  Double_t worst = 0;
  for (Int_t k = 0; k < 3; k++) {
    std::vector<Row> subset;
    for (size_t i = 0; i < rows.size(); i++) {
      Bool_t drop_hi = (k != 1) && (Int_t)i == i_hi;
      Bool_t drop_lo = (k != 0) && (Int_t)i == i_lo;
      if (!drop_hi && !drop_lo)
        subset.push_back(rows[i]);
    }
    Combination c = Combine(subset, 1, kept);
    if (!c.valid)
      continue;
    Double_t delta = c.mean - central;
    TString dropped = (k == 0)   ? rows[i_hi].label
                      : (k == 1) ? rows[i_lo].label
                                 : rows[i_hi].label + ", " + rows[i_lo].label;
    std::cout << "  " << names[k] << " (" << dropped << ") = " << std::fixed
              << std::setprecision(5) << c.mean
              << "   delta = " << std::setprecision(1) << delta * 1000.0
              << " eV" << std::endl;
    worst = std::max(worst, std::fabs(delta));
  }
  return worst;
}

// The paper-style cross-check: an error-weighted mean of the same runs with
// every term on the diagonal, no correlations.
Double_t UncorrelatedMean(const std::vector<Row> &rows,
                          const std::set<TString> &kept, Double_t &err) {
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < rows.size(); i++) {
    const Row &r = rows[i];
    if (!UsedInBoth(r, kept))
      continue;
    Double_t var = r.fit_err * r.fit_err + r.cal_err * r.cal_err +
                   r.sys_err * r.sys_err + r.model_shift * r.model_shift +
                   r.escape_shift * r.escape_shift;
    if (!(var > 0))
      continue;
    sw += 1.0 / var;
    swx += r.mu / var;
  }
  err = (sw > 0) ? std::sqrt(1.0 / sw) : 0;
  return (sw > 0) ? swx / sw : 0;
}

// A value of a "name value" line in a CdTe result file; -1 when absent.
Double_t CdTeResultValue(const TString &path, const std::string &name) {
  std::ifstream in(path.Data());
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream ss(line);
    std::string key;
    Double_t value = 0;
    if (ss >> key >> value && key == name)
      return value;
  }
  return -1;
}

// The CZT and CdTe values, error-weighted, taken as independent: different
// detectors, beam times, calibrations (Am and Pb K-alpha against the in-run
// Pb and escape lines) and escape inputs. What they share, the reference
// line energies and Geant4's Te Kb/Ka escape ratio, is below an eV here.
void PrintCdTeCombination(const TString &cdte_path, Double_t czt,
                          Double_t czt_err) {
  Double_t cdte = CdTeResultValue(cdte_path, "ge");
  Double_t cdte_err = CdTeResultValue(cdte_path, "total");
  std::cout << std::endl;
  std::cout << "========== CZT + CdTe (PSI), error-weighted =========="
            << std::endl;
  if (!(cdte > 0) || !(cdte_err > 0)) {
    std::cerr << "  no CdTe result in " << cdte_path
              << "; run cdte/GeBudget.cpp first." << std::endl;
    return;
  }
  Double_t w_czt = 1.0 / (czt_err * czt_err);
  Double_t w_cdte = 1.0 / (cdte_err * cdte_err);
  Double_t mean = (w_czt * czt + w_cdte * cdte) / (w_czt + w_cdte);
  Double_t err = 1.0 / std::sqrt(w_czt + w_cdte);
  Double_t diff = czt - cdte;
  Double_t diff_err = std::sqrt(czt_err * czt_err + cdte_err * cdte_err);
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  CZT  " << czt << " +/- " << czt_err << " keV  (weight "
            << std::setprecision(3) << w_czt / (w_czt + w_cdte) << ")"
            << std::endl;
  std::cout << "  CdTe " << std::setprecision(4) << cdte << " +/- " << cdte_err
            << " keV  (weight " << std::setprecision(3)
            << w_cdte / (w_czt + w_cdte) << ")" << std::endl;
  std::cout << "  they differ by " << std::setprecision(1) << diff * 1000.0
            << " +/- " << diff_err * 1000.0 << " eV (" << std::setprecision(2)
            << diff / diff_err << " sigma)" << std::endl;
  std::cout << "  E = " << std::setprecision(4) << mean << " +/- " << err
            << " keV" << std::endl;
}

// tag: the RESULT_TAG the rows were written with (results/ge_<tag>.result).
// alt_tag: the same analysis under the other transfer model (gain against
// offset, which no line inside the calibrated range can decide); each run's
// difference enters as a fully correlated systematic. "" leaves it out.
// escape_tag: the same analysis without the 95.902 keV escapes; each run's
// correction (its Ge with them minus without) is charged in full as a
// systematic, correlated as EscapesFromShield() says. "" leaves it out.
void CombineGeResult(const char *tag = "tau20_sigpb_esc",
                     const char *alt_tag = "tau20_sigpb_esc_oxfer",
                     const char *escape_tag = "tau20_sigpb") {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  std::vector<Row> rows =
      ReadResult(project_root + "/results/ge_" + TString(tag) + ".result");
  if (rows.empty()) {
    std::cerr << "No results for tag '" << tag
              << "'; run CalibrateAndFit.cpp first." << std::endl;
    return;
  }
  if (alt_tag && *alt_tag) {
    std::vector<Row> alt = ReadResult(project_root + "/results/ge_" +
                                      TString(alt_tag) + ".result");
    if (alt.empty())
      std::cerr << "WARNING: no results for transfer-model tag '" << alt_tag
                << "'; the transfer-model term is left out." << std::endl;
    for (size_t i = 0; i < rows.size(); i++)
      for (size_t j = 0; j < alt.size(); j++)
        if (alt[j].method == rows[i].method && alt[j].label == rows[i].label)
          rows[i].model_shift = alt[j].mu - rows[i].mu;
    std::cout << "Transfer model, " << alt_tag << " - " << tag << " [eV]:";
    for (size_t i = 0; i < rows.size(); i++)
      std::cout << "  " << rows[i].label << " " << std::fixed
                << std::setprecision(1) << rows[i].model_shift * 1000.0;
    std::cout << std::endl;
  }
  if (escape_tag && *escape_tag) {
    std::vector<Row> without = ReadResult(project_root + "/results/ge_" +
                                          TString(escape_tag) + ".result");
    if (without.empty())
      std::cerr << "WARNING: no results for escape tag '" << escape_tag
                << "'; the escape term is left out." << std::endl;
    for (size_t i = 0; i < rows.size(); i++)
      for (size_t j = 0; j < without.size(); j++)
        if (without[j].method == rows[i].method &&
            without[j].label == rows[i].label)
          rows[i].escape_shift = rows[i].mu - without[j].mu;
    std::cout << "Escape correction, " << tag << " - " << escape_tag
              << " [eV]:";
    for (size_t i = 0; i < rows.size(); i++)
      std::cout << "  " << rows[i].label << " " << std::fixed
                << std::setprecision(1) << rows[i].escape_shift * 1000.0;
    std::cout << std::endl;
  }

  // bkg_err is |cal(precal Ge) - postcal Ge|, how far the calibrated in-situ
  // fit moves Ge from the raw fit's value mapped through the calibration: a
  // diagnostic, printed, never in the budget.
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
  // with different calibrations, backgrounds and treatments of the Pb lines.
  // Their agreement is the analysis's main internal check.
  std::cout << std::endl;
  std::cout << "===== Two independent measurements of E(Ge-73m) ====="
            << std::endl;
  Combination in_only = Combine(rows, 0, kept);
  if (in_only.valid)
    std::cout << "  Method 1  in-situ simultaneous fit    (" << in_only.n
              << " datasets, 01/13): " << in_only.mean << " +/- " << in_only.err
              << " keV" << std::endl;
  if (rs.valid)
    std::cout << "  Method 2  live-time rate subtraction  (" << rs.n
              << " datasets, 01/15-16): " << rs.mean << " +/- " << rs.err
              << " keV" << std::endl;
  if (in_only.valid && rs.valid) {
    Double_t diff = std::fabs(in_only.mean - rs.mean);
    Double_t ddiff = std::sqrt(in_only.err * in_only.err + rs.err * rs.err);
    std::cout << "  ---> they agree to " << diff * 1000.0 << " eV, against "
              << ddiff * 1000.0 << " eV on the difference ("
              << std::setprecision(2) << (ddiff > 0 ? diff / ddiff : 0)
              << " sigma)" << std::setprecision(4) << std::endl;
  }
  std::cout << "  [diag] |cal(pre)-post| = " << background_sys
            << "  (post-cal drift from seed; NOT in any budget)" << std::endl;

  Combination both = Combine(rows, 1, kept);
  std::cout << std::endl;
  std::cout << "========== Both methods, BLUE ==========" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  if (both.valid) {
    std::cout << "  n=" << both.n << "  mean " << both.mean << "  err "
              << both.err << "  chi2/ndf " << std::setprecision(2)
              << both.chi2ndf << std::setprecision(5) << std::endl;
    std::cout << "  split of err with the BLUE weights: stat " << both.stat_part
              << " (+) cal offset " << both.off_part
              << " (+) cal multiplicative " << both.mult_part << " (+) sys "
              << both.sys_part << " (+) transfer model " << both.model_part
              << " (+) escapes " << both.escape_part << std::endl;
    std::cout
        << "  per run: BLUE weight, Ge - mean, and that in its own total error"
        << std::endl;
    for (Int_t i = 0; i < both.n; i++)
      std::cout << "    " << std::left << std::setw(46) << both.labels[i]
                << std::right << std::setprecision(3) << std::setw(7)
                << both.weights[i] << std::setprecision(1) << std::setw(8)
                << (both.values[i] - both.mean) * 1000.0 << " eV"
                << std::setprecision(2) << std::setw(7)
                << (both.values[i] - both.mean) / both.total_errs[i] << " sigma"
                << std::endl;
    std::cout << std::setprecision(5);
  }

  std::cout << std::endl;
  std::cout << "----- outlier removal (as in the published approach) -----"
            << std::endl;
  Double_t nuisance = both.valid ? OutlierNuisance(rows, kept, both.mean) : 0;
  Double_t quoted_err = std::sqrt(both.err * both.err + nuisance * nuisance);
  Double_t plain_err = 0;
  Double_t plain = UncorrelatedMean(rows, kept, plain_err);

  std::cout << std::endl;
  std::cout << "========== FINAL: Ge-73m gamma energy ==========" << std::endl;
  std::cout << std::setprecision(4);
  std::cout << "  QUOTED = BLUE of both methods' runs (+) the published "
               "budget's outlier removal"
            << std::endl;
  std::cout << "  E = " << both.mean << " +/- " << both.err << " (BLUE) +/- "
            << nuisance << " (outlier removal) keV" << std::endl;
  std::cout << "  E = " << both.mean << " +/- " << quoted_err
            << " keV  (combined)" << std::endl;
  std::cout << "  [uncorrelated weighted mean] " << plain << " +/- "
            << plain_err << " keV" << std::endl;
  if (in_only.valid)
    std::cout << "  [Method 1 alone] " << in_only.mean << " +/- " << in_only.err
              << " keV" << std::endl;
  std::cout << "  literature 68.752(7), collaborator 68.755" << std::endl;

  // Comparison table in the collaborator's format: Eg in keV, uncertainties in
  // eV. Offset^2 + Mult.^2 = cal^2, see the file-top note for the definition;
  // Sys. is the rate-sub rows' own systematic, Model the transfer-model shift
  // and Esc. the escape correction. Per run, Total = fit (+) cal (+) sys (+)
  // model (+) esc; the combined row is the quoted BLUE.
  std::cout << std::endl;
  std::cout << "----- Comparison table (Eg keV, uncertainties eV) -----"
            << std::endl;
  std::cout << "  (offset = cal error at its best-constrained energy; mult. = "
               "its growth to the Ge line + gain transfer)"
            << std::endl;
  std::cout << std::left << std::setw(46) << "Run" << std::right
            << std::setw(10) << "Eg[keV]" << std::setw(8) << "Fit"
            << std::setw(8) << "Offset" << std::setw(8) << "Mult."
            << std::setw(8) << "Sys." << std::setw(8) << "Model" << std::setw(8)
            << "Esc." << std::setw(8) << "Total" << std::endl;
  for (size_t i = 0; i < rows.size(); i++) {
    if (rows[i].method != "insitu" && rows[i].method != "ratesub")
      continue;
    Double_t fit_ev = rows[i].fit_err * 1000.0;
    Double_t cal_ev = rows[i].cal_err * 1000.0;
    Double_t sys_ev = rows[i].sys_err * 1000.0;
    Double_t model_ev = rows[i].model_shift * 1000.0;
    Double_t esc_ev = rows[i].escape_shift * 1000.0;
    // With the signal-run anchor, the offset is that run's own.
    Double_t off_ev = (rows[i].anchor_err > 0)
                          ? rows[i].anchor_err * 1000.0
                          : std::min(rows[i].cal_off * 1000.0, cal_ev);
    Double_t mult_ev =
        std::sqrt(std::max(0.0, cal_ev * cal_ev - off_ev * off_ev));
    Double_t tot_ev =
        std::sqrt(fit_ev * fit_ev + cal_ev * cal_ev + sys_ev * sys_ev +
                  model_ev * model_ev + esc_ev * esc_ev);
    std::cout << std::left << std::setw(46) << rows[i].label << std::right
              << std::fixed << std::setprecision(5) << std::setw(10)
              << rows[i].mu << std::setprecision(2) << std::setw(8) << fit_ev
              << std::setw(8) << off_ev << std::setw(8) << mult_ev
              << std::setw(8) << sys_ev << std::setw(8) << model_ev
              << std::setw(8) << esc_ev << std::setprecision(1) << std::setw(8)
              << tot_ev << std::endl;
  }
  std::cout << std::left << std::setw(46) << "COMBINED (both methods, BLUE)"
            << std::right << std::fixed << std::setprecision(5) << std::setw(10)
            << both.mean << std::setprecision(2) << std::setw(8)
            << both.stat_part * 1000.0 << std::setw(8) << both.off_part * 1000.0
            << std::setw(8) << both.mult_part * 1000.0 << std::setw(8)
            << both.sys_part * 1000.0 << std::setw(8)
            << both.model_part * 1000.0 << std::setw(8)
            << both.escape_part * 1000.0 << std::setprecision(1) << std::setw(8)
            << both.err * 1000.0 << std::endl;
  std::cout << "  quoted Total = sqrt(" << std::setprecision(1)
            << both.err * 1000.0 << "^2 + outlier removal " << nuisance * 1000.0
            << "^2) = " << quoted_err * 1000.0 << " eV" << std::endl;

  if (COMBINE_WITH_CDTE && both.valid)
    PrintCdTeCombination(project_root + "/results/cdte_ge.result", both.mean,
                         quoted_err);
}
