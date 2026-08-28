// Final-budget combiner for the Ge-73m gamma energy.
//
// Reads the single per-run .result file written by CalibrationLow.cpp:
//   results/ge_<RESULT_TAG>.result
// Each row: method(insitu|ratesub)  ge_mu  fit_err  cal_err  bkg_err  chi2
//           label
//
// The calibration uncertainty is ONE term. An orthogonal slope/offset pivot
// split used to be quoted alongside it, but it is undefined for the pol2
// calibration the Am anchor requires, and a single total is what we report.
//
// Central value and the stat/cal/method/background systematics all come from
// that one file. Correlated systematics are reported as floors, not 1/sqrt(N).
//
// HISTORY: a lineshape systematic was once taken as the per-run
// |shared - ka1only| in-situ spread, with the alternate arm produced by
// CalibrationLowKa1.cpp. That macro and its .result file were deleted and the
// systematic was dropped deliberately -- Table II carries no lineshape column.
// Nothing here computes one, and this comment no longer claims otherwise.
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

// Which CalibrationLow.cpp run to combine. Must match the RESULT_TAG that macro
// wrote: "shared" for the Pb-only scheme, "amxfer" for the Am-anchored one.
const TString RESULT_TAG = "amxfer";

// Lineshape variants to compare. Both files must exist; each is produced by
// CalibrationLow.cpp with USE_HIGH_EXP_TAIL set accordingly. The nominal is
// the no-high-tail fit: that component rails at both bounds when enabled
// (amplitude pinned at 0.5, decay at 100 keV across a 19 keV window, i.e. a
// flat pedestal degenerate with the background), which is not a model the data
// supports so much as a free parameter absorbing slack.
const TString PRIMARY_VARIANT = "_nohitail";
const TString ALTERNATE_VARIANT = "_hitail";

// Use the alternate lineshape variant at all.
//
// OFF. The alternate is the nominal model plus a high-side exponential tail,
// and that component is not something the data supports -- it is a null
// parameter. On Am-241, the cleanest isolated line in the analysis (78k events,
// no doublet partner), it fits to
//
//   HighExpTailAmplitude = 6.76e-09 +/- 1.11e-02
//   HighExpTailRatio     = 1.564    +/- 88.57
//
// an amplitude six orders of magnitude below its own error, and a decay
// constant with no effect on the likelihood and therefore no determination. It
// is also backwards physically: incomplete charge collection makes a LOW-energy
// tail, and a high-energy one would be pileup, which is negligible at these
// rates.
//
// Enabling it does not test an alternative model, it hands the fit a degree of
// freedom the data cannot constrain, and the fits degrade accordingly:
// chi2/ndf 2.21 against 0.57, fit errors swinging by factors of 2-3 on
// IDENTICAL data (7.1 vs 18.3 eV, 61.6 vs 26.6 eV), a 7.6 eV error attached to
// a chi2/ndf 6.55 fit, and central values dragged 45-60 eV low. Under the
// post-cal Pb lock it reaches chi2/ndf 6.91.
//
// A systematic is the spread between models the data supports, not the gap to
// a rejected one, so quoting that separation would be inventing 8.7 eV of
// uncertainty rather than measuring it.
//
// Nothing replaces it. Locking the post-cal Pb to its tabulated energies -- the
// obvious external-truth alternative -- moves the kept runs by +3.7/-2.8/+0.7
// eV and the combined value by 1.7 eV, with chi2 unchanged, so there is no
// model dependence there to quote either. The remaining budget is fit,
// calibration and method; method already spans a lineshape change, since
// in-situ and rate-sub use entirely different response functions, backgrounds
// and treatments of the Pb.
const Bool_t USE_LINESHAPE_VARIANT = kFALSE;

// Goodness-of-fit cut on the in-situ simultaneous fits. This is a cut on fit
// QUALITY, not on the fitted value: chi2/ndf is computed from the residuals and
// carries no information about where mu landed, so applying it cannot bias the
// result toward or away from any particular energy.
//
// It exists because the 01/14 in-situ fit fails it at 4.62 while every other
// in-situ fit sits at 1.5-2.8. That run is also the one place we can prove the
// in-situ fit is at fault rather than the calibration: its rate-sub row uses
// the SAME cal_func and the SAME cal_err, and lands 126 eV lower with
// chi2/ndf 1.59. Same data, same calibration, so the whole 126 eV is the fit.
// The in-situ fit couples the Ge centroid to the Pb doublet through the linked
// Ka shape parameters, and on that run the doublet is ~47 eV wide against a
// truth known to <1 eV.
//
// The rate-sub row of a rejected run is KEPT -- the cut removes a failed fit,
// not a dataset. Set to a large value to disable.
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

// Weight each run by its TOTAL uncertainty (fit (+) cal) rather than by the fit
// error alone. Weighting on the fit error lets a run with a small fit error
// dominate no matter how badly calibrated it is: Cd 25% has fit_err 4.7 eV
// against 22-28 eV for the other two and was carrying ~93% of the weight, while
// its cal_err of 16.1 eV -- the LARGEST of the three -- never entered the
// weighting at all.
const Bool_t WEIGHT_BY_TOTAL_ERR = kTRUE;

// Inverse-variance on the per-run total treats the calibration error as
// INDEPENDENT between runs, so it shrinks by ~sqrt(N). That is only right if
// each run's calibration is its own. Here all runs descend from ONE Am-241
// reference through the gain transfer, so a common component is shared and must
// not be averaged away. With this on, the base error is floored at the
// correlated calibration term (the fit-weighted mean of the per-run cal_err),
// which is the no-reduction limit. Off = full sqrt(N) reduction, i.e. runs
// treated as fully independent.
const Bool_t CAL_CORRELATED_FLOOR = kTRUE;

// Per-run total: fit and calibration terms in quadrature.
std::vector<Double_t> TotalErrors(const std::vector<Double_t> &fe,
                                  const std::vector<Double_t> &ce) {
  std::vector<Double_t> t(fe.size(), 0.0);
  for (size_t i = 0; i < fe.size() && i < ce.size(); i++)
    t[i] = std::sqrt(fe[i] * fe[i] + ce[i] * ce[i]);
  return t;
}

struct Combo {
  Double_t mean = 0, int_err = 0, ext_err = 0, chi2ndf = 0;
  Int_t n = 0;
};

// Inverse-variance combine on fit error; report internal and external (scatter)
// errors. ext = int * sqrt(chi2/ndf about the weighted mean).
Combo Combine(const std::vector<Double_t> &m, const std::vector<Double_t> &e) {
  Combo c;
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < m.size(); i++) {
    if (e[i] <= 0)
      continue;
    Double_t w = 1.0 / (e[i] * e[i]);
    sw += w;
    swx += w * m[i];
    c.n++;
  }
  if (sw <= 0)
    return c;
  c.mean = swx / sw;
  c.int_err = std::sqrt(1.0 / sw);
  Double_t chi2 = 0;
  for (size_t i = 0; i < m.size(); i++) {
    if (e[i] <= 0)
      continue;
    Double_t d = (m[i] - c.mean) / e[i];
    chi2 += d * d;
  }
  c.chi2ndf = (c.n > 1) ? chi2 / (c.n - 1) : 0;
  c.ext_err = c.int_err * std::sqrt(c.chi2ndf > 1 ? c.chi2ndf : 1.0);
  return c;
}

// Fit-error-weighted mean of a per-run quantity (for correlated-systematic
// floors: cal error, background spread).
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

// Best linear unbiased estimate with a CORRELATED calibration component.
//
// Each run carries an independent fit error and a calibration error that is
// largely COMMON to all runs -- they descend from one Am-241 reference through
// the gain transfer. Treating the calibration as independent lets it shrink by
// sqrt(N), which understates the result; treating the whole per-run total as
// correlated overstates it. Neither is what the old max(stat, cal) floor did
// either: that was an ad-hoc conservative guess standing in for the covariance.
//
// The covariance is
//
//   V_ij = sigma_fit,i^2 delta_ij + sigma_cal,i sigma_cal,j
//
// diagonal plus rank one, so Sherman-Morrison inverts it in closed form and no
// matrix code is needed:
//
//   1' V^-1 1 = S0 - S1^2 / (1 + S2)
//   1' V^-1 x = Sx - S1 Sxc / (1 + S2)
//
// with S0 = sum 1/d, S1 = sum c/d, S2 = sum c^2/d, Sx = sum x/d,
// Sxc = sum c x/d, d = sigma_fit^2, c = sigma_cal.
//
// rho = 1 for the calibration is the conservative end: part of each run's cal
// error is its own gain-transfer term, so the truth lies between rho = 1 and
// rho = 0. The returned error already CONTAINS the calibration; it must not be
// added again downstream.
struct Blue {
  Double_t mean = 0, err = 0, err_nocal = 0, cal_part = 0, chi2ndf = 0;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

Blue CombineBlue(const std::vector<Double_t> &x,
                 const std::vector<Double_t> &fit,
                 const std::vector<Double_t> &cal) {
  Blue b;
  Double_t S0 = 0, S1 = 0, S2 = 0, Sx = 0, Sxc = 0;
  for (size_t i = 0; i < x.size(); i++) {
    if (!(fit[i] > 0))
      continue;
    Double_t d = fit[i] * fit[i];
    S0 += 1.0 / d;
    S1 += cal[i] / d;
    S2 += cal[i] * cal[i] / d;
    Sx += x[i] / d;
    Sxc += cal[i] * x[i] / d;
    b.n++;
  }
  if (b.n == 0 || !(S0 > 0))
    return b;
  Double_t den = 1.0 + S2;
  Double_t A = S0 - S1 * S1 / den;
  if (!(A > 0))
    return b;
  b.mean = (Sx - S1 * Sxc / den) / A;

  Double_t r0 = 0, r1 = 0;
  for (size_t i = 0; i < x.size(); i++) {
    if (!(fit[i] > 0))
      continue;
    Double_t d = fit[i] * fit[i];
    Double_t r = x[i] - b.mean;
    r0 += r * r / d;
    r1 += cal[i] * r / d;
  }
  Double_t chi2 = r0 - r1 * r1 / den;
  b.chi2ndf = (b.n > 1) ? chi2 / (b.n - 1) : 0;
  // PDG-style scale factor: inflate by sqrt(chi2/ndf) when the runs scatter by
  // more than their errors admit. Never deflate when they scatter by less.
  Double_t scale = (b.chi2ndf > 1.0) ? std::sqrt(b.chi2ndf) : 1.0;
  b.err = std::sqrt(1.0 / A) * scale;
  // Diagnostics: the same estimate with the calibration switched off, so the
  // calibration's share of the total is visible rather than asserted.
  b.err_nocal = std::sqrt(1.0 / S0) * scale;
  b.cal_part = (b.err > b.err_nocal)
                   ? std::sqrt(b.err * b.err - b.err_nocal * b.err_nocal)
                   : 0.0;
  b.valid = kTRUE;
  return b;
}

// Invert a small symmetric matrix by Gauss-Jordan with partial pivoting.
// n <= 6 here, so an explicit inverse is simpler and more general than any
// rank-one shortcut and admits arbitrarily many correlated components.
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

// BLUE with an explicit covariance:
//
//   V_ij = sigma_diag,i^2 delta_ij + sum_k c_k,i c_k,j
//
// Each c_k is one FULLY CORRELATED component shared across runs. Here there are
// two: the Am reference calibration, which every run inherits from the single
// pol2 anchor fit, and the gain transfer, whose per-run ratios are all measured
// against that same reference's Pb centroids and are therefore correlated
// through it rather than independent. Only the fit errors are left on the
// diagonal.
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
// err ALREADY CONTAINS the calibration (see CombineBlue); never add cal again.
struct Variant {
  Double_t mean = 0, err = 0, err_nocal = 0, cal_part = 0, chi2ndf = 0;
  Int_t n = 0;
};

Variant BuildVariant(const std::vector<Row> &rows, Int_t mode,
                     const std::set<TString> &kept_insitu) {
  // Three error sources per run, from the .result columns:
  //   fit       independent, always on the diagonal
  //   transfer  gain_err; every transferred run measures its ratio against the
  //             SAME reference Pb centroids, so these are correlated with one
  //             another, not independent -> correlated block
  //   reference sqrt(cal^2 - gain^2); the single Am pol2 anchor that the 13th
  //             and 14th all inherit -> correlated block
  //
  // The 15th and 16th self-calibrate: they run no gain transfer and inherit no
  // Am reference (their gain_err is 0 because there is nothing to transfer, not
  // because they share one). Their calibration is their own, so it belongs on
  // the DIAGONAL. Putting it in the shared block would over-correlate them.
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

  std::cout << "\n========== " << title << " ==========" << std::endl;
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

  // Split the primary (shared) file by method.
  std::vector<Double_t> in_mu, in_fe, in_ce, in_be, rs_mu, rs_fe, rs_ce, rs_be;
  for (size_t i = 0; i < shared.size(); i++) {
    const Row &r = shared[i];
    if (r.method == "insitu") {
      in_mu.push_back(r.mu);
      in_fe.push_back(r.fit_err);
      in_ce.push_back(r.cal_err);
      in_be.push_back(r.bkg_err);
    } else if (r.method == "ratesub") {
      rs_mu.push_back(r.mu);
      rs_fe.push_back(r.fit_err);
      rs_ce.push_back(r.cal_err);
      rs_be.push_back(r.bkg_err);
    }
  }

  std::vector<Double_t> in_w =
      WEIGHT_BY_TOTAL_ERR ? TotalErrors(in_fe, in_ce) : in_fe;
  std::vector<Double_t> rs_w =
      WEIGHT_BY_TOTAL_ERR ? TotalErrors(rs_fe, rs_ce) : rs_fe;
  Combo m1 = Combine(in_mu, in_w);
  Combo m2 = Combine(rs_mu, rs_w);

  std::vector<Double_t> all_mu = in_mu, all_fe = in_fe, all_ce = in_ce,
                        all_be = in_be;
  all_mu.insert(all_mu.end(), rs_mu.begin(), rs_mu.end());
  all_fe.insert(all_fe.end(), rs_fe.begin(), rs_fe.end());
  all_ce.insert(all_ce.end(), rs_ce.begin(), rs_ce.end());

  all_be.insert(all_be.end(), rs_be.begin(), rs_be.end());
  std::vector<Double_t> all_w =
      WEIGHT_BY_TOTAL_ERR ? TotalErrors(all_fe, all_ce) : all_fe;
  Combo mc = Combine(all_mu, all_w);

  // STAT: combined fit error (larger of internal/external).
  Double_t stat = std::max(mc.int_err, mc.ext_err);

  // CAL: weighted-mean per-run calibration error (correlated -> floor). ONE
  // term: the full parameter covariance of the calibration curve propagated to
  // E(mu_Ge), which for pol2 is sum_ij mu^(i+j) cov(pi,pj).
  Double_t cal_sys = WeightedMean(all_ce, all_w);

  // METHOD: in-situ vs rate-subtraction central-value difference.
  Double_t method_sys =
      (m1.n > 0 && m2.n > 0) ? std::fabs(m1.mean - m2.mean) : 0;
  // Overwritten below against the lineshape-averaged value when both variants
  // are available.

  // BACKGROUND TERM REMOVED FROM THE BUDGET.
  //
  // bkg_err is |cal(precal Ge mu) - postcal Ge mu|. It was never a background
  // systematic: for most of its life it was measuring a UNIT BUG -- the
  // post-cal fit inherited sigma and the tail lengths in RAW units and held
  // them fixed while fitting calibrated data, so the two quantities differed by
  // the ~1% calibration derivative. With that conversion fixed, the post-cal
  // fit starts from a self-consistent seed and simply returns it, and bkg_err
  // collapses to ~0.2 eV. That is a tautology, not a measurement: it is no
  // longer an independent check of anything, so quoting it would understate the
  // budget rather than pad it.
  //
  // The column is still written and printed per run as a diagnostic -- a
  // non-zero value now means the post-cal fit genuinely moved away from its
  // seed, which is worth seeing -- but it does not enter the systematic sum.
  Double_t background_sys = WeightedMean(all_be, all_fe);

  // LINESHAPE: the in-situ central value under the two lineshape variants,
  // half-separated. Both variants fit comparably (summed chi2 differs by under
  // one unit across three runs), so neither is rejectable and the spread
  // between them is a genuine model dependence rather than a bad fit. Quoted as
  // a floor, like the other correlated terms.
  // LINESHAPE COMBINATION.
  //
  // The two variants are fits to the SAME DATA with different lineshapes, so
  // they are ~100% correlated. The central value is their uncertainty-weighted
  // average, but the base uncertainty must NOT be reduced in quadrature the way
  // independent measurements would be -- re-fitting one spectrum with a second
  // model adds no information, and 1/sqrt(sum 1/sigma^2) would return an error
  // smaller than either input. The weighted MEAN of the two totals is carried
  // instead, and the full separation between the variants is then added as the
  // lineshape systematic.
  Double_t lineshape_sys = 0;
  Bool_t have_lineshape = kFALSE;
  Double_t combined_mean = mc.mean;
  Double_t combined_base_err = std::max(mc.int_err, mc.ext_err);
  {
    std::vector<Double_t> a_mu, a_fe, a_ce;
    for (size_t i = 0; i < alt.size(); i++)
      if (alt[i].method == "insitu") {
        a_mu.push_back(alt[i].mu);
        a_fe.push_back(alt[i].fit_err);
        a_ce.push_back(alt[i].cal_err);
      }
    if (!a_mu.empty() && m1.n > 0) {
      std::vector<Double_t> a_w =
          WEIGHT_BY_TOTAL_ERR ? TotalErrors(a_fe, a_ce) : a_fe;
      Combo alt_insitu = Combine(a_mu, a_w);
      Double_t alt_stat = std::max(alt_insitu.int_err, alt_insitu.ext_err);
      Double_t alt_cal = WeightedMean(a_ce, a_w);
      Double_t alt_tot =
          WEIGHT_BY_TOTAL_ERR
              ? (CAL_CORRELATED_FLOOR ? std::max(alt_stat, alt_cal) : alt_stat)
              : std::sqrt(alt_stat * alt_stat + alt_cal * alt_cal);

      Double_t pri_stat = std::max(m1.int_err, m1.ext_err);
      Double_t pri_cal = WeightedMean(in_ce, in_w);
      Double_t pri_tot =
          WEIGHT_BY_TOTAL_ERR
              ? (CAL_CORRELATED_FLOOR ? std::max(pri_stat, pri_cal) : pri_stat)
              : std::sqrt(pri_stat * pri_stat + pri_cal * pri_cal);

      if (pri_tot > 0 && alt_tot > 0) {
        Double_t wp = 1.0 / (pri_tot * pri_tot);
        Double_t wa = 1.0 / (alt_tot * alt_tot);
        combined_mean = (m1.mean * wp + alt_insitu.mean * wa) / (wp + wa);
        // Correlated inputs: weighted mean of the totals, NOT 1/sqrt(wp+wa).
        combined_base_err = (pri_tot * wp + alt_tot * wa) / (wp + wa);
        lineshape_sys = std::fabs(m1.mean - alt_insitu.mean);
        have_lineshape = kTRUE;

        std::cout << "\n  [lineshape] " << PRIMARY_VARIANT << " = " << m1.mean
                  << " +/- " << pri_tot << "   " << ALTERNATE_VARIANT << " = "
                  << alt_insitu.mean << " +/- " << alt_tot << std::endl;
        std::cout << "  [lineshape] weighted mean = " << combined_mean
                  << "   base err = " << combined_base_err
                  << "   separation = " << lineshape_sys << std::endl;
        // Method term measured against the LINESHAPE-AVERAGED in-situ value.
        // Against the primary variant alone it would also carry the lineshape
        // shift, counting the same effect in two places.
        if (m2.n > 0)
          method_sys = std::fabs(combined_mean - m2.mean);
      }
    }
  }

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
  std::cout << "\n===== Two independent measurements of E(Ge-73m) ====="
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

  combined_mean = sA.valid ? sA.central : combined_mean;
  Double_t sys = std::sqrt(sA.method * sA.method + sA.lineshape * sA.lineshape);
  Double_t total = sA.valid ? sA.total : 0.0;
  method_sys = sA.method;
  lineshape_sys = sA.lineshape;
  combined_base_err = sA.base;

  std::cout << "\n========== FINAL: Ge-73m gamma energy =========="
            << std::endl;
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
  std::cout << "\n----- Comparison table (Eg keV, uncertainties eV) -----"
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
  Double_t stat_ev = stat * 1000.0;
  Double_t cal_ev = cal_sys * 1000.0;
  Double_t bkg_ev = (background_sys > 0) ? background_sys * 1000.0 : 0;
  Double_t comb_ev =
      std::sqrt(stat_ev * stat_ev + cal_ev * cal_ev + bkg_ev * bkg_ev);
  std::cout << std::left << std::setw(34) << "COMBINED (this work, CZT)"
            << std::right << std::fixed << std::setprecision(4) << std::setw(10)
            << combined_mean << std::setprecision(1) << std::setw(7) << stat_ev
            << std::setw(8) << cal_ev << std::setw(7) << bkg_ev << std::setw(9)
            << comb_ev << std::endl;
}
