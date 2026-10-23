// Paper-style combiner for the Ge-73m gamma energy.
//
// This is the *simple* counterpart to CombineGeResult.cpp. It deliberately
// does NOT use BLUE and does NOT model the calibration correlations. It mimics
// the treatment in the published measurement:
//
//   1. one entry per dataset (runs measured by both methods are counted once,
//      keeping the rate-subtracted row);
//   2. error-weighted mean with sigma = sqrt(fit_err^2 + cal_err^2);
//   3. spread of the dataset means;
//   4. re-calculate the mean with the upper extreme removed, the lower extreme
//      removed, and both removed;
//   5. the largest |difference| from the all-data mean is the nuisance;
//   6. all terms are summed in quadrature into one quoted uncertainty.
//
// An optional second argument adds the calibration-peak-selection term
// (nominal Am+Pb line vs Pb-only line, evaluated at the Ge centroid). It is a
// property of the calibration model, so it is passed in rather than derived
// from these rows.
//
// bkg_err is printed but never enters sigma: it is a diagnostic (|cal(precal
// Ge) - postcal Ge|), not an independent uncertainty.
//
// Usage (ROOT auto-calls with the default):
//   root -b -q 'CombineGeResultAlt.cpp++'
//   root -b -q 'CombineGeResultAlt.cpp++("gaintransfer_hard")'
//   root -b -q 'CombineGeResultAlt.cpp++("tau20", 5.79)'
//   root -b -q 'CombineGeResultAlt.cpp++("tau20", 5.79, kFALSE)'
//       (the kFALSE keeps only the simultaneous in-situ rows; the rate-sub
//        method then enters solely as the measured M1-M2 systematic)

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
#include <sstream>
#include <vector>

struct Row {
  TString method;
  TString label;
  Double_t mu = 0, fit_err = 0, cal_err = 0, bkg_err = 0, chi2 = 0;
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
    rows.push_back(r);
  }
  return rows;
}

static Double_t Sigma(const Row &r) {
  return std::sqrt(r.fit_err * r.fit_err + r.cal_err * r.cal_err);
}

struct WMean {
  Int_t n = 0;
  Double_t mean = 0, err = 0, err_scaled = 0, chi2 = 0;
  Int_t ndf = 0, scale = 1;
};

// Weighted mean over the rows flagged in `use`.
WMean Weighted(const std::vector<Row> &rows, const std::vector<bool> &use) {
  WMean w;
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < rows.size(); i++) {
    if (!use[i])
      continue;
    Double_t s = Sigma(rows[i]);
    if (s <= 0)
      continue;
    Double_t wt = 1.0 / (s * s);
    sw += wt;
    swx += wt * rows[i].mu;
    w.n++;
  }
  if (sw <= 0)
    return w;
  w.mean = swx / sw;
  w.err = std::sqrt(1.0 / sw);
  for (size_t i = 0; i < rows.size(); i++) {
    if (!use[i])
      continue;
    Double_t s = Sigma(rows[i]);
    if (s <= 0)
      continue;
    Double_t d = rows[i].mu - w.mean;
    w.chi2 += d * d / (s * s);
  }
  w.ndf = w.n - 1;
  w.scale = (w.ndf > 0 && w.chi2 > w.ndf) ? (Int_t)1 : 0;
  Double_t sc = (w.ndf > 0) ? std::sqrt(w.chi2 / w.ndf) : 1.0;
  if (sc < 1.0)
    sc = 1.0;
  w.err_scaled = w.err * sc;
  return w;
}

void CombineGeResultAlt(const char *tag = "tau20", Double_t peak_sel_eV = 0.0,
                        Bool_t keep_ratesub = kTRUE) {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  TString path = project_root + "/results/ge_" + TString(tag) + ".result";
  std::vector<Row> all = ReadResult(path);
  if (all.empty()) {
    std::cerr << "no rows in " << path << std::endl;
    return;
  }

  // One entry per dataset: drop an in-situ row when the same run also has a
  // rate-subtracted row (the rate-sub row is the kept one). With
  // keep_ratesub = false the rate-subtracted rows are dropped instead, leaving
  // only the simultaneous in-situ fits -- i.e. the same datasets the BLUE
  // combiner quotes, with the rate-sub method then entering only as the
  // measured M1-M2 systematic.
  std::vector<Row> rows;
  std::vector<bool> use;
  for (size_t i = 0; i < all.size(); i++) {
    if (!keep_ratesub && all[i].method == "ratesub") {
      std::cout << "drop rate-sub row " << all[i].label
                << "  (in-situ only requested)" << std::endl;
      continue;
    }
    bool drop = false;
    if (keep_ratesub && all[i].method == "insitu") {
      for (size_t j = 0; j < all.size(); j++) {
        if (all[j].method != "ratesub")
          continue;
        if (all[j].label.BeginsWith(all[i].label)) {
          std::cout << "dedup: dropping in-situ " << all[i].label
                    << "  (keeping rate-sub " << all[j].label << ")"
                    << std::endl;
          drop = true;
          break;
        }
      }
    }
    if (drop)
      continue;
    rows.push_back(all[i]);
    use.push_back(true);
  }

  std::cout << std::fixed << std::setprecision(5);
  std::cout << std::endl
            << "===== " << tag << ": " << rows.size()
            << " datasets =====" << std::endl;
  std::cout << std::left << std::setw(46) << "Run" << std::right
            << std::setw(10) << "Eg[keV]" << std::setw(8) << "fit"
            << std::setw(8) << "cal" << std::setw(8) << "sigma" << std::setw(9)
            << "weight" << std::setw(7) << "chi2" << std::endl;
  for (size_t i = 0; i < rows.size(); i++) {
    Double_t s = Sigma(rows[i]);
    std::cout << std::left << std::setw(46) << rows[i].label << std::right
              << std::fixed << std::setprecision(5) << std::setw(10)
              << rows[i].mu << std::setprecision(1) << std::setw(8)
              << rows[i].fit_err * 1000 << std::setw(8)
              << rows[i].cal_err * 1000 << std::setw(8) << s * 1000
              << std::setw(9) << std::setprecision(2) << 1.0 / (s * s)
              << std::setprecision(2) << std::setw(7) << rows[i].chi2
              << std::endl;
  }

  WMean w = Weighted(rows, use);

  // Spread of the dataset means.
  Double_t lo = rows[0].mu, hi = rows[0].mu;
  size_t ilo = 0, ihi = 0;
  Double_t sum = 0;
  for (size_t i = 0; i < rows.size(); i++) {
    sum += rows[i].mu;
    if (rows[i].mu < lo) {
      lo = rows[i].mu;
      ilo = i;
    }
    if (rows[i].mu > hi) {
      hi = rows[i].mu;
      ihi = i;
    }
  }
  Double_t unweighted = sum / rows.size();
  Double_t s2 = 0;
  for (size_t i = 0; i < rows.size(); i++)
    s2 += (rows[i].mu - unweighted) * (rows[i].mu - unweighted);
  Double_t sd = std::sqrt(s2 / (rows.size() - 1));

  std::cout << std::setprecision(5);
  std::cout << std::endl
            << "spread (max - min)      = " << (hi - lo) * 1000 << " eV"
            << std::endl;
  std::cout << "unweighted mean         = " << unweighted
            << "   (sd = " << sd * 1000 << " eV)" << std::endl;
  std::cout << "error-weighted mean     = " << w.mean << " +/- " << w.err
            << " keV  (chi2/ndf = " << w.chi2 << "/" << w.ndf << ")"
            << std::endl;
  std::cout << "  scaled error (sqrt(chi2/ndf) when >1) = " << w.err_scaled
            << " keV" << std::endl;

  // Leave-one-out extremes.
  std::vector<bool> no_hi = use, no_lo = use, no_both = use;
  no_hi[ihi] = false;
  no_lo[ilo] = false;
  no_both[ihi] = false;
  no_both[ilo] = false;
  WMean wh = Weighted(rows, no_hi);
  WMean wl = Weighted(rows, no_lo);
  WMean wb = Weighted(rows, no_both);

  std::cout << std::endl
            << "----- outlier removal (as in the published approach) -----"
            << std::endl;
  std::cout << "  drop upper extreme (" << rows[ihi].label << ")" << std::endl;
  std::cout << "     = " << wh.mean
            << "   delta = " << (wh.mean - w.mean) * 1000 << " eV" << std::endl;
  std::cout << "  drop lower extreme (" << rows[ilo].label << ")" << std::endl;
  std::cout << "     = " << wl.mean
            << "   delta = " << (wl.mean - w.mean) * 1000 << " eV" << std::endl;
  std::cout << "  drop both" << std::endl;
  std::cout << "     = " << wb.mean
            << "   delta = " << (wb.mean - w.mean) * 1000 << " eV" << std::endl;
  Double_t d1 = std::fabs(wh.mean - w.mean) * 1000;
  Double_t d2 = std::fabs(wl.mean - w.mean) * 1000;
  Double_t d3 = std::fabs(wb.mean - w.mean) * 1000;
  Double_t nuisance = std::max(d1, std::max(d2, d3));

  // The published budget sums every term in quadrature into one number, so do
  // the same here: the weighted mean carries fit (+) calibration, the
  // outlier-removal nuisance is a separate term (not a re-weighting), and the
  // calibration-peak-selection term (nominal Am+Pb line vs Pb-only line) is
  // passed in by the caller because it is a property of the calibration model,
  // not of these rows.
  Double_t stat_ev = w.err * 1000.0;
  Double_t sel_ev = std::fabs(peak_sel_eV);
  Double_t total_ev =
      std::sqrt(stat_ev * stat_ev + nuisance * nuisance + sel_ev * sel_ev);

  std::cout << std::endl << "----- result -----" << std::endl;
  std::cout << std::setprecision(5);
  std::cout << "  weighted mean            = " << w.mean << " +/- " << w.err
            << " keV   (fit (+) cal, chi2/ndf = " << w.chi2 << "/" << w.ndf
            << ")" << std::endl;
  std::cout << "  outlier-removal nuisance = " << std::setprecision(5)
            << nuisance / 1000.0 << " keV   (max |delta|, "
            << std::setprecision(1) << nuisance << " eV)" << std::endl;
  if (sel_ev > 0)
    std::cout << "  calibration peaks        = " << std::setprecision(5)
              << sel_ev / 1000.0 << " keV   (Am+Pb line vs Pb-only, "
              << std::setprecision(1) << sel_ev << " eV)" << std::endl;
  std::cout << std::setprecision(5);
  std::cout << "  total (quadrature)       = " << total_ev / 1000.0 << " keV"
            << std::endl;
  std::cout << std::endl
            << "  E = " << w.mean << " +/- " << total_ev / 1000.0
            << " keV   (terms summed in quadrature, as in the published budget)"
            << std::endl;
  std::cout << std::setprecision(1);
  std::cout << "  E = " << std::setprecision(5) << w.mean << " +/- " << total_ev
            << " eV" << std::endl;
}
