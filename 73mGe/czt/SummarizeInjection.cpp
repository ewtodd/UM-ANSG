// Aggregate the injection-test CSVs, print the bias/pull tables, classify the
// bias versus energy, and draw the bias-vs-energy figures in the standard
// PlottingUtils style.
//
// Usage (no arguments; edit SUM_RUN_TAG near the bottom):
//   root -q -l SummarizeInjection.cpp+
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TCollection.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Row {
  std::string pair, shape, method;
  Double_t e = 0, n_inj = 0, jitter = 0, mu = 0, mu_err = 0, residual = 0;
  Int_t trial = 0, status = 0, cov_qual = 0, valid = 0, rail = 0;
};

struct GroupKey {
  std::string pair, shape, method;
  Double_t e = 0;
  bool operator<(const GroupKey &o) const {
    if (pair != o.pair)
      return pair < o.pair;
    if (shape != o.shape)
      return shape < o.shape;
    if (method != o.method)
      return method < o.method;
    return e < o.e;
  }
};

struct AggKey {
  std::string pair, shape, method;
  bool operator<(const AggKey &o) const {
    if (pair != o.pair)
      return pair < o.pair;
    if (shape != o.shape)
      return shape < o.shape;
    return method < o.method;
  }
};

// One trial's in-situ fits with and without the escape lines, on the same
// data, for the escapes' effect trial by trial.
struct TrialKey {
  std::string pair, shape;
  Double_t e = 0;
  Int_t trial = 0;
  bool operator<(const TrialKey &o) const {
    if (pair != o.pair)
      return pair < o.pair;
    if (shape != o.shape)
      return shape < o.shape;
    if (e != o.e)
      return e < o.e;
    return trial < o.trial;
  }
};

struct PairedTrial {
  Bool_t has_with = kFALSE, has_without = kFALSE;
  Double_t with_escapes = 0, without_escapes = 0;
};

struct Acc {
  Int_t n = 0, nvalid = 0, nfail = 0, nrail = 0, nclean = 0;
  Double_t sres = 0, sres2 = 0, spull = 0, spull2 = 0;
  Double_t sj = 0, sj2 = 0, sjr = 0;
};

void Add(Acc &a, const Row &r) {
  a.n++;
  // The AU's valid flag accepts status-0 fits with a broken covariance; two
  // 70.5 keV trials returned mu errors of ~7 keV. Require a sane error so a
  // single degenerate fit cannot dominate a 20-trial mean.
  Bool_t sane =
      r.valid && std::isfinite(r.mu_err) && r.mu_err > 0 && r.mu_err < 0.5;
  if (!sane) {
    a.nfail++;
    return;
  }
  a.nvalid++;
  a.nrail += r.rail;
  // The simultaneous fit often ends at status 3 / covQual 2 in production;
  // count those separately from hard failures and from clean fits.
  Bool_t in_situ = r.method.compare(0, 3, "sim") == 0;
  if (!in_situ || (r.status == 0 && r.cov_qual == 3))
    a.nclean++;
  a.sres += r.residual;
  a.sres2 += r.residual * r.residual;
  if (r.mu_err > 0) {
    Double_t pull = r.residual / (r.mu_err * 1000.0);
    a.spull += pull;
    a.spull2 += pull * pull;
  }
  a.sj += r.jitter;
  a.sj2 += r.jitter * r.jitter;
  a.sjr += r.jitter * r.residual;
}

Double_t Bias(const Acc &a) { return a.nvalid ? a.sres / a.nvalid : 0.0; }
Double_t BiasSE(const Acc &a) {
  if (a.nvalid < 2)
    return 0.0;
  Double_t var = (a.sres2 - a.sres * a.sres / a.nvalid) / (a.nvalid - 1.0);
  return var > 0 ? std::sqrt(var / a.nvalid) : 0.0;
}
Double_t PullMean(const Acc &a) { return a.nvalid ? a.spull / a.nvalid : 0.0; }
Double_t PullRMS(const Acc &a) {
  return a.nvalid ? std::sqrt(a.spull2 / a.nvalid) : 0.0;
}
Double_t JitterCorr(const Acc &a) {
  if (a.nvalid < 3)
    return 0.0;
  Double_t num = a.nvalid * a.sjr - a.sj * a.sres;
  Double_t den1 = a.nvalid * a.sj2 - a.sj * a.sj;
  Double_t den2 = a.nvalid * a.sres2 - a.sres * a.sres;
  if (!(den1 > 0) || !(den2 > 0))
    return 0.0;
  return num / std::sqrt(den1 * den2);
}

std::vector<std::string> SplitCsv(const std::string &line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ','))
    out.push_back(field);
  return out;
}

Bool_t ParseRow(const std::string &line, Row &r) {
  std::vector<std::string> f = SplitCsv(line);
  if (f.size() < 15 || f[0] == "pair")
    return kFALSE;
  try {
    r.pair = f[0];
    r.shape = f[1];
    r.e = std::stod(f[2]);
    r.trial = std::stoi(f[3]);
    r.n_inj = std::stod(f[4]);
    r.jitter = std::stod(f[5]);
    r.method = f[6];
    r.mu = std::stod(f[7]);
    r.mu_err = std::stod(f[8]);
    r.residual = std::stod(f[9]);
    r.status = std::stoi(f[10]);
    r.cov_qual = std::stoi(f[11]);
    r.valid = std::stoi(f[13]);
    r.rail = std::stoi(f[14]);
  } catch (...) {
    return kFALSE;
  }
  return kTRUE;
}

struct LineFit {
  Int_t n = 0;
  Double_t mean = 0, mean_err = 0, chi2_const = 0, ndf_const = 0;
  Double_t slope = 0, slope_err = 0, intercept = 0, chi2_line = 0, ndf_line = 0;
};

LineFit FitBias(const std::vector<Double_t> &x, const std::vector<Double_t> &y,
                const std::vector<Double_t> &ey) {
  LineFit f;
  Double_t S0 = 0, S1 = 0, S2 = 0, T0 = 0, T1 = 0;
  for (size_t i = 0; i < y.size(); i++) {
    Double_t w = (ey[i] > 0) ? 1.0 / (ey[i] * ey[i]) : 1.0;
    S0 += w;
    S1 += w * x[i];
    S2 += w * x[i] * x[i];
    T0 += w * y[i];
    T1 += w * x[i] * y[i];
  }
  f.n = (Int_t)y.size();
  if (!(S0 > 0))
    return f;
  f.mean = T0 / S0;
  f.mean_err = std::sqrt(1.0 / S0);
  f.ndf_const = f.n - 1;
  for (size_t i = 0; i < y.size(); i++) {
    Double_t w = (ey[i] > 0) ? 1.0 / (ey[i] * ey[i]) : 1.0;
    Double_t d = y[i] - f.mean;
    f.chi2_const += w * d * d;
  }
  Double_t D = S0 * S2 - S1 * S1;
  if (D > 0) {
    f.slope = (S0 * T1 - S1 * T0) / D;
    f.intercept = (S2 * T0 - S1 * T1) / D;
    f.slope_err = std::sqrt(S0 / D);
    f.ndf_line = f.n - 2;
    for (size_t i = 0; i < y.size(); i++) {
      Double_t w = (ey[i] > 0) ? 1.0 / (ey[i] * ey[i]) : 1.0;
      Double_t d = y[i] - (f.intercept + f.slope * x[i]);
      f.chi2_line += w * d * d;
    }
  }
  return f;
}

std::string Classify(const std::vector<Double_t> &y,
                     const std::vector<Double_t> &ey, const LineFit &fit) {
  if (y.empty())
    return "no usable trials";
  if (fit.n < 2)
    return "insufficient energy points";
  Bool_t all_zero = kTRUE;
  for (size_t i = 0; i < y.size(); i++) {
    if (!(ey[i] > 0) || std::fabs(y[i]) > 2.0 * ey[i])
      all_zero = kFALSE;
  }
  if (all_zero)
    return "zero-consistent";
  // A linear slope alone misses a V-shaped or otherwise non-monotonic energy
  // dependence: the Cu 01/13 scan flips sign at both ends while its line fit
  // has a slope consistent with zero. The constant-fit chi2/ndf catches it.
  if (fit.ndf_const > 0 && fit.chi2_const / fit.ndf_const > 2.0)
    return "ENERGY-DEPENDENT (non-linear): stop and report; do not correct";
  if (fit.slope_err > 0 && std::fabs(fit.slope) > 2.0 * fit.slope_err)
    return "ENERGY-DEPENDENT: stop and report; do not correct";
  return "energy-independent shift: correct and propagate";
}

std::string PairLabel(const std::string &key) {
  if (key == "CuShield10")
    return "Cu 10% (01/13)";
  if (key == "CdShield10")
    return "Cd 10% (01/13)";
  if (key == "CdShield25")
    return "Cd 25% (01/13)";
  return key;
}

std::vector<std::string> PairOrder() {
  return {"CuShield10", "CdShield10", "CdShield25"};
}

// The in-situ fit's bias per pair with the escape lines and without, on the
// same trials: what leaving out the escapes the data carry costs.
void DrawEscapeBiasByPair(const std::map<AggKey, Acc> &aggregates) {
  std::vector<std::string> pairs = PairOrder();
  TGraphErrors *with_escapes = new TGraphErrors();
  TGraphErrors *without_escapes = new TGraphErrors();
  Double_t lo = 0, hi = 0;
  for (size_t pi = 0; pi < pairs.size(); pi++) {
    const char *methods[2] = {"sim", "sim_noesc"};
    TGraphErrors *graphs[2] = {with_escapes, without_escapes};
    for (Int_t m = 0; m < 2; m++) {
      std::map<AggKey, Acc>::const_iterator it =
          aggregates.find(AggKey{pairs[pi], "model", methods[m]});
      if (it == aggregates.end() || it->second.nvalid == 0)
        continue;
      Double_t x = pi + 1 + (m == 0 ? -0.1 : 0.1);
      Int_t i = graphs[m]->GetN();
      graphs[m]->SetPoint(i, x, Bias(it->second));
      graphs[m]->SetPointError(i, 0, BiasSE(it->second));
      lo = std::min(lo, Bias(it->second) - BiasSE(it->second));
      hi = std::max(hi, Bias(it->second) + BiasSE(it->second));
    }
  }
  if (with_escapes->GetN() == 0 || without_escapes->GetN() == 0)
    return;
  Double_t pad = 0.25 * (hi - lo);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas();
  TH1F *frame =
      new TH1F("escape_bias_frame", ";;Bias (#mu_{fit} - E_{inj}) [eV]",
               (Int_t)pairs.size(), 0.5, pairs.size() + 0.5);
  for (size_t pi = 0; pi < pairs.size(); pi++)
    frame->GetXaxis()->SetBinLabel(pi + 1, PairLabel(pairs[pi]).c_str());
  frame->SetMinimum(lo - pad);
  frame->SetMaximum(hi + 2.0 * pad);
  frame->SetStats(0);
  frame->Draw("AXIS");
  TLine *zero = new TLine(0.5, 0.0, pairs.size() + 0.5, 0.0);
  zero->SetLineStyle(2);
  zero->SetLineColor(kGray + 2);
  zero->Draw("SAME");
  PlottingUtils::ConfigureGraph(with_escapes, kBlue + 1);
  PlottingUtils::ConfigureGraph(without_escapes, kRed + 1);
  with_escapes->SetMarkerStyle(20);
  without_escapes->SetMarkerStyle(21);
  with_escapes->Draw("P SAME");
  without_escapes->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.75, 0.74, 0.9);
  leg->AddEntry(with_escapes, "Fit with the escape lines (production)", "pe");
  leg->AddEntry(without_escapes, "Same trials, fit without them", "pe");
  leg->Draw();
  PlottingUtils::SaveFigure(c, "escape_bias_by_pair", "",
                            PlotSaveOptions::kLINEAR);
}

// A pair's trials, each one's residual without the escape lines against the
// same events' residual with them: the trial noise moves both together, the
// escapes move every trial off the diagonal by the same amount.
void DrawEscapePaired(const std::string &pair,
                      const std::vector<Double_t> &with_escapes,
                      const std::vector<Double_t> &without_escapes) {
  const Int_t n = (Int_t)with_escapes.size();
  if (n < 5)
    return;
  Double_t lo = 1e30, hi = -1e30, sum = 0, sum2 = 0;
  for (Int_t i = 0; i < n; i++) {
    lo = std::min(lo, std::min(with_escapes[i], without_escapes[i]));
    hi = std::max(hi, std::max(with_escapes[i], without_escapes[i]));
    Double_t d = without_escapes[i] - with_escapes[i];
    sum += d;
    sum2 += d * d;
  }
  Double_t shift = sum / n;
  Double_t rms = std::sqrt(std::max(sum2 / n - shift * shift, 0.0));
  Double_t pad = 0.1 * (hi - lo);
  lo -= pad;
  hi += pad;
  TGraph *trials = new TGraph(n, with_escapes.data(), without_escapes.data());
  TCanvas *c = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      trials, kRed + 1,
      ";#mu_{fit} - E_{inj}, fit with the escape lines [eV];"
      "Same trial, fit without them [eV]");
  trials->SetMarkerStyle(20);
  trials->SetMarkerColor(kRed + 1);
  trials->GetXaxis()->SetLimits(lo, hi);
  trials->SetMinimum(lo);
  trials->SetMaximum(hi);
  trials->Draw("AP");
  TLine *diagonal = new TLine(lo, lo, hi, hi);
  diagonal->SetLineStyle(2);
  diagonal->SetLineColor(kGray + 2);
  diagonal->Draw("SAME");
  TLine *shifted = new TLine(lo, lo + shift, hi - shift, hi);
  shifted->SetLineColor(kRed + 1);
  shifted->Draw("SAME");
  trials->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.42, 0.9, 0.2, 0.42);
  leg->AddEntry(trials, Form("%s, %d trials", PairLabel(pair).c_str(), n), "p");
  leg->AddEntry(diagonal, "Escapes without effect", "l");
  leg->AddEntry(shifted,
                Form("Shift %.1f #pm %.1f eV (rms %.1f eV)", shift,
                     rms / std::sqrt(n - 1.0), rms),
                "l");
  leg->Draw();
  PlottingUtils::SaveFigure(c, TString("escape_paired_") + pair.c_str(), "",
                            PlotSaveOptions::kLINEAR);
}

void PrintRow(std::ostream &out, const std::string &pair,
              const std::string &shape, const std::string &method, Double_t e,
              const Acc &a) {
  out << std::left << std::setw(18) << PairLabel(pair) << std::setw(10) << shape
      << std::setw(9) << method << std::fixed << std::setprecision(3)
      << std::setw(10) << e << std::setw(6) << a.nvalid << std::setw(10)
      << std::setprecision(2) << Bias(a) << " +/- " << std::setw(7) << BiasSE(a)
      << std::setw(9) << std::setprecision(3) << PullMean(a) << std::setw(9)
      << PullRMS(a) << std::setw(8) << std::setprecision(1)
      << (a.n ? 100.0 * a.nfail / a.n : 0.0) << std::setw(8)
      << (a.nvalid ? 100.0 * a.nclean / a.nvalid : 0.0) << std::setw(8)
      << (a.n ? 100.0 * a.nrail / a.n : 0.0) << std::setw(9)
      << std::setprecision(3) << JitterCorr(a) << std::endl;
}

} // namespace

// ---------------------------------------------------------------- run controls
// Edit and re-run. No command-line arguments by design. Must match the
// INJ_RUN_TAG that produced the CSVs.
const TString SUM_RUN_TAG = "closure";
// ------------------------------------------------------------------------------

void SummarizeInjection() {
  const TString &run_tag = SUM_RUN_TAG;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  const TString base = project_root + "/plots/testbed/injection_" + run_tag;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, base,
                                project_root + "/root_files");

  std::map<GroupKey, Acc> groups;
  std::map<AggKey, Acc> aggregates;
  std::map<TrialKey, PairedTrial> paired;

  TSystemDirectory dir(base + "/results", base + "/results");
  TList *files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "ERROR: no results directory " << base + "/results"
              << std::endl;
    return;
  }
  Int_t n_rows = 0;
  TIter next(files);
  while (TSystemFile *file = static_cast<TSystemFile *>(next())) {
    TString name = file->GetName();
    if (!name.EndsWith(".csv"))
      continue;
    std::ifstream in((base + "/results/" + name).Data());
    std::string line;
    while (std::getline(in, line)) {
      Row r;
      if (!ParseRow(line, r))
        continue;
      n_rows++;
      GroupKey gk{r.pair, r.shape, r.method, r.e};
      Add(groups[gk], r);
      AggKey ak{r.pair, r.shape, r.method};
      Add(aggregates[ak], r);
      Bool_t sane =
          r.valid && std::isfinite(r.mu_err) && r.mu_err > 0 && r.mu_err < 0.5;
      if (sane && (r.method == "sim" || r.method == "sim_noesc")) {
        PairedTrial &p = paired[TrialKey{r.pair, r.shape, r.e, r.trial}];
        if (r.method == "sim") {
          p.has_with = kTRUE;
          p.with_escapes = r.residual;
        } else {
          p.has_without = kTRUE;
          p.without_escapes = r.residual;
        }
      }
    }
  }
  if (n_rows == 0) {
    std::cerr << "ERROR: no rows found under " << base + "/results"
              << std::endl;
    return;
  }

  std::ofstream md((base + "/injection_summary.md").Data());
  md << "# Injection (closure) test — " << run_tag << std::endl << std::endl;
  md << "Bias is `mu_fit - E_inj` in eV; SE is the standard error of the mean "
        "over valid trials. Pull is `(mu_fit - E_inj) / sigma_fit`; its RMS "
        "should be ~1 if the fit errors are correct."
     << std::endl
     << std::endl;

  std::cout << "=== Aggregate per pair / lineshape / method ===" << std::endl;
  std::cout << std::left << std::setw(18) << "Pair" << std::setw(10) << "Shape"
            << std::setw(9) << "Method" << std::setw(6) << "N" << std::setw(18)
            << "Bias [eV]" << std::setw(9) << "Pull" << std::setw(9) << "RMS"
            << std::setw(8) << "Fail%" << std::setw(8) << "Clean%"
            << std::setw(8) << "Rail%" << std::setw(9) << "corr" << std::endl;
  for (std::map<AggKey, Acc>::const_iterator it = aggregates.begin();
       it != aggregates.end(); ++it) {
    const AggKey &k = it->first;
    const Acc &a = it->second;
    std::cout << std::left << std::setw(18) << PairLabel(k.pair)
              << std::setw(10) << k.shape << std::setw(9) << k.method
              << std::setw(6) << a.nvalid << std::fixed << std::setprecision(2)
              << std::setw(10) << Bias(a) << " +/- " << std::setw(7)
              << BiasSE(a) << std::setw(9) << std::setprecision(3)
              << PullMean(a) << std::setw(9) << PullRMS(a) << std::setw(8)
              << std::setprecision(1) << (a.n ? 100.0 * a.nfail / a.n : 0.0)
              << std::setw(8) << (a.nvalid ? 100.0 * a.nclean / a.nvalid : 0.0)
              << std::setw(8) << (a.n ? 100.0 * a.nrail / a.n : 0.0)
              << std::setw(9) << std::setprecision(3) << JitterCorr(a)
              << std::endl;
  }

  std::cout << std::endl << "=== Per energy ===" << std::endl;
  std::cout << std::left << std::setw(18) << "Pair" << std::setw(10) << "Shape"
            << std::setw(9) << "Method" << std::setw(10) << "E_inj"
            << std::setw(6) << "N" << std::setw(18) << "Bias [eV]"
            << std::setw(9) << "Pull" << std::setw(9) << "RMS" << std::setw(8)
            << "Fail%" << std::setw(8) << "Clean%" << std::setw(8) << "Rail%"
            << std::setw(9) << "corr" << std::endl;
  for (std::map<GroupKey, Acc>::const_iterator it = groups.begin();
       it != groups.end(); ++it) {
    const GroupKey &k = it->first;
    PrintRow(std::cout, k.pair, k.shape, k.method, k.e, it->second);
  }

  // Markdown: aggregate table.
  md << "## Aggregate" << std::endl << std::endl;
  md << "| Pair | Shape | Method | N | Bias [eV] | SE [eV] | Pull mean | Pull "
        "RMS | Fail % | Clean % | Rail % | corr(jitter, residual) |"
     << std::endl;
  md << "|---|---|---|---|---|---|---|---|---|---|---|---|" << std::endl;
  for (std::map<AggKey, Acc>::const_iterator it = aggregates.begin();
       it != aggregates.end(); ++it) {
    const AggKey &k = it->first;
    const Acc &a = it->second;
    md << "| " << PairLabel(k.pair) << " | " << k.shape << " | " << k.method
       << " | " << a.nvalid << " | " << std::fixed << std::setprecision(2)
       << Bias(a) << " | " << BiasSE(a) << " | " << std::setprecision(3)
       << PullMean(a) << " | " << PullRMS(a) << " | " << std::setprecision(1)
       << (a.n ? 100.0 * a.nfail / a.n : 0.0) << " | "
       << (a.nvalid ? 100.0 * a.nclean / a.nvalid : 0.0) << " | "
       << (a.n ? 100.0 * a.nrail / a.n : 0.0) << " | " << std::setprecision(3)
       << JitterCorr(a) << " |" << std::endl;
  }

  // Markdown: per-energy table.
  md << std::endl << "## Per injection energy" << std::endl << std::endl;
  md << "| Pair | Shape | Method | E_inj [keV] | N | Bias [eV] | SE [eV] | "
        "Pull mean | Pull RMS | Fail % | Clean % | Rail % | corr |"
     << std::endl;
  md << "|---|---|---|---|---|---|---|---|---|---|---|---|---|" << std::endl;
  for (std::map<GroupKey, Acc>::const_iterator it = groups.begin();
       it != groups.end(); ++it) {
    const GroupKey &k = it->first;
    const Acc &a = it->second;
    md << "| " << PairLabel(k.pair) << " | " << k.shape << " | " << k.method
       << " | " << std::fixed << std::setprecision(3) << k.e << " | "
       << a.nvalid << " | " << std::setprecision(2) << Bias(a) << " | "
       << BiasSE(a) << " | " << std::setprecision(3) << PullMean(a) << " | "
       << PullRMS(a) << " | " << std::setprecision(1)
       << (a.n ? 100.0 * a.nfail / a.n : 0.0) << " | "
       << (a.nvalid ? 100.0 * a.nclean / a.nvalid : 0.0) << " | "
       << (a.n ? 100.0 * a.nrail / a.n : 0.0) << " | " << std::setprecision(3)
       << JitterCorr(a) << " |" << std::endl;
  }

  // Classification and bias-vs-energy figures.
  // The escapes' effect: the fit without them minus the fit with them, trial
  // by trial on the same data, per cell.
  std::map<GroupKey, std::vector<Double_t>> escape_effect;
  for (std::map<TrialKey, PairedTrial>::const_iterator it = paired.begin();
       it != paired.end(); ++it)
    if (it->second.has_with && it->second.has_without)
      escape_effect[GroupKey{it->first.pair, it->first.shape, "escapes",
                             it->first.e}]
          .push_back(it->second.without_escapes - it->second.with_escapes);
  if (!escape_effect.empty()) {
    std::cout << "=== Escape effect: fit without them - fit with them, paired "
                 "by trial [eV] ==="
              << std::endl;
    md << std::endl
       << "## Escape effect (fit without the escape lines minus with, paired "
          "by trial)"
       << std::endl
       << std::endl;
    md << "| Pair | Shape | E_inj [keV] | N | Mean [eV] | SE [eV] | RMS [eV] |"
       << std::endl;
    md << "|---|---|---|---|---|---|---|" << std::endl;
    for (std::map<GroupKey, std::vector<Double_t>>::const_iterator it =
             escape_effect.begin();
         it != escape_effect.end(); ++it) {
      const std::vector<Double_t> &d = it->second;
      Double_t sum = 0, sum2 = 0;
      for (size_t i = 0; i < d.size(); i++) {
        sum += d[i];
        sum2 += d[i] * d[i];
      }
      Double_t n = d.size();
      Double_t mean = sum / n;
      Double_t rms = std::sqrt(std::max(sum2 / n - mean * mean, 0.0));
      Double_t se = (n > 1) ? rms / std::sqrt(n - 1) : 0.0;
      std::cout << std::left << std::setw(18) << PairLabel(it->first.pair)
                << std::setw(10) << it->first.shape << std::fixed
                << std::setprecision(3) << std::setw(10) << it->first.e
                << std::setw(6) << (Int_t)n << std::setprecision(2)
                << std::setw(10) << mean << " +/- " << std::setw(7) << se
                << "  rms " << rms << std::endl;
      md << "| " << PairLabel(it->first.pair) << " | " << it->first.shape
         << " | " << std::setprecision(3) << it->first.e << " | " << (Int_t)n
         << " | " << std::setprecision(2) << mean << " | " << se << " | " << rms
         << " |" << std::endl;
    }
    DrawEscapeBiasByPair(aggregates);
    md << std::endl
       << "![escape bias by pair](escape_bias_by_pair.png)" << std::endl;
    std::vector<std::string> pairs = PairOrder();
    for (size_t pi = 0; pi < pairs.size(); pi++) {
      std::vector<Double_t> with_escapes, without_escapes;
      for (std::map<TrialKey, PairedTrial>::const_iterator it = paired.begin();
           it != paired.end(); ++it)
        if (it->first.pair == pairs[pi] && it->second.has_with &&
            it->second.has_without) {
          with_escapes.push_back(it->second.with_escapes);
          without_escapes.push_back(it->second.without_escapes);
        }
      DrawEscapePaired(pairs[pi], with_escapes, without_escapes);
    }
  }

  md << std::endl << "## Interpretation" << std::endl << std::endl;
  const char *shapes[2] = {"model", "empirical"};
  const Int_t colors[4] = {kBlue, kRed, kGreen + 2, kMagenta + 1};
  std::vector<std::string> pair_order = PairOrder();

  // Distinct method names actually present, in first-seen order. Each gets
  // the classification table and a bias-vs-energy figure.
  std::vector<std::string> method_order;
  for (std::map<GroupKey, Acc>::const_iterator m_it = groups.begin();
       m_it != groups.end(); ++m_it) {
    bool seen = false;
    for (size_t m = 0; m < method_order.size(); m++)
      if (method_order[m] == m_it->first.method)
        seen = true;
    if (!seen)
      method_order.push_back(m_it->first.method);
  }

  for (Int_t si = 0; si < 2; si++) {
    for (size_t mi = 0; mi < method_order.size(); mi++) {
      std::string shape = shapes[si], method = method_order[mi];
      const Bool_t make_plot =
          (method == "sim" || method == "sim_noesc" || method == "ratesub");
      TString name =
          TString("bias_vs_energy_") + shape.c_str() + "_" + method.c_str();

      if (make_plot) {
        std::vector<TString> legend_labels;
        std::vector<TGraphErrors *> graphs;
        std::vector<Double_t> y_lo, y_hi;

        for (size_t pi = 0; pi < pair_order.size(); pi++) {
          std::vector<Double_t> xs, ys, eys, exs;
          for (std::map<GroupKey, Acc>::const_iterator it = groups.begin();
               it != groups.end(); ++it) {
            const GroupKey &k = it->first;
            if (k.pair != pair_order[pi] || k.shape != shape ||
                k.method != method)
              continue;
            if (it->second.nvalid == 0)
              continue;
            xs.push_back(k.e);
            ys.push_back(Bias(it->second));
            eys.push_back(BiasSE(it->second));
            exs.push_back(0.0);
            y_lo.push_back(Bias(it->second) - BiasSE(it->second));
            y_hi.push_back(Bias(it->second) + BiasSE(it->second));
          }
          if (xs.empty())
            continue;
          Acc agg;
          AggKey ak{pair_order[pi], shape, method};
          if (aggregates.count(ak))
            agg = aggregates[ak];
          legend_labels.push_back(
              TString::Format("%s  (mean %.1f eV)",
                              PairLabel(pair_order[pi]).c_str(), Bias(agg)));
          TGraphErrors *g = new TGraphErrors((Int_t)xs.size(), xs.data(),
                                             ys.data(), exs.data(), eys.data());
          PlottingUtils::ConfigureGraph(
              g, colors[pi],
              ";Injection energy [keV];Bias (#mu_{fit} - E_{inj}) [eV]");
          graphs.push_back(g);
        }

        if (!graphs.empty()) {
          Double_t lo = 1e30, hi = -1e30;
          for (size_t i = 0; i < y_lo.size(); i++) {
            if (y_lo[i] < lo)
              lo = y_lo[i];
            if (y_hi[i] > hi)
              hi = y_hi[i];
          }
          lo = std::min(lo, 0.0);
          hi = std::max(hi, 0.0);
          if (hi - lo < 5.0) {
            Double_t mid = 0.5 * (hi + lo);
            lo = mid - 2.5;
            hi = mid + 2.5;
          }
          Double_t pad = 0.15 * (hi - lo);
          lo -= pad;
          hi += pad;

          Double_t legend_w =
              PlottingUtils::LegendWidthPx(legend_labels, "Injection bias");
          TPad *plot = nullptr, *legend_pad = nullptr;
          TCanvas *c = PlottingUtils::GetConfiguredCanvasWithSideLegend(
              plot, legend_pad, (Int_t)legend_w);
          plot->cd();
          TLegend *leg = PlottingUtils::AddSideLegend(
              plot, legend_pad, (Int_t)legend_labels.size());
          for (size_t i = 0; i < graphs.size(); i++) {
            TGraphErrors *g = graphs[i];
            g->SetMinimum(lo);
            g->SetMaximum(hi);
            g->GetXaxis()->SetLimits(66.0, 71.2);
            g->Draw(i == 0 ? "AP" : "P SAME");
            leg->AddEntry(g, legend_labels[i], "lp");
          }
          TLine *zero = new TLine(66.0, 0.0, 71.2, 0.0);
          zero->SetLineStyle(2);
          zero->SetLineColor(kGray + 2);
          zero->Draw("SAME");
          PlottingUtils::SaveFigure(c, name, "", PlotSaveOptions::kLINEAR);
        }
      }

      // Classification per pair.
      md << std::endl
         << "### " << shape << " / " << method << std::endl
         << std::endl;
      if (make_plot)
        md << "![bias vs energy](" << name << ".png)" << std::endl << std::endl;
      md << "| Pair | Bias fit | Slope [eV/keV] | chi2/ndf (const) | "
            "chi2/ndf (line) | Classification |"
         << std::endl;
      md << "|---|---|---|---|---|---|" << std::endl;
      for (size_t pi = 0; pi < pair_order.size(); pi++) {
        std::vector<Double_t> xs, ys, eys;
        for (std::map<GroupKey, Acc>::const_iterator it = groups.begin();
             it != groups.end(); ++it) {
          const GroupKey &k = it->first;
          if (k.pair != pair_order[pi] || k.shape != shape ||
              k.method != method || it->second.nvalid == 0)
            continue;
          xs.push_back(k.e);
          ys.push_back(Bias(it->second));
          eys.push_back(BiasSE(it->second));
        }
        if (ys.empty())
          continue;
        LineFit fit = FitBias(xs, ys, eys);
        std::string cls = Classify(ys, eys, fit);
        md << "| " << PairLabel(pair_order[pi]) << " | " << std::fixed
           << std::setprecision(2) << fit.mean << " +/- " << fit.mean_err
           << " eV | " << std::setprecision(4) << fit.slope << " +/- "
           << fit.slope_err << " | " << std::setprecision(3)
           << (fit.ndf_const > 0 ? fit.chi2_const / fit.ndf_const : 0.0)
           << " | " << (fit.ndf_line > 0 ? fit.chi2_line / fit.ndf_line : 0.0)
           << " | " << cls << " |" << std::endl;
      }
    }
  }

  md << std::endl << "## Pre-fixed interpretation" << std::endl << std::endl;
  md << "- Zero-consistent bias: unbiased; quote the bias uncertainty as the "
        "systematic."
     << std::endl;
  md << "- Significant energy-independent bias: correct and propagate."
     << std::endl;
  md << "- Energy-dependent bias (line slope > 2 sigma): stop and report; do "
        "not correct."
     << std::endl;
  md.close();
  std::cout << std::endl
            << "Wrote " << base + "/injection_summary.md" << std::endl;
}
