// Aggregate the injection-test CSVs, print the bias/pull tables, classify the
// bias versus energy, and draw the bias-vs-energy figures in the standard
// PlottingUtils style.
//
// Usage (no arguments; edit SUM_RUN_TAG near the bottom):
//   root -q -l SummarizeInjection.cpp+
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCollection.h>
#include <TGraphErrors.h>
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
  if (r.method != "sim" || (r.status == 0 && r.cov_qual == 3))
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
  if (key == "CuShield10_0114")
    return "Cu 10% (01/14)";
  return key;
}

std::vector<std::string> PairOrder() {
  return {"CuShield10", "CdShield10", "CdShield25", "CuShield10_0114"};
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
const TString SUM_RUN_TAG = "pilot";
// ------------------------------------------------------------------------------

void SummarizeInjection() {
  const TString &run_tag = SUM_RUN_TAG;
  const TString project_root =
      TString(gSystem->Getenv("HOME")) + "/UM-ANSG/73mGe";
  const TString base = project_root + "/plots/injection_" + run_tag;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, base,
                                project_root + "/root_files");

  std::map<GroupKey, Acc> groups;
  std::map<AggKey, Acc> aggregates;

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
  md << std::endl << "## Interpretation" << std::endl << std::endl;
  const char *shapes[2] = {"model", "empirical"};
  const Int_t colors[4] = {kBlue, kRed, kGreen + 2, kMagenta + 1};
  std::vector<std::string> pair_order = PairOrder();

  // Distinct method names actually present, in first-seen order: "sim" plus
  // every rate-sub variant. The classification table below covers all of them;
  // only "sim" gets a figure (the rate-sub numbers stay in the tables).
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
      const Bool_t make_plot = (method == "sim");
      TString name =
          TString("bias_vs_energy_") + shape.c_str() + "_" + method.c_str();

      // Only the simultaneous-fit figure is drawn; the rate-sub numbers stay
      // in the tables. This keeps the plot count down to one figure per
      // lineshape.
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
          // Same style as PlottingUtils::AddLegend: boxed, white fill,
          // font 43.
          TLegend *leg = new TLegend(0.05, 0.92, 0.98, 0.72);
          leg->SetBorderSize(1);
          leg->SetFillColor(kWhite);
          leg->SetFillStyle(1001);
          leg->SetTextSize(30);
          leg->SetTextFont(43);
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
          legend_pad->cd();
          leg->Draw();
          PlottingUtils::SaveFigure(c, name, "", PlotSaveOptions::kLINEAR);
        }
      }

      // Classification per pair; rate-sub keeps its table but no figure.
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
