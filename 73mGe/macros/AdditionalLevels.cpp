#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

struct LRTResult {
  TString name;
  Double_t min_nll;
  Int_t n_sample_float;
  Bool_t valid;
};

// Energy scale. Verified against Ba-133 (81-384), the Cd anchor (558.456,
// exact by construction) and Cd 805.89: residuals scatter +/-0.6 keV about
// zero with no trend, so the through-origin per-pixel gain is adequate here for
// a first pass and 0.6 keV is the systematic to quote on any centroid.
const Double_t E_SCALE_SYST_KEV = 0.6;

// A sample line with a tabulated energy is a reference, not an unknown, so it
// is pinned the same way the common Cd lines are. EGAF supplies tabulated
// energies for the prompt capture lines of every Ge isotope in the target.
struct Line {
  Double_t mu;
  TString tag;
  Bool_t fixed;
};

struct Region {
  TString name;
  Double_t lo, hi;
  std::vector<Line> common; // in BOTH runs: detector Cd, room background
  std::vector<Line> sample; // only with the Ge target: 73Ge, 77Ge, 77As
};

// Events from the interaction-summed tree, windowed on load. The whole tree is
// 375M entries for the 01/15 signal run; keeping only the fit window holds
// memory to tens of MB.
//
// Cached on (run, window), because competing models of the same window would
// otherwise each pay a full pass over the tree.
std::map<TString, std::vector<Double_t>> g_window_cache;

std::vector<Double_t> LoadSummed(const TString &run, Double_t lo, Double_t hi) {
  const TString key = TString::Format("%s_%.3f_%.3f", run.Data(), lo, hi);
  if (g_window_cache.count(key) > 0)
    return g_window_cache[key];

  std::vector<Double_t> out;
  TFile *f = IO::OpenForReading("filtered/" + run + ".root");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open filtered/" << run << ".root" << std::endl;
    return out;
  }
  TTree *t = static_cast<TTree *>(f->Get("total_energy_tree"));
  if (!t) {
    std::cerr << "ERROR: no total_energy_tree in " << run << std::endl;
    f->Close();
    return out;
  }
  Float_t e = 0;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("energykeV", 1);
  t->SetBranchAddress("energykeV", &e);
  Long64_t n = t->GetEntries();
  out.reserve(n / 50);
  for (Long64_t i = 0; i < n; i++) {
    t->GetEntry(i);
    if (e >= lo && e <= hi)
      out.push_back(e);
  }
  f->Close();
  g_window_cache[key] = out;
  return out;
}

LRTResult FitRegion(const TString &sig_run, const TString &bkg_run,
                    const Region &r, const TString &project_root) {
  LRTResult lrt;
  lrt.name = r.name;
  lrt.min_nll = 0;
  lrt.n_sample_float = 0;
  lrt.valid = kFALSE;
  std::cout << std::endl;
  std::cout << "================ " << r.name << "  (" << r.lo << "-" << r.hi
            << " keV) ================" << std::endl;
  std::vector<Double_t> sig = LoadSummed(sig_run, r.lo, r.hi);
  std::vector<Double_t> bkg = LoadSummed(bkg_run, r.lo, r.hi);
  if (sig.empty() || bkg.empty()) {
    std::cerr << "  no events" << std::endl;
    return lrt;
  }
  std::cout << "  sig " << sig.size() << " events   bkg " << bkg.size()
            << " events" << std::endl;

  Int_t n_common = (Int_t)r.common.size();
  Int_t n_sample = (Int_t)r.sample.size();
  std::vector<Double_t> bkg_mus, sig_mus;
  for (Int_t i = 0; i < n_common; i++) {
    bkg_mus.push_back(r.common[i].mu);
    sig_mus.push_back(r.common[i].mu);
  }
  // Sample lines are appended AFTER the common ones so LinkPeakShape can tie
  // signal peak i to background peak i for i < n_common.
  for (Int_t i = 0; i < n_sample; i++)
    sig_mus.push_back(r.sample[i].mu);

  RooFitUtils sim;
  sim.SetTailRatioMax(8.0);
  sim.SetRefitAfterLoad();

  // LINEAR background, not flat. The low-energy analysis can use a flat
  // background because its window sits on a locally flat region; here the
  // source Compton continuum falls visibly across an 18 keV window, and
  // forcing it flat gave chi2/ndf 22-452 and peaks wandering +/-3.5 keV.
  const Bool_t kFlat = kFALSE, kStep = kFALSE, kLowExp = kTRUE, kLowLin = kTRUE,
               kHighExp = kFALSE;
  // The common lines are 114Cd transitions with tabulated energies, so
  // their centroids are FIXED: they are reference points, not unknowns.
  // Letting four peaks float freely in 18 keV let them swap places.
  std::vector<Bool_t> bkg_fixed(n_common, kTRUE);
  std::vector<Bool_t> sig_fixed(n_common + n_sample, kFALSE);
  for (Int_t i = 0; i < n_common; i++)
    sig_fixed[i] = kTRUE;
  for (Int_t i = 0; i < n_sample; i++)
    sig_fixed[n_common + i] = r.sample[i].fixed;
  sim.AddChannel("bkg", bkg, r.lo, r.hi, Constants::BIN_WIDTH_KEV, n_common,
                 bkg_mus, kFlat, kStep, kLowExp, kLowLin, kHighExp, bkg_fixed,
                 kFALSE, kFALSE);
  sim.AddChannel("sig", sig, r.lo, r.hi, Constants::BIN_WIDTH_KEV,
                 n_common + n_sample, sig_mus, kFlat, kStep, kLowExp, kLowLin,
                 kHighExp, sig_fixed, kFALSE, kFALSE);
  for (Int_t i = 0; i < n_common; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  std::vector<FitResult> res = sim.FitSimultaneous(sig_run, "AddLev_" + r.name);
  const TString csv_base =
      project_root + "/plots/additional_levels/fits/csv_" + r.name;
  sim.DumpChannelCSV("bkg", csv_base + "_bkg");
  sim.DumpChannelCSV("sig", csv_base + "_sig");
  if (res.size() < 2) {
    std::cerr << "  SIMULTANEOUS FIT RETURNED NO RESULTS" << std::endl;
    return lrt;
  }
  if (!res[1].valid)
    std::cerr
        << "  WARNING: fit did not converge (results below may be unreliable)"
        << std::endl;
  const FitResult &S = res[1];
  std::cout << "  chi2/ndf  bkg " << res[0].reduced_chi2 << "   sig "
            << S.reduced_chi2 << std::endl;
  std::cout << "  " << std::left << std::setw(10) << "expected" << std::setw(12)
            << "fitted" << std::setw(10) << "+/-fit" << std::setw(10) << "resid"
            << std::setw(9) << "sigma" << "line" << std::endl;
  for (size_t i = 0; i < S.peaks.size() && i < sig_mus.size(); i++) {
    const PeakFitResult &p = S.peaks[i];
    TString tag = (i < (size_t)n_common)
                      ? ("[common] " + r.common[i].tag)
                      : ("[sample] " + r.sample[i - n_common].tag);
    std::cout << "  " << std::fixed << std::setprecision(3) << std::left
              << std::setw(10) << sig_mus[i] << std::setw(12) << p.mu
              << std::setw(10) << p.mu_error << std::showpos << std::setw(10)
              << (p.mu - sig_mus[i]) << std::noshowpos << std::setw(9)
              << p.sigma << tag << std::endl;
  }
  std::cout << "  (energy-scale systematic on every centroid: +/- "
            << E_SCALE_SYST_KEV << " keV)" << std::endl;

  lrt.min_nll = res[0].min_nll;
  // Count floating params contributed by sample peaks. Per peak:
  // sigma(1) + gaus_amplitude(1) + lowExp(2) + lowLin(2) = 6 always,
  // plus mu(1) if the peak is free (fixed=kFALSE).
  Int_t nf = 0;
  for (Int_t i = 0; i < n_sample; i++)
    nf += r.sample[i].fixed ? 6 : 7;
  lrt.n_sample_float = nf;
  lrt.valid = res[0].has_fit_diagnostics && std::isfinite(res[0].min_nll);
  std::cout << "  min_nll = " << std::setprecision(6) << lrt.min_nll
            << "   n_sample_float = " << nf
            << "   fit_status = " << res[0].fit_status
            << "   cov_qual = " << res[0].cov_qual << std::endl;
  return lrt;
}

void AdditionalLevels() {
  // IO::OpenForReading resolves against a base directory that defaults to a
  // RELATIVE "root_files", i.e. macros/root_files when run from here. Every
  // other macro in this project sets it explicitly; without this the loader
  // silently finds nothing and every region reports "no events".
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/additional_levels",
                                project_root + "/root_files");
  const TString sig = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;
  const TString bkg = Constants::NOSHIELDBACKGROUND_5PERCENT_20260115;

  // Common-line lists for each energy window.
  std::vector<Line> common_297 = {{300.87, "Cd-114", kTRUE},
                                  {304.86, "Cd-114", kTRUE}};
  std::vector<Line> common_351 = {{342.20, "Cd-110(n,g)", kTRUE},
                                  {345.07, "Cd-114", kTRUE},
                                  {359.20, "Cd-114", kTRUE},
                                  {361.50, "Cd-114", kTRUE}};
  std::vector<Line> common_709 = {{701.30, "Cd-110(n,g)", kTRUE},
                                  {707.42, "Cd-113(n,g)", kTRUE},
                                  {713.79, "Te-123(n,g)", kTRUE}};

  std::vector<Region> regions = {
      // --- 297 keV: VALIDATION (known 73Ge line) ---
      {"null_297", 290, 308, common_297, {}},
      {"validation_297",
       290,
       308,
       common_297,
       {{297.30, "73Ge 364.03->66.73 CONFIRMED", kFALSE}}},

      // --- 351 keV: adopted-only candidate ---
      {"null_351", 340, 362, common_351, {}},
      {"candidate_351",
       340,
       362,
       common_351,
       {{351.0, "73Ge 364.03->13.28 ADOPTED-ONLY", kFALSE}}},

      // --- 709 keV: 73Ge 708.8 vs 70Ge 708.15 ---
      {"null_709", 700, 718, common_709, {}},
      {"candidate_709_73Ge",
       700,
       718,
       common_709,
       {{708.8, "73Ge 776.66->68.75 ADOPTED-ONLY", kFALSE}}},
      {"candidate_709_70Ge",
       700,
       718,
       common_709,
       {{708.15, "70Ge(n,g) 708.15 TABULATED", kTRUE}}},
      {"candidate_709_both",
       700,
       718,
       common_709,
       {{708.15, "70Ge(n,g) 708.15 TABULATED", kTRUE},
        {708.8, "73Ge 776.66->68.75 ADOPTED-ONLY", kFALSE}}}};

  std::vector<LRTResult> results;
  for (size_t i = 0; i < regions.size(); i++)
    results.push_back(FitRegion(sig, bkg, regions[i], project_root));

  // Likelihood ratio tests. Each pair: (null, alternative).
  // Test statistic: lambda = 2*(NLL_null - NLL_alt) ~ chi2(delta_k)
  // delta_k estimated from sample-peak parameter count difference.
  struct LRTPair {
    TString null_name, alt_name, description;
  };
  std::vector<LRTPair> tests = {
      {"null_297", "validation_297", "73Ge 297.30 (validation)"},
      {"null_351", "candidate_351", "73Ge 351.0 (adopted-only)"},
      {"null_709", "candidate_709_73Ge", "73Ge 708.8 vs common-only"},
      {"null_709", "candidate_709_70Ge", "70Ge 708.15 vs common-only"},
      {"candidate_709_70Ge", "candidate_709_both",
       "adding 73Ge 708.8 to 70Ge model"}};

  std::cout << std::endl;
  std::cout << "======== LIKELIHOOD RATIO TESTS ========" << std::endl;
  std::cout << std::left << std::setw(45) << "Test" << std::setw(14)
            << "NLL_null" << std::setw(14) << "NLL_alt" << std::setw(12)
            << "lambda" << std::setw(8) << "delta_k" << std::setw(12)
            << "p-value" << "sigma" << std::endl;
  std::cout << std::string(115, '-') << std::endl;

  for (size_t t = 0; t < tests.size(); t++) {
    const LRTResult *rNull = nullptr;
    const LRTResult *rAlt = nullptr;
    for (size_t i = 0; i < results.size(); i++) {
      if (results[i].name == tests[t].null_name)
        rNull = &results[i];
      if (results[i].name == tests[t].alt_name)
        rAlt = &results[i];
    }
    if (!rNull || !rAlt || !rNull->valid || !rAlt->valid) {
      std::cout << std::setw(45) << tests[t].description << "  INVALID"
                << std::endl;
      continue;
    }
    Double_t lambda = 2.0 * (rNull->min_nll - rAlt->min_nll);
    Int_t delta_k = rAlt->n_sample_float - rNull->n_sample_float;
    if (delta_k <= 0)
      delta_k = 1;
    Double_t pval = (lambda > 0) ? TMath::Prob(lambda, delta_k) : 1.0;
    Double_t nsigma = (pval > 0 && pval < 1.0)
                          ? TMath::ErfcInverse(pval) * TMath::Sqrt(2.0)
                          : 0.0;
    std::cout << std::setw(45) << tests[t].description << std::setw(14)
              << std::fixed << std::setprecision(2) << rNull->min_nll
              << std::setw(14) << rAlt->min_nll << std::setw(12)
              << std::setprecision(2) << lambda << std::setw(8) << delta_k
              << std::setw(12) << std::scientific << std::setprecision(3)
              << pval << std::fixed << std::setprecision(1) << "  " << nsigma
              << "σ" << std::endl;
  }
}
