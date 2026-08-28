// Search for 73Ge gamma transitions above 150 keV in interaction-summed CZT
// data (section "Additional 73Ge Levels?" of the PRC draft).
//
// WHY A SIMULTANEOUS FIT RATHER THAN BACKGROUND SUBTRACTION.
//
// The obvious approach -- subtract a live-time-normalised background run and
// fit the residual -- does not work here, for two reasons found the hard way:
//
//  1. The Ge sample adds ~2.85x continuum across the WHOLE spectrum (its own
//     Compton), so sig-minus-bkg is positive everywhere and "net > 0" is not
//     evidence of a line. Scaling the background to match that continuum
//     instead (k = 5.15, the counts ratio) over-subtracts every line whose RATE
//     is common to both runs: S - kB = R*t_s(1 - 2.85), a deep negative trench.
//     One scale factor cannot both cancel common lines and match the continuum.
//  2. Even with live-time scaling, the region is dense with lines and a
//     Gaussian-on-a-line fit cannot separate blends: 430.3/432.7 (2.4 keV
//     apart) both converged on the same peak, and 561.6 was swallowed by the
//     Cd anchor 3.1 keV away. CZT peaks also carry low-energy tails, so
//     gaus+pol1 returns chi2/ndf 2-4 even where it "works".
//
// Instead both runs are fitted SIMULTANEOUSLY with one shared lineshape, the
// same structure the low-energy 68.75 analysis uses. Lines present in BOTH runs
// (113Cd(n,g) capture from the detector itself, room background) go in the
// background channel and are shape-linked into the signal channel, where they
// are constrained by the background data rather than free to absorb signal.
// Lines that appear only with the Ge sample -- 77Ge/77As activation and any
// 73Ge transition -- are extra peaks in the signal channel alone.
//
// THREE CONTAMINANT LAYERS, all established from data:
//   detector    113Cd(n,g), notably 558.456 (the pixel-calibration anchor) and
//               805.89; present in every run, scales with neutron flux
//   sample      77Ge (11.2 h) and 77As (38.8 h) from 76Ge, which ACCUMULATE
//               across a multi-day campaign; 75Ge from 74Ge
//   room        the usual, though 214Bi looks weak here -- its 609.3 line is
//               absent even where 351.9 would demand it
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <vector>

// Energy scale. Verified against Ba-133 (81-384), the Cd anchor (558.456,
// exact by construction) and Cd 805.89: residuals scatter +/-0.6 keV about
// zero with no trend, so the through-origin per-pixel gain is adequate here and
// 0.6 keV is the systematic to quote on any centroid.
const Double_t E_SCALE_SYST_KEV = 0.6;

struct Line {
  Double_t mu;
  TString tag;
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
std::vector<Double_t> LoadSummed(const TString &run, Double_t lo, Double_t hi) {
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
  t->SetBranchAddress("energykeV", &e);
  Long64_t n = t->GetEntries();
  out.reserve(n / 50);
  for (Long64_t i = 0; i < n; i++) {
    t->GetEntry(i);
    if (e >= lo && e <= hi)
      out.push_back(e);
  }
  f->Close();
  return out;
}

void FitRegion(const TString &sig_run, const TString &bkg_run,
               const Region &r) {
  std::cout << "\n================ " << r.name << "  (" << r.lo << "-" << r.hi
            << " keV) ================" << std::endl;
  std::vector<Double_t> sig = LoadSummed(sig_run, r.lo, r.hi);
  std::vector<Double_t> bkg = LoadSummed(bkg_run, r.lo, r.hi);
  if (sig.empty() || bkg.empty()) {
    std::cerr << "  no events" << std::endl;
    return;
  }
  std::cout << "  sig " << sig.size() << " events   bkg " << bkg.size()
            << " events" << std::endl;

  Int_t n_common = (Int_t)r.common.size();
  Int_t n_sample = (Int_t)r.sample.size();
  std::vector<Double_t> bkg_mus, sig_mus;
  for (const Line &l : r.common) {
    bkg_mus.push_back(l.mu);
    sig_mus.push_back(l.mu);
  }
  // Sample lines are appended AFTER the common ones so LinkPeakShape can tie
  // signal peak i to background peak i for i < n_common.
  for (const Line &l : r.sample)
    sig_mus.push_back(l.mu);

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
  sim.AddChannel("bkg", bkg, r.lo, r.hi, Constants::BIN_WIDTH_KEV, n_common,
                 bkg_mus, kFlat, kStep, kLowExp, kLowLin, kHighExp, bkg_fixed,
                 kFALSE, kFALSE);
  sim.AddChannel("sig", sig, r.lo, r.hi, Constants::BIN_WIDTH_KEV,
                 n_common + n_sample, sig_mus, kFlat, kStep, kLowExp, kLowLin,
                 kHighExp, sig_fixed, kFALSE, kFALSE);
  for (Int_t i = 0; i < n_common; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  std::vector<FitResult> res = sim.FitSimultaneous(sig_run, "AddLev_" + r.name);
  if (res.size() < 2 || !res[1].valid) {
    std::cerr << "  SIMULTANEOUS FIT FAILED" << std::endl;
    return;
  }
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

  std::vector<Region> regions = {
      // VALIDATION. 297.30 is confirmed in both the adopted and (n,g) columns,
      // so this region measures whether the method recovers a known 73Ge line.
      {"validation_297",
       290,
       308,
       {{300.87, "Cd-114"}, {304.86, "Cd-114"}},
       {{297.30, "73Ge 364.03->66.73 CONFIRMED"}}},
      // The 351.0 candidate. Adopted-only, and its own level scheme predicts
      // 350.98 from the confirmed 297.30 sibling. A single free peak here lands
      // at ~352.5 in the residual analysis, 1.5 keV away.
      {"candidate_351",
       344,
       362,
       {{345.07, "Cd-114"}, {359.20, "Cd-114"}, {361.50, "Cd-114"}},
       {{351.0, "73Ge 364.03->13.28 ADOPTED-ONLY"}}},
      // The 708.8 candidate. Cd 707.42 sits 1.38 keV away but is common to both
      // runs, so the background channel pins it -- this is the blend the
      // residual method could not separate.
      {"candidate_709",
       700,
       718,
       {{706.60, "Cd-114"}, {707.42, "Cd-114"}},
       {{708.8, "73Ge 776.66->68.75 ADOPTED-ONLY"},
        {714.37, "77Ge activation"}}}};

  for (const Region &r : regions)
    FitRegion(sig, bkg, r);
}
