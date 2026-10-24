// The CdTe Ge result moved +5.3 eV (68.7448 -> 68.7501 keV) when the waveform
// processing moved from the study's own trapezoid code (the `trap` trees) to
// Analysis-Utilities' PulseHeightProcessor (the `features` trees). Whether
// that is the scatter of reprocessing the same events, event by event on the
// five Ge beam runs:
//   match     pair the two trees' records by timestamp and keep the pairs
//             near the fit window in a cache: which records each processing
//             keeps, and why the new one rejects some;
//   events    per-event differences: the new/old height ratio and the rms of
//             what remains, against energy and rise time; the rise-time
//             difference; the migration between the budget's rise slices;
//   budget N  the nominal budget (GeBudget.cpp's nominal and curvature-free
//             fits per slice, its weighted mean, statistics scale, outlier and
//             extrapolation terms) of configuration N: old or new energies
//             with old or new rise slicing, on all or the common records;
//   replicas  pseudo-reprocessings of the old common records: each energy
//             moved by the measured residual and each record moved to a new
//             slice with the measured migration probabilities, then the
//             nominal budget; `replicas_joint` instead gives each record the
//             (residual, new slice) pair of a random record of the same old
//             slice and energy band;
//   summary   the replica distributions against the observed shift, and the
//             per-slice values of every configuration.
// The cache, the budgets and the replicas go to RESULT_DIR; the figures to
// plots/studies/psi_cdte/reprocessing_comparison.
//
// Usage (in order; the budgets and replicas can run in parallel):
//   root -b -q 'ReprocessingComparison.cpp+("match")'
//   root -b -q 'ReprocessingComparison.cpp+("events")'
//   root -b -q 'ReprocessingComparison.cpp+("budget", <0..6>)'
//   root -b -q 'ReprocessingComparison.cpp+("replicas", <first>, <n>)'
//   root -b -q 'ReprocessingComparison.cpp+("replicas_joint", <first>, <n>)'
//   root -b -q 'ReprocessingComparison.cpp+("summary")'
#include "CdTeAnalysis.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>

namespace CdTeReprocessing {

using namespace CdTe;

const char *const TRAP_DIR = "/labdata/ANSG/PSI/processed/CdTe/trap/";
const char *const RESULT_DIR =
    "/labdata/ANSG/PSI/processed/CdTe/reprocessing_comparison/";
// The old trees' trapezoid of rise 200 and flat top 100 samples.
const Int_t OLD_SHAPE = 1;

// Records kept in the cache: either energy in this range (keV, Am-only
// scale), which holds the in-run fit (62-89 keV) and the capture line's
// (90-99 keV).
const Double_t KEEP_LO = 55.0, KEEP_HI = 105.0;
const Double_t FIT_LO = 62.0, FIT_HI = 89.0;
// The budget combines the slices up to 320 ns, its statistics scale those up
// to 300 ns (GeBudget.cpp).
const Int_t N_BUDGET_SLICES = 5;
const Int_t N_SCALE_SLICES = 4;
const Double_t BUDGET_MAX_RISE_NS = 320.0;

// The escape tie (Te Ka / Cd Kb) each processing's NoGeEscapes.cpp gave;
// Cd Kb/Ka is Geant4's in both.
const Double_t TE_CDKB_OLD = 0.333, TE_CDKB_NEW = 0.353;
const Double_t KB_KA_NOMINAL = 0.275;

// The new processing's crop (pre and post samples around its trigger) and
// the first sample its trapezoid reads, for the rejection categories.
const Int_t NEW_PRE_SAMPLES = 400;
const Int_t NEW_PEAKING_AFTER_TRIGGER = 290;
const Int_t RECORD_SAMPLES = 2000;

// A record of the cache: heights in each tree's own units, rises in samples
// (each tree's own value, so a slice is decided exactly as FillRise does).
struct Record {
  Int_t run;
  Float_t h_old, rise_old, h_new, rise_new;
  Bool_t in_new;
  Short_t trigger_old;
};

Double_t EnergyOld(const Record &r) { return r.h_old * KEV_PER_UNIT; }
Double_t EnergyNew(const Record &r) { return r.h_new * KEV_PER_UNIT; }

// FillRise's selection: the slice of a record with this height and rise, -1
// when it enters none.
Int_t SliceOf(Float_t height, Float_t rise) {
  if (height <= 0 || rise <= 0)
    return -1;
  Double_t rise_ns = SAMPLE_NS * rise;
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    if (rise_ns >= RISE_SLICE_LO[s] && rise_ns < RISE_SLICE_HI[s])
      return s;
  return -1;
}

Int_t SliceOld(const Record &r) { return SliceOf(r.h_old, r.rise_old); }
Int_t SliceNew(const Record &r) {
  if (!r.in_new)
    return -1;
  return SliceOf(r.h_new, r.rise_new);
}

// In the fit window and one of the budget's slices.
Bool_t InBudgetWindow(Double_t e, Int_t slice) {
  return e >= FIT_LO && e < FIT_HI && slice >= 0 &&
         RISE_SLICE_HI[slice] <= BUDGET_MAX_RISE_NS;
}

// Both processings keep the record and give it a slice.
Bool_t IsCommon(const Record &r) {
  return SliceOld(r) >= 0 && SliceNew(r) >= 0;
}

TString CacheFile() { return TString(RESULT_DIR) + "events.root"; }

TString SliceName(Int_t s) {
  return Form("%.0f-%.0f", RISE_SLICE_LO[s], RISE_SLICE_HI[s]);
}

// ------------------------------------------------------------------ match

struct MatchCounts {
  Long64_t old_records = 0, new_records = 0, matched = 0, new_unmatched = 0;
  Long64_t old_window = 0, new_window = 0, both_window = 0;
  Long64_t old_only_rejected_early = 0, old_only_rejected_late = 0;
  Long64_t old_only_rejected_other = 0, old_only_energy = 0, old_only_rise = 0;
  Long64_t new_only_invalid = 0, new_only_energy = 0, new_only_rise = 0;

  void Add(const MatchCounts &o) {
    old_records += o.old_records;
    new_records += o.new_records;
    matched += o.matched;
    new_unmatched += o.new_unmatched;
    old_window += o.old_window;
    new_window += o.new_window;
    both_window += o.both_window;
    old_only_rejected_early += o.old_only_rejected_early;
    old_only_rejected_late += o.old_only_rejected_late;
    old_only_rejected_other += o.old_only_rejected_other;
    old_only_energy += o.old_only_energy;
    old_only_rise += o.old_only_rise;
    new_only_invalid += o.new_only_invalid;
    new_only_energy += o.new_only_energy;
    new_only_rise += o.new_only_rise;
  }
};

// Where a record of one processing in the budget window went in the other:
// nowhere (rejected), out of the fit window in energy, or out of the
// budget's rise slices.
void CountWindow(const Record &r, MatchCounts &c) {
  Int_t s_old = SliceOld(r), s_new = SliceNew(r);
  Bool_t old_in = InBudgetWindow(EnergyOld(r), s_old);
  Bool_t new_in = r.in_new && InBudgetWindow(EnergyNew(r), s_new);
  if (old_in)
    c.old_window++;
  if (new_in)
    c.new_window++;
  if (old_in && new_in)
    c.both_window++;
  if (old_in && !new_in) {
    if (!r.in_new) {
      if (r.trigger_old < NEW_PRE_SAMPLES)
        c.old_only_rejected_early++;
      else if (r.trigger_old > RECORD_SAMPLES - NEW_PEAKING_AFTER_TRIGGER)
        c.old_only_rejected_late++;
      else
        c.old_only_rejected_other++;
    } else if (EnergyNew(r) < FIT_LO || EnergyNew(r) >= FIT_HI) {
      c.old_only_energy++;
    } else {
      c.old_only_rise++;
    }
  }
  if (new_in && !old_in) {
    if (r.h_old <= 0 || r.rise_old <= 0)
      c.new_only_invalid++;
    else if (EnergyOld(r) < FIT_LO || EnergyOld(r) >= FIT_HI)
      c.new_only_energy++;
    else
      c.new_only_rise++;
  }
}

void PrintCounts(const TString &label, const MatchCounts &c) {
  std::cout << "MATCH " << std::left << std::setw(40) << label << std::right
            << " records old " << c.old_records << " new " << c.new_records
            << ", new matched " << c.matched << " (unmatched "
            << c.new_unmatched << ")" << std::endl;
  std::cout << "MATCH " << std::left << std::setw(40) << label << std::right
            << " window 62-89 keV, rise < 320 ns: old " << c.old_window
            << " new " << c.new_window << " both " << c.both_window << " ("
            << std::fixed << std::setprecision(2)
            << 100.0 * c.both_window / c.old_window << "% of old, "
            << 100.0 * c.both_window / c.new_window << "% of new)" << std::endl;
  std::cout << "MATCH " << std::left << std::setw(40) << label << std::right
            << "   old only: rejected by the new processing "
            << c.old_only_rejected_early + c.old_only_rejected_late +
                   c.old_only_rejected_other
            << " (old trigger before sample " << NEW_PRE_SAMPLES << ": "
            << c.old_only_rejected_early << ", after "
            << RECORD_SAMPLES - NEW_PEAKING_AFTER_TRIGGER << ": "
            << c.old_only_rejected_late << ", other "
            << c.old_only_rejected_other << "), new energy outside "
            << c.old_only_energy << ", new rise outside " << c.old_only_rise
            << std::endl;
  std::cout << "MATCH " << std::left << std::setw(40) << label << std::right
            << "   new only: old invalid " << c.new_only_invalid
            << ", old energy outside " << c.new_only_energy
            << ", old rise outside " << c.new_only_rise << std::endl;
}

// The new tree is the old one's records in the same order less those the
// new processing rejects, so one pass pairs them.
MatchCounts MatchRun(Int_t run_index, const TString &run, TTree *cache,
                     Record &out) {
  MatchCounts c;
  TFile old_file(TString(TRAP_DIR) + FILE_PREFIX + run + ".root");
  TTree *old_tree = static_cast<TTree *>(old_file.Get("trap"));
  TFile new_file(FeaturesFile(run));
  TTree *new_tree = static_cast<TTree *>(new_file.Get("features"));

  ULong64_t new_ts = 0;
  Float_t new_height = 0, new_rise = 0;
  new_tree->SetBranchStatus("*", 0);
  new_tree->SetBranchStatus("timestamp", 1);
  new_tree->SetBranchStatus("trapezoid_height", 1);
  new_tree->SetBranchStatus("rise_time", 1);
  new_tree->SetBranchAddress("timestamp", &new_ts);
  new_tree->SetBranchAddress("trapezoid_height", &new_height);
  new_tree->SetBranchAddress("rise_time", &new_rise);

  ULong64_t old_ts = 0;
  Float_t old_energy[4], old_rise = 0;
  Short_t old_trigger = 0;
  old_tree->SetBranchStatus("*", 0);
  old_tree->SetBranchStatus("timestamp", 1);
  old_tree->SetBranchStatus("energy", 1);
  old_tree->SetBranchStatus("rise", 1);
  old_tree->SetBranchStatus("trigger", 1);
  old_tree->SetBranchAddress("timestamp", &old_ts);
  old_tree->SetBranchAddress("energy", old_energy);
  old_tree->SetBranchAddress("rise", &old_rise);
  old_tree->SetBranchAddress("trigger", &old_trigger);

  c.old_records = old_tree->GetEntries();
  c.new_records = new_tree->GetEntries();
  Long64_t j = 0;
  if (c.new_records > 0)
    new_tree->GetEntry(0);
  for (Long64_t i = 0; i < c.old_records; i++) {
    old_tree->GetEntry(i);
    out.run = run_index;
    out.h_old = old_energy[OLD_SHAPE];
    out.rise_old = old_rise;
    out.trigger_old = old_trigger;
    out.in_new = (j < c.new_records && new_ts == old_ts);
    out.h_new = out.in_new ? new_height : -1;
    out.rise_new = out.in_new ? new_rise : -1;
    if (out.in_new) {
      c.matched++;
      j++;
      if (j < c.new_records)
        new_tree->GetEntry(j);
    }
    CountWindow(out, c);
    Double_t e_old = EnergyOld(out), e_new = EnergyNew(out);
    Bool_t keep = (e_old >= KEEP_LO && e_old < KEEP_HI) ||
                  (out.in_new && e_new >= KEEP_LO && e_new < KEEP_HI);
    if (keep)
      cache->Fill();
  }
  c.new_unmatched = c.new_records - j;
  return c;
}

void Match() {
  gSystem->mkdir(RESULT_DIR, kTRUE);
  TFile file(CacheFile(), "RECREATE");
  TTree *cache = new TTree("events", "old and new CdTe processing, paired");
  Record r;
  cache->Branch("run", &r.run, "run/I");
  cache->Branch("h_old", &r.h_old, "h_old/F");
  cache->Branch("rise_old", &r.rise_old, "rise_old/F");
  cache->Branch("h_new", &r.h_new, "h_new/F");
  cache->Branch("rise_new", &r.rise_new, "rise_new/F");
  cache->Branch("in_new", &r.in_new, "in_new/O");
  cache->Branch("trigger_old", &r.trigger_old, "trigger_old/S");
  std::vector<TString> runs = GeRuns();
  MatchCounts total;
  for (size_t i = 0; i < runs.size(); i++) {
    MatchCounts c = MatchRun((Int_t)i, runs[i], cache, r);
    PrintCounts(runs[i], c);
    total.Add(c);
  }
  PrintCounts("all five Ge runs", total);
  file.cd();
  cache->Write();
  std::cout << "Wrote " << CacheFile() << " (" << cache->GetEntries()
            << " records)" << std::endl;
}

std::vector<Record> LoadCache() {
  std::vector<Record> records;
  TFile file(CacheFile());
  TTree *cache = static_cast<TTree *>(file.Get("events"));
  Record r;
  cache->SetBranchAddress("run", &r.run);
  cache->SetBranchAddress("h_old", &r.h_old);
  cache->SetBranchAddress("rise_old", &r.rise_old);
  cache->SetBranchAddress("h_new", &r.h_new);
  cache->SetBranchAddress("rise_new", &r.rise_new);
  cache->SetBranchAddress("in_new", &r.in_new);
  cache->SetBranchAddress("trigger_old", &r.trigger_old);
  const Long64_t n = cache->GetEntries();
  records.reserve(n);
  for (Long64_t i = 0; i < n; i++) {
    cache->GetEntry(i);
    records.push_back(r);
  }
  return records;
}

// ------------------------------------------------------------- statistics

Double_t Median(std::vector<Double_t> v) {
  if (v.empty())
    return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

// The rms of a Gaussian core: 1.4826 times the median absolute deviation.
Double_t RobustSigma(const std::vector<Double_t> &v, Double_t centre) {
  std::vector<Double_t> deviations;
  for (size_t i = 0; i < v.size(); i++)
    deviations.push_back(TMath::Abs(v[i] - centre));
  return 1.4826 * Median(deviations);
}

Double_t Rms(const std::vector<Double_t> &v, Double_t centre) {
  if (v.empty())
    return 0;
  Double_t sum2 = 0;
  for (size_t i = 0; i < v.size(); i++)
    sum2 += (v[i] - centre) * (v[i] - centre);
  return TMath::Sqrt(sum2 / v.size());
}

Double_t FractionBeyond(const std::vector<Double_t> &v, Double_t centre,
                        Double_t limit) {
  if (v.empty())
    return 0;
  Long64_t n = 0;
  for (size_t i = 0; i < v.size(); i++)
    if (TMath::Abs(v[i] - centre) > limit)
      n++;
  return Double_t(n) / v.size();
}

// The energy bands the residual is measured in (keV, old energy).
const Int_t N_BANDS = 5;
const Double_t BAND_LO[N_BANDS] = {55, 62, 70, 80, 89};
const Double_t BAND_HI[N_BANDS] = {62, 70, 80, 89, 105};

Int_t BandOf(Double_t e) {
  for (Int_t b = 0; b < N_BANDS; b++)
    if (e >= BAND_LO[b] && e < BAND_HI[b])
      return b;
  return -1;
}

// The new/old gain of each old slice: the median height ratio in the fit
// window, over the common records.
std::vector<Double_t> SliceGains(const std::vector<Record> &records) {
  std::vector<std::vector<Double_t>> ratios(N_RISE_SLICES);
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Double_t e = EnergyOld(r);
    if (e < FIT_LO || e >= FIT_HI)
      continue;
    ratios[SliceOld(r)].push_back(r.h_new / r.h_old);
  }
  std::vector<Double_t> gains;
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    gains.push_back(ratios[s].empty() ? 1.0 : Median(ratios[s]));
  return gains;
}

// What the replicas draw from: per old slice and energy band, the residual
// (new minus gain times old, keV) of every common record and the new slice
// it went to (N_RISE_SLICES: above the last slice's range).
struct ReprocessingModel {
  std::vector<Double_t> gain;
  std::vector<std::vector<Double_t>> residual;  // [slice * N_BANDS + band]
  std::vector<std::vector<Int_t>> new_slice;    // same index
  std::vector<std::vector<Double_t>> migration; // [old][new], row-normalised
  std::vector<Double_t> sigma;                  // robust, keV, same index
};

Int_t ModelIndex(Int_t slice, Int_t band) { return slice * N_BANDS + band; }

ReprocessingModel BuildModel(const std::vector<Record> &records) {
  ReprocessingModel m;
  m.gain = SliceGains(records);
  m.residual.resize(N_RISE_SLICES * N_BANDS);
  m.new_slice.resize(N_RISE_SLICES * N_BANDS);
  m.migration.assign(N_RISE_SLICES, std::vector<Double_t>(N_RISE_SLICES, 0));
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Int_t band = BandOf(EnergyOld(r));
    if (band < 0)
      continue;
    Int_t s_old = SliceOld(r), s_new = SliceNew(r);
    Int_t k = ModelIndex(s_old, band);
    m.residual[k].push_back(EnergyNew(r) - m.gain[s_old] * EnergyOld(r));
    m.new_slice[k].push_back(s_new);
    m.migration[s_old][s_new] += 1;
  }
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    Double_t row = 0;
    for (Int_t t = 0; t < N_RISE_SLICES; t++)
      row += m.migration[s][t];
    for (Int_t t = 0; t < N_RISE_SLICES; t++)
      m.migration[s][t] = (row > 0) ? m.migration[s][t] / row : 0;
  }
  for (size_t k = 0; k < m.residual.size(); k++)
    m.sigma.push_back(RobustSigma(m.residual[k], Median(m.residual[k])));
  return m;
}

// ----------------------------------------------------------------- events

TH1D *NewProfileFrame(const TString &name, Double_t lo, Double_t hi) {
  TH1D *frame = new TH1D(name, "", 10, lo, hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(0);
  return frame;
}

const Int_t SLICE_COLORS[N_BUDGET_SLICES] = {kBlue + 1, kRed + 1, kGreen + 2,
                                             kOrange + 7, kMagenta + 1};

// One series per budget slice against old energy: the gain or the robust
// residual width per 2 keV.
void DrawPerSliceVsEnergy(const std::vector<Record> &records,
                          const std::vector<Double_t> &gains, Bool_t width,
                          const TString &name) {
  const Double_t e_lo = KEEP_LO, e_hi = KEEP_HI, step = 2.5;
  const Int_t n_bins = (Int_t)((e_hi - e_lo) / step + 0.5);
  std::vector<std::vector<std::vector<Double_t>>> values(
      N_BUDGET_SLICES, std::vector<std::vector<Double_t>>(n_bins));
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Int_t s = SliceOld(r);
    if (s >= N_BUDGET_SLICES)
      continue;
    Double_t e = EnergyOld(r);
    Int_t b = (Int_t)((e - e_lo) / step);
    if (b < 0 || b >= n_bins)
      continue;
    if (width)
      values[s][b].push_back(1000.0 * (EnergyNew(r) - gains[s] * e));
    else
      values[s][b].push_back(r.h_new / r.h_old);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TH1D *frame = NewProfileFrame("frame_" + name, e_lo, e_hi);
  PlottingUtils::ConfigureHistogram(
      frame, kBlack,
      width ? ";Old energy [keV];Residual width (robust) [eV]"
            : ";Old energy [keV];New / old pulse height (median)");
  frame->SetMinimum(width ? 0 : 0.9975);
  frame->SetMaximum(width ? 200 : 1.0035);
  frame->Draw("AXIS");
  TLegend *leg = PlottingUtils::AddLegend(0.55, 0.92, 0.65, 0.92);
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
    TGraphErrors *g = new TGraphErrors();
    for (Int_t b = 0; b < n_bins; b++) {
      const std::vector<Double_t> &v = values[s][b];
      if (v.size() < 200)
        continue;
      Double_t centre = Median(v);
      Double_t sigma = RobustSigma(v, centre);
      Double_t x = e_lo + (b + 0.5) * step + 0.2 * (s - 2);
      Int_t k = g->GetN();
      if (width) {
        g->SetPoint(k, x, sigma);
        g->SetPointError(k, 0, sigma / TMath::Sqrt(2.0 * v.size()));
      } else {
        g->SetPoint(k, x, centre);
        g->SetPointError(k, 0, 1.2533 * sigma / TMath::Sqrt(v.size()));
      }
    }
    PlottingUtils::ConfigureGraph(g, SLICE_COLORS[s]);
    g->SetMarkerStyle(20 + s);
    g->Draw("P SAME");
    leg->AddEntry(g, "Old rise " + SliceName(s) + " ns", "pe");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

void DrawResidualPlanes(const std::vector<Record> &records,
                        const std::vector<Double_t> &gains) {
  TH2D *vs_energy =
      new TH2D("residual_vs_energy", "", 250, KEEP_LO, KEEP_HI, 160, -400, 400);
  TH2D *vs_rise =
      new TH2D("residual_vs_rise", "", 90, 200, 470, 160, -400, 400);
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Int_t s = SliceOld(r);
    Double_t e = EnergyOld(r);
    Double_t d = 1000.0 * (EnergyNew(r) - gains[s] * e);
    if (s < N_BUDGET_SLICES)
      vs_energy->Fill(e, d);
    if (e >= FIT_LO && e < FIT_HI)
      vs_rise->Fill(SAMPLE_NS * r.rise_old, d);
  }
  TCanvas *c1 = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      vs_energy, c1,
      ";Old energy [keV];New #minus gain #times old [eV];Counts / 200 eV / "
      "5 eV");
  PlottingUtils::SaveFigure(c1, "residual_vs_energy", "",
                            PlotSaveOptions::kLINEAR);
  TCanvas *c2 = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      vs_rise, c2,
      ";Old rise time 10-90% [ns];New #minus gain #times old [eV];Counts / "
      "3 ns / 5 eV");
  PlottingUtils::SaveFigure(c2, "residual_vs_rise", "",
                            PlotSaveOptions::kLINEAR);
}

Double_t GaussPlusFlat(Double_t *x, Double_t *p) {
  Double_t z = (x[0] - p[1]) / p[2];
  return p[0] * TMath::Exp(-0.5 * z * z) + p[3];
}

// The residual distribution in the fit window and the budget's slices, with
// a Gaussian on a flat level over the core.
void DrawResidualFit(const std::vector<Record> &records,
                     const std::vector<Double_t> &gains) {
  const Double_t lo = -300, hi = 300, bin_ev = 5;
  TH1D *h =
      new TH1D("residual_fit_window", "", (Int_t)((hi - lo) / bin_ev), lo, hi);
  h->SetDirectory(nullptr);
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Int_t s = SliceOld(r);
    Double_t e = EnergyOld(r);
    if (!InBudgetWindow(e, s))
      continue;
    h->Fill(1000.0 * (EnergyNew(r) - gains[s] * e));
  }
  PlottingUtils::ConfigureHistogram(
      h, kBlack,
      Form(";New #minus gain #times old energy [eV];Counts / %.0f eV", bin_ev));
  TF1 *f = new TF1("residual_gauss", GaussPlusFlat, -150, 150, 4);
  f->SetParameters(h->GetMaximum(), 0, 30, 1);
  f->SetParLimits(3, 0, 1e9);
  h->Fit(f, "QLRN0");
  std::cout << "EVENTS residual in the budget window: Gaussian core mean "
            << std::fixed << std::setprecision(1) << f->GetParameter(1)
            << " eV, sigma " << f->GetParameter(2) << " +/- "
            << f->GetParError(2) << " eV, chi2/ndf " << std::setprecision(1)
            << f->GetChisquare() << "/" << f->GetNDF() << ", entries "
            << std::setprecision(0) << h->GetEntries() << std::endl;
  const Int_t n_pts = 600;
  TGraph *total = new TGraph(n_pts);
  for (Int_t i = 0; i < n_pts; i++) {
    Double_t x = -150 + 300.0 * i / (n_pts - 1.0);
    total->SetPoint(i, x, f->Eval(x));
  }
  total->SetLineColor(kRed + 1);
  total->SetLineWidth(3);
  std::vector<TGraph *> components;
  PlottingUtils::PlotFitWithResiduals(h, total, components, -150, 150,
                                      "residual_fit_window", "",
                                      "62-89 keV, old rise < 320 ns", kTRUE);
}

void DrawRiseDifference(const std::vector<Record> &records) {
  const Double_t lo = -40, hi = 40, bin_ns = 0.5;
  TH1D *h =
      new TH1D("rise_difference", "", (Int_t)((hi - lo) / bin_ns), lo, hi);
  h->SetDirectory(nullptr);
  std::vector<Double_t> d;
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Double_t e = EnergyOld(r);
    if (e < FIT_LO || e >= FIT_HI)
      continue;
    Double_t diff = SAMPLE_NS * (r.rise_new - r.rise_old);
    h->Fill(diff);
    d.push_back(diff);
  }
  Double_t centre = Median(d);
  std::cout << "EVENTS rise new - old (62-89 keV): median " << std::fixed
            << std::setprecision(2) << centre << " ns, robust sigma "
            << RobustSigma(d, centre) << " ns, rms " << Rms(d, centre)
            << " ns, beyond 10 ns " << std::setprecision(4)
            << 100.0 * FractionBeyond(d, centre, 10.0) << "%" << std::endl;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  PlottingUtils::ConfigureAndDrawHistogram(
      h, kBlack,
      Form(";New #minus old rise time 10-90%% [ns];Counts / %.1f ns", bin_ns));
  PlottingUtils::SaveFigure(canvas, "rise_difference", "",
                            PlotSaveOptions::kLINEAR);
}

// Each processing's rise-time distribution in the fit window, against the
// slice edges.
void DrawRiseDistributions(const std::vector<Record> &records) {
  const Double_t lo = 240, hi = 400, bin_ns = 0.5;
  const Int_t n_bins = (Int_t)((hi - lo) / bin_ns + 0.5);
  TH1D *h_old = new TH1D("rise_old", "", n_bins, lo, hi);
  TH1D *h_new = new TH1D("rise_new", "", n_bins, lo, hi);
  h_old->SetDirectory(nullptr);
  h_new->SetDirectory(nullptr);
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    Double_t e_old = EnergyOld(r), e_new = EnergyNew(r);
    if (r.h_old > 0 && r.rise_old > 0 && e_old >= FIT_LO && e_old < FIT_HI)
      h_old->Fill(SAMPLE_NS * r.rise_old);
    if (r.in_new && r.h_new > 0 && r.rise_new > 0 && e_new >= FIT_LO &&
        e_new < FIT_HI)
      h_new->Fill(SAMPLE_NS * r.rise_new);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(
      h_old, kBlack, Form(";Rise time 10-90%% [ns];Counts / %.1f ns", bin_ns));
  PlottingUtils::ConfigureHistogram(h_new, kRed + 1);
  h_old->SetMaximum(1.25 *
                    TMath::Max(h_old->GetMaximum(), h_new->GetMaximum()));
  h_old->Draw("HIST");
  h_new->Draw("HIST SAME");
  for (Int_t s = 1; s < N_RISE_SLICES; s++) {
    if (RISE_SLICE_LO[s] >= hi)
      continue;
    TLine *edge =
        new TLine(RISE_SLICE_LO[s], 0, RISE_SLICE_LO[s], h_old->GetMaximum());
    edge->SetLineStyle(2);
    edge->SetLineColor(kGray + 2);
    edge->Draw();
  }
  TLegend *leg = PlottingUtils::AddLegend(0.6, 0.92, 0.75, 0.92);
  leg->AddEntry(h_old, "Old processing", "l");
  leg->AddEntry(h_new, "New processing", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "rise_distributions", "",
                            PlotSaveOptions::kLINEAR);
}

// Old slice against new slice over the records in 62-99 keV (old energy):
// rows the old slices, then old invalid; columns the new slices, then new
// invalid or above 450 ns, then rejected. Row-normalised, in percent.
void DrawMigration(const std::vector<Record> &records) {
  const Int_t n_rows = N_RISE_SLICES + 1, n_cols = N_RISE_SLICES + 2;
  TH2D *h = new TH2D("migration", "", n_cols, 0, n_cols, n_rows, 0, n_rows);
  h->SetDirectory(nullptr);
  std::vector<Double_t> row_total(n_rows, 0);
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    Double_t e = (r.h_old > 0) ? EnergyOld(r) : EnergyNew(r);
    if (e < FIT_LO || e >= 99.0)
      continue;
    Int_t s_old = SliceOld(r);
    Int_t row = (s_old >= 0) ? s_old : N_RISE_SLICES;
    Int_t col;
    if (!r.in_new)
      col = N_RISE_SLICES + 1;
    else if (SliceNew(r) < 0)
      col = N_RISE_SLICES;
    else
      col = SliceNew(r);
    if (row == N_RISE_SLICES && col >= N_RISE_SLICES)
      continue;
    h->Fill(col + 0.5, row + 0.5);
    row_total[row] += 1;
  }
  std::cout << "MIGRATION rows old slice, columns new slice (% of the row), "
               "62-99 keV:"
            << std::endl;
  for (Int_t row = 0; row < n_rows; row++) {
    TString line =
        (row < N_RISE_SLICES) ? SliceName(row) + " ns" : TString("old none");
    line = Form("%-12s n %9.0f:", line.Data(), row_total[row]);
    for (Int_t col = 0; col < n_cols; col++) {
      Double_t count = h->GetBinContent(col + 1, row + 1);
      Double_t pct = (row_total[row] > 0) ? 100.0 * count / row_total[row] : 0;
      h->SetBinContent(col + 1, row + 1, pct);
      line += Form(" %7.3f", pct);
    }
    std::cout << "MIGRATION " << line << std::endl;
  }
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    h->GetXaxis()->SetBinLabel(s + 1, SliceName(s));
    h->GetYaxis()->SetBinLabel(s + 1, SliceName(s));
  }
  h->GetXaxis()->SetBinLabel(N_RISE_SLICES + 1, "none");
  h->GetXaxis()->SetBinLabel(N_RISE_SLICES + 2, "rejected");
  h->GetYaxis()->SetBinLabel(N_RISE_SLICES + 1, "none");
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::Configure2DHistogram(
      h, canvas, ";New rise slice [ns];Old rise slice [ns];% of the old slice");
  canvas->SetLogz(1);
  canvas->SetLeftMargin(0.16);
  canvas->SetBottomMargin(0.16);
  h->SetMinimum(0.01);
  h->SetMaximum(100);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetPaintTextFormat(".2f");
  h->SetMarkerSize(1.3);
  h->Draw("COLZ TEXT");
  PlottingUtils::SaveFigure(canvas, "migration_matrix", "",
                            PlotSaveOptions::kLINEAR);
}

void Events() {
  std::vector<Record> records = LoadCache();
  ReprocessingModel m = BuildModel(records);
  std::cout << "EVENTS gain (new/old median height ratio, 62-89 keV) per old "
               "slice:";
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    std::cout << " " << SliceName(s) << " " << std::fixed
              << std::setprecision(5) << m.gain[s];
  std::cout << std::endl;
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    for (Int_t b = 0; b < N_BANDS; b++) {
      const std::vector<Double_t> &v = m.residual[ModelIndex(s, b)];
      if (v.empty())
        continue;
      Double_t centre = Median(v);
      Double_t sigma = RobustSigma(v, centre);
      std::cout << "EVENTS residual old rise " << std::setw(7) << SliceName(s)
                << " ns, " << std::setprecision(0) << BAND_LO[b] << "-"
                << BAND_HI[b] << " keV: n " << v.size() << ", median "
                << std::setprecision(1) << 1000 * centre << " eV, robust sigma "
                << 1000 * sigma << " eV, rms " << 1000 * Rms(v, centre)
                << " eV, beyond 5 sigma " << std::setprecision(3)
                << 100.0 * FractionBeyond(v, centre, 5 * sigma) << "%"
                << std::endl;
    }
  }
  DrawPerSliceVsEnergy(records, m.gain, kFALSE, "gain_vs_energy");
  DrawPerSliceVsEnergy(records, m.gain, kTRUE, "residual_width_vs_energy");
  DrawResidualPlanes(records, m.gain);
  DrawResidualFit(records, m.gain);
  DrawRiseDifference(records);
  DrawRiseDistributions(records);
  DrawMigration(records);
}

// ---------------------------------------------------------------- budgets

// A record's energy (keV) and slice under one configuration; slice -1 leaves
// it out.
struct Entry {
  Float_t e;
  Int_t slice;
};

// The nominal budget of one configuration: the per-slice nominal and
// curvature-free fits, and GeBudget.cpp's combination of the slices up to
// 320 ns.
struct NominalBudget {
  Double_t value[N_BUDGET_SLICES], error[N_BUDGET_SLICES];
  Double_t curved[N_BUDGET_SLICES], curved_error[N_BUDGET_SLICES];
  Double_t mean = 0, weighted_err = 0, scale = 1, stat = 0, outlier = 0;
  Double_t extrapolation = 0, curvature_shift = 0;
  Int_t i_hi = 0, i_lo = 0;
};

Double_t WeightedMean(const Double_t *v, const Double_t *e, Int_t n,
                      Int_t skip1, Int_t skip2, Double_t &err) {
  Double_t sw = 0, swx = 0;
  for (Int_t i = 0; i < n; i++) {
    if (i == skip1 || i == skip2)
      continue;
    sw += 1.0 / (e[i] * e[i]);
    swx += v[i] / (e[i] * e[i]);
  }
  err = 1.0 / TMath::Sqrt(sw);
  return swx / sw;
}

void Combine(NominalBudget &b) {
  const Int_t n = N_BUDGET_SLICES;
  b.mean = WeightedMean(b.value, b.error, n, -1, -1, b.weighted_err);
  Double_t scale_err = 0;
  Double_t scale_mean =
      WeightedMean(b.value, b.error, N_SCALE_SLICES, -1, -1, scale_err);
  Double_t chi2 = 0;
  for (Int_t i = 0; i < N_SCALE_SLICES; i++)
    chi2 += TMath::Power((b.value[i] - scale_mean) / b.error[i], 2);
  b.scale = TMath::Sqrt(TMath::Max(chi2 / (N_SCALE_SLICES - 1), 1.0));
  b.stat = b.weighted_err * b.scale;
  b.i_hi = 0;
  b.i_lo = 0;
  for (Int_t i = 1; i < n; i++) {
    if (b.value[i] > b.value[b.i_hi])
      b.i_hi = i;
    if (b.value[i] < b.value[b.i_lo])
      b.i_lo = i;
  }
  Double_t unused = 0;
  b.outlier = TMath::Max(
      TMath::Abs(WeightedMean(b.value, b.error, n, b.i_hi, -1, unused) -
                 b.mean),
      TMath::Abs(WeightedMean(b.value, b.error, n, b.i_lo, -1, unused) -
                 b.mean));
  b.outlier = TMath::Max(
      b.outlier,
      TMath::Abs(WeightedMean(b.value, b.error, n, b.i_hi, b.i_lo, unused) -
                 b.mean));
  b.curvature_shift =
      WeightedMean(b.curved, b.error, n, -1, -1, unused) - b.mean;
  Double_t curved_err = 0;
  WeightedMean(b.curved, b.curved_error, n, -1, -1, curved_err);
  Double_t curvature_cost = TMath::Sqrt(TMath::Max(
      curved_err * curved_err - b.weighted_err * b.weighted_err, 0.0));
  b.extrapolation = TMath::Max(TMath::Abs(b.curvature_shift), curvature_cost);
}

TH1D *FillSlice(const std::vector<Entry> &entries, Int_t slice,
                const TString &name) {
  TH1D *h = new TH1D(name, Form(";Energy [keV];Counts / %.0f eV", 1000.0 * BIN),
                     (Int_t)(120.0 / BIN + 0.5), 0, 120);
  h->SetDirectory(nullptr);
  for (size_t i = 0; i < entries.size(); i++)
    if (entries[i].slice == slice)
      h->Fill(entries[i].e);
  return h;
}

// GeBudget.cpp's nominal settings, with the given escape tie, and the
// curvature-free variant.
void SetNominal(Double_t te_cdkb, Bool_t free_curvature) {
  g_lo = 62.0;
  g_hi = 89.0;
  g_ge_own_tail = kFALSE;
  g_free_scale2 = free_curvature;
  SetEscapeTie(KB_KA_NOMINAL, te_cdkb);
  ResetKbGroups();
}

NominalBudget RunBudget(const std::vector<Entry> &entries, Double_t te_cdkb,
                        const TString &tag, Bool_t verbose) {
  NominalBudget b;
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
    TH1D *h = FillSlice(entries, s, Form("h_%s_%d", tag.Data(), s));
    g_capture_meas = FitCaptureLine(h, Form("%s_%d", tag.Data(), s)).mu;
    SetNominal(te_cdkb, kFALSE);
    TF1 *f = FitInRun(h, Form("%s_%d_nominal", tag.Data(), s));
    b.value[s] = f->GetParameter(Q_GE_E);
    b.error[s] = f->GetParError(Q_GE_E);
    delete f;
    SetNominal(te_cdkb, kTRUE);
    f = FitInRun(h, Form("%s_%d_curved", tag.Data(), s));
    b.curved[s] = f->GetParameter(Q_GE_E);
    b.curved_error[s] = f->GetParError(Q_GE_E);
    delete f;
    SetNominal(te_cdkb, kFALSE);
    if (verbose)
      std::cout << "BUDGET " << tag << " slice " << SliceName(s)
                << " ns: entries " << std::setprecision(0) << std::fixed
                << h->Integral(h->FindBin(FIT_LO), h->FindBin(FIT_HI))
                << ", Ge " << std::setprecision(4) << b.value[s] << " +/- "
                << std::setprecision(1) << 1000 * b.error[s]
                << " eV, curvature free " << std::setprecision(4) << b.curved[s]
                << " +/- " << std::setprecision(1) << 1000 * b.curved_error[s]
                << " eV" << std::endl;
    delete h;
  }
  Combine(b);
  return b;
}

void WriteBudget(std::ofstream &out, const TString &label,
                 const NominalBudget &b) {
  out << std::fixed << std::setprecision(7) << label;
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
    out << " " << b.value[s] << " " << b.error[s] << " " << b.curved[s] << " "
        << b.curved_error[s];
  out << std::endl;
}

Bool_t ReadBudget(std::ifstream &in, TString &label, NominalBudget &b) {
  std::string word;
  if (!(in >> word))
    return kFALSE;
  label = word.c_str();
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
    in >> b.value[s] >> b.error[s] >> b.curved[s] >> b.curved_error[s];
  Combine(b);
  return kTRUE;
}

void PrintBudget(const TString &label, const NominalBudget &b) {
  std::cout << "BUDGET " << std::left << std::setw(24) << label << std::right
            << std::fixed << std::setprecision(4) << " Ge " << b.mean
            << " keV, stat " << std::setprecision(1) << 1000 * b.stat
            << " eV (x" << std::setprecision(2) << b.scale << "), outlier "
            << std::setprecision(1) << 1000 * b.outlier << " eV (highest "
            << SliceName(b.i_hi) << ", lowest " << SliceName(b.i_lo)
            << "), extrapolation " << 1000 * b.extrapolation
            << " eV (curvature-free shift " << std::showpos
            << 1000 * b.curvature_shift << std::noshowpos << ")" << std::endl;
}

// The configurations: which energies, which slicing, which records, which
// escape tie.
enum {
  C_OLD_TIE_OLD,
  C_OLD,
  C_NEW,
  C_COMMON_OLD,
  C_COMMON_NEW,
  C_NEW_E_OLD_RISE,
  C_OLD_E_NEW_RISE,
  N_CONFIGS
};
const char *const CONFIG_LABELS[N_CONFIGS] = {
    "old_tie0.333", "old",          "new",         "common_old",
    "common_new",   "newE_oldrise", "oldE_newrise"};
const char *const CONFIG_TITLES[N_CONFIGS] = {
    "Old processing, its escape tie (0.333)",
    "Old processing",
    "New processing",
    "Common records, old",
    "Common records, new",
    "New energies, old rise slices",
    "Old energies, new rise slices"};

std::vector<Entry> ConfigEntries(const std::vector<Record> &records,
                                 Int_t config) {
  std::vector<Entry> entries;
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    Entry en;
    en.slice = -1;
    en.e = 0;
    Bool_t common = IsCommon(r);
    switch (config) {
    case C_OLD_TIE_OLD:
    case C_OLD:
      en.e = EnergyOld(r);
      en.slice = SliceOld(r);
      break;
    case C_NEW:
      en.e = EnergyNew(r);
      en.slice = SliceNew(r);
      break;
    case C_COMMON_OLD:
      en.e = EnergyOld(r);
      en.slice = common ? SliceOld(r) : -1;
      break;
    case C_COMMON_NEW:
      en.e = EnergyNew(r);
      en.slice = common ? SliceNew(r) : -1;
      break;
    case C_NEW_E_OLD_RISE:
      en.e = EnergyNew(r);
      en.slice = common ? SliceOld(r) : -1;
      break;
    case C_OLD_E_NEW_RISE:
      en.e = EnergyOld(r);
      en.slice = common ? SliceNew(r) : -1;
      break;
    }
    if (en.slice >= 0 && en.slice < N_BUDGET_SLICES)
      entries.push_back(en);
  }
  return entries;
}

TString BudgetFile(Int_t config) {
  return Form("%sbudget_%d.txt", RESULT_DIR, config);
}

void Budget(Int_t config) {
  std::vector<Record> records = LoadCache();
  std::vector<Entry> entries = ConfigEntries(records, config);
  Double_t te = (config == C_OLD_TIE_OLD) ? TE_CDKB_OLD : TE_CDKB_NEW;
  NominalBudget b = RunBudget(entries, te, CONFIG_LABELS[config], kTRUE);
  PrintBudget(CONFIG_LABELS[config], b);
  std::ofstream out(BudgetFile(config).Data());
  WriteBudget(out, CONFIG_LABELS[config], b);
  std::cout << "Wrote " << BudgetFile(config) << std::endl;
}

// --------------------------------------------------------------- replicas

// One pseudo-reprocessing of the old common records. Gaussian: the old
// energy plus a Gaussian of its slice and band's robust residual width, the
// new slice drawn from the old slice's migration row. Joint: the residual
// and new slice of a random common record of the same old slice and band.
std::vector<Entry> ReplicaEntries(const std::vector<Record> &records,
                                  const ReprocessingModel &m, Bool_t joint,
                                  TRandom3 &random) {
  std::vector<Entry> entries;
  for (size_t i = 0; i < records.size(); i++) {
    const Record &r = records[i];
    if (!IsCommon(r))
      continue;
    Double_t e = EnergyOld(r);
    Int_t band = BandOf(e);
    if (band < 0)
      continue;
    Int_t s_old = SliceOld(r);
    Int_t k = ModelIndex(s_old, band);
    Entry en;
    if (joint) {
      Int_t pick = (Int_t)(random.Rndm() * m.residual[k].size());
      en.e = e + m.residual[k][pick];
      en.slice = m.new_slice[k][pick];
    } else {
      en.e = e + random.Gaus(0, m.sigma[k]);
      Double_t u = random.Rndm(), cumulative = 0;
      en.slice = N_RISE_SLICES - 1;
      for (Int_t t = 0; t < N_RISE_SLICES; t++) {
        cumulative += m.migration[s_old][t];
        if (u < cumulative) {
          en.slice = t;
          break;
        }
      }
    }
    if (en.slice < N_BUDGET_SLICES)
      entries.push_back(en);
  }
  return entries;
}

TString ReplicaFile(Bool_t joint, Int_t first) {
  return Form("%sreplicas_%s_%04d.txt", RESULT_DIR,
              joint ? "joint" : "gaussian", first);
}

void Replicas(Bool_t joint, Int_t first, Int_t count) {
  std::vector<Record> records = LoadCache();
  ReprocessingModel m = BuildModel(records);
  std::ofstream out(ReplicaFile(joint, first).Data());
  for (Int_t k = first; k < first + count; k++) {
    TRandom3 random(1000 + k);
    std::vector<Entry> entries = ReplicaEntries(records, m, joint, random);
    TString tag = Form("replica%s%d", joint ? "j" : "g", k);
    NominalBudget b = RunBudget(entries, TE_CDKB_NEW, tag, kFALSE);
    PrintBudget(tag, b);
    WriteBudget(out, tag, b);
  }
  std::cout << "Wrote " << ReplicaFile(joint, first) << std::endl;
}

std::vector<NominalBudget> ReadReplicas(Bool_t joint) {
  std::vector<NominalBudget> replicas;
  TString pattern = joint ? "replicas_joint_" : "replicas_gaussian_";
  void *dir = gSystem->OpenDirectory(RESULT_DIR);
  const char *entry;
  std::vector<TString> files;
  while ((entry = gSystem->GetDirEntry(dir))) {
    TString name = entry;
    if (name.BeginsWith(pattern) && name.EndsWith(".txt"))
      files.push_back(TString(RESULT_DIR) + name);
  }
  gSystem->FreeDirectory(dir);
  std::sort(files.begin(), files.end());
  for (size_t i = 0; i < files.size(); i++) {
    std::ifstream in(files[i].Data());
    TString label;
    NominalBudget b;
    while (ReadBudget(in, label, b))
      replicas.push_back(b);
  }
  return replicas;
}

// ---------------------------------------------------------------- summary

void MeanAndRms(const std::vector<Double_t> &v, Double_t &mean, Double_t &rms) {
  mean = 0;
  for (size_t i = 0; i < v.size(); i++)
    mean += v[i] / v.size();
  rms = 0;
  for (size_t i = 0; i < v.size(); i++)
    rms += (v[i] - mean) * (v[i] - mean) / (v.size() - 1.0);
  rms = TMath::Sqrt(rms);
}

Double_t FractionAtLeast(const std::vector<Double_t> &v, Double_t x) {
  Int_t n = 0;
  for (size_t i = 0; i < v.size(); i++)
    if (TMath::Abs(v[i]) >= TMath::Abs(x))
      n++;
  return Double_t(n) / v.size();
}

// A replica distribution (eV) with the observed values marked.
void DrawReplicaHistogram(const std::vector<Double_t> &values,
                          const std::vector<Double_t> &marks,
                          const std::vector<TString> &mark_labels,
                          const TString &x_title, const TString &kind_label,
                          Double_t lo, Double_t hi, Double_t bin_ev,
                          const TString &name) {
  TH1D *h =
      new TH1D("h_" + name, "", (Int_t)((hi - lo) / bin_ev + 0.5), lo, hi);
  h->SetDirectory(nullptr);
  for (size_t i = 0; i < values.size(); i++)
    h->Fill(values[i]);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(
      h, kBlack, Form(";%s [eV];Replicas / %.1f eV", x_title.Data(), bin_ev));
  const Double_t peak = h->GetMaximum();
  h->SetMaximum(1.9 * peak);
  h->SetMinimum(0);
  h->SetFillColor(kGray);
  h->SetFillStyle(1001);
  h->Draw("HIST");
  TLegend *leg = PlottingUtils::AddLegend(0.42, 0.92, 0.66, 0.92);
  leg->AddEntry(
      h, Form("%d replicas, %s", (Int_t)values.size(), kind_label.Data()), "f");
  const Int_t colors[3] = {kRed + 1, kBlue + 1, kGreen + 2};
  const Int_t styles[3] = {1, 2, 7};
  for (size_t k = 0; k < marks.size() && k < 3; k++) {
    TLine *l = new TLine(marks[k], 0, marks[k], 1.2 * peak);
    l->SetLineColor(colors[k]);
    l->SetLineStyle(styles[k]);
    l->SetLineWidth(3);
    l->Draw();
    leg->AddEntry(l, mark_labels[k], "l");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

void DrawPerSlice(const std::vector<NominalBudget> &budgets,
                  const std::vector<Int_t> &configs) {
  std::vector<TString> labels;
  for (size_t k = 0; k < configs.size(); k++)
    labels.push_back(CONFIG_TITLES[configs[k]]);
  Double_t legend_w = PlottingUtils::LegendWidthPx(labels);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)legend_w);
  plot->cd();
  TH1D *frame = new TH1D("per_slice_frame", "", N_BUDGET_SLICES, -0.5,
                         N_BUDGET_SLICES - 0.5);
  frame->SetDirectory(nullptr);
  frame->SetStats(0);
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
    frame->GetXaxis()->SetBinLabel(s + 1, SliceName(s));
  PlottingUtils::ConfigureHistogram(frame, kBlack,
                                    ";Rise slice [ns];E_{#gamma} [keV]");
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->SetMinimum(68.60);
  frame->SetMaximum(68.88);
  frame->Draw("AXIS");
  const Int_t colors[7] = {kGray + 2,  kBlack,      kRed + 1,  kBlue + 1,
                           kAzure + 8, kOrange + 7, kGreen + 2};
  const Int_t markers[7] = {24, 20, 21, 22, 26, 23, 33};
  std::vector<TGraphErrors *> graphs;
  for (size_t k = 0; k < configs.size(); k++) {
    TGraphErrors *g = new TGraphErrors();
    Double_t offset = 0.08 * (k - 0.5 * (configs.size() - 1.0));
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
      g->SetPoint(s, s + offset, budgets[k].value[s]);
      g->SetPointError(s, 0, budgets[k].error[s]);
    }
    PlottingUtils::ConfigureGraph(g, colors[k % 7]);
    g->SetMarkerStyle(markers[k % 7]);
    g->SetMarkerSize(1.4);
    g->Draw("P SAME");
    graphs.push_back(g);
  }
  TLine *ensdf = new TLine(-0.5, ReferenceLines::GE73M_ENSDF,
                           N_BUDGET_SLICES - 0.5, ReferenceLines::GE73M_ENSDF);
  ensdf->SetLineStyle(3);
  ensdf->SetLineColor(kGray + 1);
  ensdf->Draw();
  legend_pad->cd();
  TLegend *leg =
      new TLegend(0.05, 0.92, 0.98, 0.92 - 0.065 * (configs.size() + 1));
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetTextSize(30);
  leg->SetTextFont(43);
  for (size_t k = 0; k < graphs.size(); k++)
    leg->AddEntry(graphs[k], labels[k], "pe");
  leg->AddEntry(ensdf, "ENSDF 68.752 keV", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "ge_per_slice_configurations", "",
                            PlotSaveOptions::kLINEAR);
}

void SummarizeReplicas(Bool_t joint, const NominalBudget &reference,
                       const NominalBudget &observed, const NominalBudget &old,
                       const NominalBudget &now) {
  std::vector<NominalBudget> replicas = ReadReplicas(joint);
  if (replicas.size() < 2) {
    std::cout << "SUMMARY no " << (joint ? "joint" : "Gaussian") << " replicas"
              << std::endl;
    return;
  }
  TString kind = joint ? "joint" : "gaussian";
  std::vector<Double_t> shift, outlier, extrapolation,
      slice_shift[N_BUDGET_SLICES];
  for (size_t i = 0; i < replicas.size(); i++) {
    shift.push_back(1000 * (replicas[i].mean - reference.mean));
    outlier.push_back(1000 * replicas[i].outlier);
    extrapolation.push_back(1000 * replicas[i].extrapolation);
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
      slice_shift[s].push_back(1000 *
                               (replicas[i].value[s] - reference.value[s]));
  }
  Double_t mean = 0, rms = 0;
  MeanAndRms(shift, mean, rms);
  Double_t observed_shift = 1000 * (observed.mean - reference.mean);
  Double_t full_shift = 1000 * (now.mean - old.mean);
  std::cout << "SUMMARY " << kind << " replicas (" << replicas.size()
            << "): Ge shift from the old common records " << std::fixed
            << std::setprecision(1) << std::showpos << mean << std::noshowpos
            << " eV, rms " << rms << " eV; reprocessing the common records "
            << std::showpos << observed_shift << std::noshowpos << " eV ("
            << std::setprecision(2) << (observed_shift - mean) / rms
            << " rms, |shift| >= it in " << std::setprecision(0)
            << 100 * FractionAtLeast(shift, observed_shift)
            << "% of replicas), old to new in full " << std::showpos
            << std::setprecision(1) << full_shift << std::noshowpos << " eV ("
            << std::setprecision(2) << (full_shift - mean) / rms << " rms)"
            << std::endl;
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
    Double_t m_s = 0, r_s = 0;
    MeanAndRms(slice_shift[s], m_s, r_s);
    Double_t obs_s = 1000 * (observed.value[s] - reference.value[s]);
    std::cout << "SUMMARY " << kind << " slice " << std::setw(7) << SliceName(s)
              << " ns: replicas " << std::showpos << std::setprecision(1) << m_s
              << std::noshowpos << " rms " << r_s
              << " eV; common records new - old " << std::showpos << obs_s
              << std::noshowpos << " eV (" << std::setprecision(2)
              << (obs_s - m_s) / r_s << " rms)" << std::endl;
  }
  // Which slice the outlier rule drops as the lowest, and how much each
  // slice's curvature-free fit moves it.
  std::vector<Double_t> lowest(N_BUDGET_SLICES, 0);
  for (size_t i = 0; i < replicas.size(); i++)
    lowest[replicas[i].i_lo] += 100.0 / replicas.size();
  std::cout << "SUMMARY " << kind << " replicas' lowest slice:";
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
    std::cout << " " << SliceName(s) << " " << std::setprecision(0) << lowest[s]
              << "%";
  std::cout << std::endl;
  for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
    std::vector<Double_t> curvature_shift;
    for (size_t i = 0; i < replicas.size(); i++)
      curvature_shift.push_back(1000 *
                                (replicas[i].curved[s] - replicas[i].value[s]));
    Double_t c_mean = 0, c_rms = 0;
    MeanAndRms(curvature_shift, c_mean, c_rms);
    std::cout << "SUMMARY " << kind << " slice " << std::setw(7) << SliceName(s)
              << " ns: curvature-free minus nominal " << std::showpos
              << std::setprecision(1) << c_mean << std::noshowpos << " rms "
              << c_rms << " eV (old common " << std::showpos
              << 1000 * (reference.curved[s] - reference.value[s])
              << ", new common "
              << 1000 * (observed.curved[s] - observed.value[s])
              << std::noshowpos << ")" << std::endl;
  }
  Double_t o_mean = 0, o_rms = 0, x_mean = 0, x_rms = 0;
  MeanAndRms(outlier, o_mean, o_rms);
  MeanAndRms(extrapolation, x_mean, x_rms);
  Double_t o_frac = 0, x_frac = 0;
  for (size_t i = 0; i < outlier.size(); i++) {
    if (outlier[i] >= 1000 * now.outlier)
      o_frac += 1.0 / outlier.size();
    if (extrapolation[i] >= 1000 * now.extrapolation)
      x_frac += 1.0 / extrapolation.size();
  }
  std::cout << "SUMMARY " << kind << " replicas' outlier term: mean "
            << std::setprecision(1) << o_mean << " rms " << o_rms
            << " eV (old common " << 1000 * reference.outlier << ", new common "
            << 1000 * observed.outlier << ", old " << 1000 * old.outlier
            << ", new " << 1000 * now.outlier << " eV; >= new in "
            << std::setprecision(0) << 100 * o_frac << "%)" << std::endl;
  std::cout << "SUMMARY " << kind << " replicas' extrapolation term: mean "
            << std::setprecision(1) << x_mean << " rms " << x_rms
            << " eV (old common " << 1000 * reference.extrapolation
            << ", new common " << 1000 * observed.extrapolation << ", old "
            << 1000 * old.extrapolation << ", new " << 1000 * now.extrapolation
            << " eV; >= new in " << std::setprecision(0) << 100 * x_frac << "%)"
            << std::endl;

  TString title = joint ? "joint residual and slice" : "Gaussian residual";
  std::vector<Double_t> marks;
  std::vector<TString> mark_labels;
  marks.push_back(observed_shift);
  mark_labels.push_back("Common records, new #minus old");
  marks.push_back(full_shift);
  mark_labels.push_back("New #minus old result");
  DrawReplicaHistogram(shift, marks, mark_labels, "Replica #minus old Ge",
                       title, -30, 30, 2.0, "replica_ge_shift_" + kind);
  marks.clear();
  mark_labels.clear();
  marks.push_back(1000 * observed.outlier);
  mark_labels.push_back("Common records, new");
  marks.push_back(1000 * old.outlier);
  mark_labels.push_back("Old result");
  marks.push_back(1000 * now.outlier);
  mark_labels.push_back("New result");
  DrawReplicaHistogram(outlier, marks, mark_labels, "Outlier term", title, 0,
                       30, 1.0, "replica_outlier_" + kind);
  marks.clear();
  marks.push_back(1000 * observed.extrapolation);
  marks.push_back(1000 * old.extrapolation);
  marks.push_back(1000 * now.extrapolation);
  DrawReplicaHistogram(extrapolation, marks, mark_labels, "Extrapolation term",
                       title, 0, 45, 1.5, "replica_extrapolation_" + kind);
}

void Summary() {
  std::vector<NominalBudget> budgets(N_CONFIGS);
  for (Int_t c = 0; c < N_CONFIGS; c++) {
    std::ifstream in(BudgetFile(c).Data());
    TString label;
    if (!ReadBudget(in, label, budgets[c])) {
      std::cout << "SUMMARY missing " << BudgetFile(c) << std::endl;
      return;
    }
    PrintBudget(label, budgets[c]);
  }
  for (Int_t c = 0; c < N_CONFIGS; c++) {
    std::cout << "SUMMARY " << std::left << std::setw(16) << CONFIG_LABELS[c]
              << std::right << " per slice:";
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
      std::cout << " " << std::fixed << std::setprecision(4)
                << budgets[c].value[s] << "(" << std::setprecision(1)
                << 1000 * budgets[c].error[s] << ")";
    std::cout << "  vs old:";
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
      std::cout << " " << std::showpos << std::setprecision(1)
                << 1000 * (budgets[c].value[s] - budgets[C_OLD].value[s])
                << std::noshowpos;
    std::cout << "  mean " << std::showpos
              << 1000 * (budgets[c].mean - budgets[C_OLD].mean)
              << std::noshowpos << " eV" << std::endl;
  }
  // The outlier term's removal shifts: what dropping each slice does.
  for (Int_t c = 0; c < N_CONFIGS; c++) {
    std::cout << "SUMMARY " << std::left << std::setw(16) << CONFIG_LABELS[c]
              << std::right << " without slice:";
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++) {
      Double_t err = 0;
      Double_t m = WeightedMean(budgets[c].value, budgets[c].error,
                                N_BUDGET_SLICES, s, -1, err);
      std::cout << " " << SliceName(s) << " " << std::showpos
                << std::setprecision(1) << 1000 * (m - budgets[c].mean)
                << std::noshowpos;
    }
    std::cout << "  | curvature-free minus nominal per slice:";
    for (Int_t s = 0; s < N_BUDGET_SLICES; s++)
      std::cout << " " << std::showpos
                << 1000 * (budgets[c].curved[s] - budgets[c].value[s])
                << std::noshowpos << "(" << 1000 * budgets[c].curved_error[s]
                << ")";
    std::cout << std::endl;
  }
  std::vector<Int_t> shown;
  shown.push_back(C_OLD);
  shown.push_back(C_NEW);
  shown.push_back(C_COMMON_OLD);
  shown.push_back(C_COMMON_NEW);
  shown.push_back(C_NEW_E_OLD_RISE);
  shown.push_back(C_OLD_E_NEW_RISE);
  std::vector<NominalBudget> shown_budgets;
  for (size_t k = 0; k < shown.size(); k++)
    shown_budgets.push_back(budgets[shown[k]]);
  DrawPerSlice(shown_budgets, shown);
  SummarizeReplicas(kFALSE, budgets[C_COMMON_OLD], budgets[C_COMMON_NEW],
                    budgets[C_OLD_TIE_OLD], budgets[C_NEW]);
  SummarizeReplicas(kTRUE, budgets[C_COMMON_OLD], budgets[C_COMMON_NEW],
                    budgets[C_OLD_TIE_OLD], budgets[C_NEW]);
}

} // namespace CdTeReprocessing

void ReprocessingComparison(const char *stage = "summary", Int_t first = 0,
                            Int_t count = 1) {
  using namespace CdTeReprocessing;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG,
      project_root + "/plots/studies/psi_cdte/reprocessing_comparison",
      project_root + "/root_files", kFALSE);
  TGaxis::SetMaxDigits(3);
  const TString which = stage;
  if (which == "match")
    Match();
  else if (which == "events")
    Events();
  else if (which == "budget")
    Budget(first);
  else if (which == "replicas")
    Replicas(kFALSE, first, count);
  else if (which == "replicas_joint")
    Replicas(kTRUE, first, count);
  else if (which == "summary")
    Summary();
  else
    std::cout << "Unknown stage " << which << std::endl;
}
