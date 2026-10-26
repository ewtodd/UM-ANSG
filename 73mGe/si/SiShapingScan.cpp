// The trapezoid shaping for the PSI Si waveforms, chosen by resolution: the
// Am-241 59.54 keV line in two source runs and Pb K-alpha in two beam runs
// (one of each with 2000- and one with 4500-sample records), for a grid of
// rise and flat-top times and two baselines (the record's first 350 samples,
// or every sample up to 75 before the nominal trigger), each processed
// through Analysis-Utilities' PulseHeightProcessor. The widths are fitted as
// in the Ge analysis (SiAnalysis.hpp), in digitiser channels.
//
// Usage: root -b -q 'SiShapingScan.cpp+(<rise>, <flat top>, <near baseline>)'
//        processes one setting; root -b -q 'SiShapingScan.cpp+' fits and
//        draws every setting processed.
#include "SiAnalysis.hpp"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>

namespace SiScan {

using namespace Si;

const Int_t N_RISE = 6, N_FLAT = 3;
const Int_t RISE_SAMPLES[N_RISE] = {100, 200, 300, 400, 500, 600};
const Int_t FLAT_SAMPLES[N_FLAT] = {50, 100, 150};
const char *const SCAN_DIR = "/labdata/ANSG/PSI/processed/Si/shaping_scan";
const Long64_t BEAM_EVENTS = 1500000;
const Int_t FAR_BASELINE_SAMPLES = 350;

struct ScanRun {
  TString tag, run;
  Bool_t beam;
};

std::vector<ScanRun> ScanRuns() {
  std::vector<ScanRun> runs;
  runs.push_back(
      {"am10", "PSI_Calibration_CdTe_Am_Si_Am_20251010_1422_1", kFALSE});
  runs.push_back(
      {"am13", "PSI_Calibration_CdTe_Am_Si_Am_20251013_0652_1", kFALSE});
  runs.push_back({"beam11", "PSI_Beam_Ge_20251011_1715_2", kTRUE});
  runs.push_back(
      {"beam12", "PSI_Beam_LongWf_CdTe_NoGe_Si_YesGe_20251012_1517_1", kTRUE});
  return runs;
}

TString SettingTag(Int_t rise, Int_t flat, Bool_t near_baseline) {
  return Form("k%d_m%d_%s", rise, flat, near_baseline ? "near" : "far");
}

TString ScanFile(const ScanRun &run, const TString &setting) {
  return TString(SCAN_DIR) + "/" + run.tag + "_" + setting + ".root";
}

void Process(Int_t rise, Int_t flat, Bool_t near_baseline) {
  PulseHeightSettings filter = Filter();
  filter.rise_samples = rise;
  filter.flat_top_samples = flat;
  TString setting = SettingTag(rise, flat, near_baseline);
  std::vector<ScanRun> runs = ScanRuns();
  for (size_t i = 0; i < runs.size(); i++) {
    FileProcessingConfig config = Config(RecordLength(runs[i].run));
    // The trapezoid at its peaking point reaches rise + flat top / 10
    // samples before the trigger.
    config.pre_samples = TMath::Max(config.pre_samples, rise + flat / 10 + 10);
    if (!near_baseline)
      config.num_samples_baseline = FAR_BASELINE_SAMPLES;
    config.verbose = kFALSE;
    if (runs[i].beam)
      config.max_events = BEAM_EVENTS;
    PulseHeightProcessor processor(config, filter);
    processor.ProcessFile(RawFile(runs[i].run), runs[i].tag + "_" + setting);
  }
}

// The trapezoid spectrum of one scan file in digitiser channels.
TH1D *ScanSpectrum(const TString &file, const TString &name) {
  TH1D *h = NewSpectrum(name);
  TFile in(file);
  TTree *tree = static_cast<TTree *>(in.Get("features"));
  if (!tree)
    return h;
  Float_t trapezoid_height = 0, baseline_rms = 0;
  tree->SetBranchAddress("trapezoid_height", &trapezoid_height);
  tree->SetBranchAddress("baseline_rms", &baseline_rms);
  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (baseline_rms < FLOOR_BASELINE_RMS)
      continue;
    h->Fill(trapezoid_height / TRAPEZOID_PER_CHANNEL);
  }
  return h;
}

// The line width in eV: Am alone, or Pb K-alpha2 and K-alpha1 with one width,
// the windows the Ge analysis uses.
Double_t LineWidth(TH1D *h, Bool_t beam, Double_t &width_err) {
  std::vector<LineSpec> lines;
  FitOutput fit;
  if (beam) {
    lines.resize(2);
    lines[0].mean_guess = 1745;
    lines[1].mean_guess = 1797;
    fit = FitLines(h, 1710, 1840, lines, 1, 1, 7.0);
  } else {
    lines.resize(1);
    lines[0].mean_guess = 1428;
    fit = FitLines(h, 1395, 1465, lines, 1, 1);
  }
  width_err = 1000.0 * fit.sigma_err[0] / DPP_PER_KEV;
  return 1000.0 * fit.sigma[0] / DPP_PER_KEV;
}

void Summarise() {
  std::vector<ScanRun> runs = ScanRuns();
  Int_t n_runs = (Int_t)runs.size();
  std::vector<TGraphErrors *> graphs[2][N_FLAT];
  for (Int_t b = 0; b < 2; b++)
    for (Int_t f = 0; f < N_FLAT; f++)
      for (Int_t r = 0; r < n_runs; r++)
        graphs[b][f].push_back(new TGraphErrors());
  for (Int_t b = 0; b < 2; b++) {
    Bool_t near_baseline = b == 1;
    for (Int_t f = 0; f < N_FLAT; f++) {
      for (Int_t k = 0; k < N_RISE; k++) {
        TString setting =
            SettingTag(RISE_SAMPLES[k], FLAT_SAMPLES[f], near_baseline);
        for (Int_t r = 0; r < n_runs; r++) {
          TString file = ScanFile(runs[r], setting);
          if (gSystem->AccessPathName(file))
            continue;
          TH1D *h = ScanSpectrum(file, "scan_" + setting);
          Double_t err = 0;
          Double_t width = LineWidth(h, runs[r].beam, err);
          TGraphErrors *g = graphs[b][f][r];
          Int_t point = g->GetN();
          g->SetPoint(point, SAMPLE_NS * RISE_SAMPLES[k] / 1000.0, width);
          g->SetPointError(point, 0, err);
          std::cout << "SI-SHAPING " << std::left << std::setw(20) << setting
                    << std::right << " " << std::setw(7) << runs[r].tag
                    << (runs[r].beam ? " Pb Ka" : " Am 59") << " sigma "
                    << std::fixed << std::setprecision(1) << width << " +- "
                    << err << " eV  entries " << std::setprecision(0)
                    << h->GetEntries() << std::endl;
          delete h;
        }
      }
    }
  }
  const Int_t colors[N_FLAT] = {kAzure + 2, kRed + 1, kGreen + 2};
  for (Int_t r = 0; r < n_runs; r++) {
    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
    TLegend *legend = PlottingUtils::AddLegend(0.45, 0.93, 0.62, 0.92);
    Double_t y_lo = 1e9, y_hi = -1e9;
    for (Int_t b = 0; b < 2; b++) {
      for (Int_t f = 0; f < N_FLAT; f++) {
        TGraphErrors *g = graphs[b][f][r];
        for (Int_t i = 0; i < g->GetN(); i++) {
          if (g->GetX()[i] < 1.5)
            continue;
          y_lo = TMath::Min(y_lo, g->GetY()[i] - g->GetEY()[i]);
          y_hi = TMath::Max(y_hi, g->GetY()[i] + g->GetEY()[i]);
        }
      }
    }
    y_lo -= 10;
    y_hi += 0.8 * (y_hi - y_lo);
    TH1D *frame = new TH1D(Form("frame_%d", r), "", 10, 1.5, 6.5);
    frame->SetDirectory(nullptr);
    PlottingUtils::ConfigureHistogram(
        frame, kBlack,
        Form(";Trapezoid rise [#mus];#sigma, %s [eV]",
             runs[r].beam ? "Pb K#alpha" : "^{241}Am 59.54 keV"));
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->SetStats(0);
    frame->Draw("AXIS");
    for (Int_t b = 0; b < 2; b++) {
      for (Int_t f = 0; f < N_FLAT; f++) {
        TGraphErrors *g = graphs[b][f][r];
        if (g->GetN() == 0)
          continue;
        PlottingUtils::ConfigureGraph(g, colors[f]);
        g->SetMarkerStyle(b == 1 ? 20 : 24);
        g->SetLineStyle(b == 1 ? 1 : 2);
        g->Draw("PL SAME");
        legend->AddEntry(g,
                         Form("flat top %.1f #mus, %s baseline",
                              SAMPLE_NS * FLAT_SAMPLES[f] / 1000.0,
                              b == 1 ? "near" : "far"),
                         "lp");
      }
    }
    legend->Draw();
    PlottingUtils::SaveFigure(canvas, "shaping_" + runs[r].tag, "",
                              PlotSaveOptions::kLINEAR);
  }
}

} // namespace SiScan

void SiShapingScan(Int_t rise = 0, Int_t flat = 0, Int_t near_baseline = 1) {
  using namespace SiScan;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/si/shaping", SCAN_DIR);
  if (rise > 0)
    Process(rise, flat, near_baseline != 0);
  else
    Summarise();
}
