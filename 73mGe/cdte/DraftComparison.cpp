// Our CdTe spectra drawn the way the paper draft's CdTe-Spec-V2 figure is:
// the DT5780's own DPP energy on the linear Am + Ba source calibration,
// 64-76 keV, 20 eV bins, log y. Each count is spread uniformly over its
// channel so the 20 eV bins carry no channel comb.
//   - all Ge beam runs together;
//   - each run normalised to its own Pb Ka1 peak, to see which look like the
//     draft's figure;
//   - with and without a rise-time cut on our 10-90% rise.
//
// Usage: root -b -q 'DraftComparison.cpp+'
#include "CdTeAnalysis.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>

namespace CdTeDraft {

using namespace CdTe;

// The linear Am + Ba source calibration of the DPP energy, channels per keV
// and the offset in channels.
const Double_t DPP_GAIN = 22.413, DPP_OFFSET = 6.1;

Double_t DppEnergy(Double_t dpp, TRandom3 &rng) {
  return (dpp + rng.Uniform(-0.5, 0.5) - DPP_OFFSET) / DPP_GAIN;
}

void AllRuns() {
  std::vector<TString> runs = GeRuns();
  runs.push_back("PSI_Running_BeamTurnOn_20251009_1914_5");
  runs.push_back("PSI_Running_BeamTurnOn_20251009_1914_6");
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  Double_t dpp = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("digitizer_energy", 1);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  TH1D *h =
      new TH1D("like_draft", ";Energy [keV];Counts / 20 eV", 600, 64.0, 76.0);
  TRandom3 rng(7);
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    h->Fill(DppEnergy(dpp, rng));
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  PlottingUtils::ConfigureHistogram(h, kBlack);
  h->Draw("HIST");
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.62, 0.8, 0.9);
  leg->AddEntry(h, "CdTe Ge beam runs, DPP energy", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_cdte_dpp_like_draft", "",
                            PlotSaveOptions::kLOG);
}

void PerRun() {
  const std::vector<TString> runs = GeRuns();
  const char *labels[5] = {"10/09 BeamTurnOn_4", "10/10 GeTarget_1",
                           "10/10 GeTarget_2", "10/11 Ge_1", "10/11 Ge_2"};
  const Int_t colors[5] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2,
                           kOrange + 7};
  TRandom3 rng(11);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.55, 0.62, 0.9);
  for (Int_t r = 0; r < 5; r++) {
    TChain chain("features");
    chain.Add(FeaturesFile(runs[r]));
    Double_t dpp = 0;
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("digitizer_energy", 1);
    chain.SetBranchAddress("digitizer_energy", &dpp);
    TH1D *h = new TH1D(Form("run%d", r), "", 600, 64.0, 76.0);
    const Long64_t n_entries = chain.GetEntries();
    for (Long64_t i = 0; i < n_entries; i++) {
      chain.GetEntry(i);
      h->Fill(DppEnergy(dpp, rng));
    }
    h->Rebin(5);
    h->GetXaxis()->SetRangeUser(74.0, 75.5);
    Double_t peak = h->GetMaximum();
    h->GetXaxis()->SetRange(0, 0);
    h->Scale(1.0 / peak);
    h->SetTitle(";Energy [keV];Counts / 100 eV, Pb K#alpha_{1} peak = 1");
    PlottingUtils::ConfigureHistogram(h, colors[r]);
    h->SetMinimum(0.1);
    h->SetMaximum(2.0);
    h->Draw(r == 0 ? "HIST" : "HIST SAME");
    leg->AddEntry(
        h, Form("%s (%.0f counts at Pb K#alpha_{1})", labels[r], peak), "l");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_cdte_dpp_per_run_like_draft", "",
                            PlotSaveOptions::kLOG);
}

void RiseCut() {
  const std::vector<TString> runs = GeRuns();
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  Double_t dpp = 0;
  Float_t rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("digitizer_energy", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  chain.SetBranchAddress("rise_time", &rise);
  const Int_t n_cuts = 3;
  const Double_t cut_ns[n_cuts] = {1e9, 290.0, 280.0};
  const char *labels[n_cuts] = {"No rise-time cut", "Rise time < 290 ns",
                                "Rise time < 280 ns"};
  const Int_t colors[n_cuts] = {kBlack, kBlue + 1, kRed + 1};
  TH1D *h[n_cuts];
  for (Int_t k = 0; k < n_cuts; k++)
    h[k] = new TH1D(Form("cut%d", k), ";Energy [keV];Counts / 20 eV", 600, 64.0,
                    76.0);
  TRandom3 rng(5);
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    Double_t e = DppEnergy(dpp, rng);
    for (Int_t k = 0; k < n_cuts; k++)
      if (rise > 0 && SAMPLE_NS * rise < cut_ns[k])
        h[k]->Fill(e);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.55, 0.7, 0.9);
  for (Int_t k = 0; k < n_cuts; k++) {
    PlottingUtils::ConfigureHistogram(h[k], colors[k]);
    h[k]->SetMinimum(300);
    h[k]->SetMaximum(3.0e4);
    h[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
    leg->AddEntry(h[k],
                  Form("%s (%.0f%% of events)", labels[k],
                       100.0 * h[k]->GetEntries() / h[0]->GetEntries()),
                  "l");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_cdte_dpp_risecut_like_draft", "",
                            PlotSaveOptions::kLOG);
}

} // namespace CdTeDraft

void DraftComparison() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/draft_comparison",
                                project_root + "/root_files");
  CdTeDraft::AllRuns();
  CdTeDraft::PerRun();
  CdTeDraft::RiseCut();
}
