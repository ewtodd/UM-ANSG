// The CdTe's rise-time figures: why the spectra are sliced in rise and why
// each height is divided by the rise gain (CdTeAnalysis.hpp).
//   - Rise time (10-90% of the deconvolved step) against energy, source runs
//     against beam runs, as measured. Hole trapping lowers the pulse height
//     of slow events; the beam-run lines are about three times broader than
//     the source-run lines for that reason.
//   - Each source line's core position over its true energy against rise,
//     as measured and over the rise gain, in bins of equal precision; and
//     each line's rise distribution. The lines sample different rise ranges,
//     so inside a slice they would sample different gains.
//   - The Ge beam runs' rise-energy plane as measured and over the rise gain,
//     the latter with the result's slice edges; and those edges recomputed
//     as equal-count quantiles (EqualCountEdges) next to CdTeAnalysis.hpp's.
//   - The 77.4 keV check: the slices the result uses (up to 320 ns), as
//     measured against over the rise gain, for the Ge runs and the no-Ge
//     runs. An earlier correction, applied at every rise time, put a line
//     near 77.4 keV into the beam spectra.
//
// Usage: root -b -q 'RiseTime.cpp+'
#include "CdTeAnalysis.hpp"

#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>

namespace CdTeRiseTime {

using namespace CdTe;

struct Event {
  Float_t e, rise; // keV on the Am-only scale, as measured; samples
};

std::vector<Event> LoadEvents(const std::vector<TString> &runs) {
  TChain chain("features");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(FeaturesFile(runs[i]));
  Float_t height = 0, rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("trapezoid_height", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("rise_time", &rise);
  std::vector<Event> out;
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    if (height <= 0 || rise <= 0)
      continue;
    Event ev;
    ev.e = height * KEV_PER_UNIT;
    ev.rise = rise;
    out.push_back(ev);
  }
  return out;
}

void PlotRiseVsEnergy(const std::vector<TString> &runs, const TString &name) {
  std::vector<Event> events = LoadEvents(runs);
  const Double_t e_lo = 20, e_hi = 110, bin_kev = 0.1;
  Int_t n_e = (Int_t)((e_hi - e_lo) / bin_kev + 0.5);
  TH2D *h = new TH2D("rise_" + name, "", n_e, e_lo, e_hi, 90, 150, 600);
  for (size_t i = 0; i < events.size(); i++)
    h->Fill(events[i].e, SAMPLE_NS * events[i].rise);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      Form(";Energy [keV];Rise time 10-90%% [ns];Counts / %.0f eV / 5 ns",
           1000.0 * bin_kev));
  PlottingUtils::SaveFigure(canvas, "rise_vs_energy_" + name, "",
                            PlotSaveOptions::kLINEAR);
}

// The core of the Ba 81.00 keV line's Cd K-alpha escape doublet, by its
// sub-lines' shares.
Double_t EscapeCoreEnergy() {
  Double_t share = 0, weighted = 0;
  for (Int_t i = 0; i < ReferenceLines::N_ESCAPE_LINES; i++) {
    if (ReferenceLines::ESCAPE_SERIES[i] != 0)
      continue;
    share += GEANT4_SERIES_SHARE[i];
    weighted += GEANT4_SERIES_SHARE[i] * ReferenceLines::ESCAPE_XRAY_KEV[i];
  }
  return ReferenceLines::BA133_81 - weighted / share;
}

const Int_t N_LINES = 3;
const char *const LINE_LABEL[N_LINES] = {"Am 59.54 keV (Am runs)",
                                         "Ba 81.00 keV (Ba runs)",
                                         "its Cd K#alpha escape (Ba runs)"};
const Int_t LINE_COLOR[N_LINES] = {kRed + 1, kBlue + 1, kGreen + 2};

// Each line's core position over its true energy against rise, the heights
// as measured or over the rise gain; each point's x error spans its bin.
void PlotGainByLine(TH2D *am, TH2D *ba, const TString &name,
                    const TString &header, Double_t y_lo, Double_t y_hi) {
  const Double_t e[N_LINES] = {ReferenceLines::AM241, ReferenceLines::BA133_81,
                               EscapeCoreEnergy()};
  TH2D *source[N_LINES] = {am, ba, ba};
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.62, 0.18, 0.42);
  leg->SetHeader(header);
  for (Int_t l = 0; l < N_LINES; l++) {
    TGraphAsymmErrors *g = CoreRatioByRise(source[l], e[l], 0.02, 0.04);
    PlottingUtils::ConfigureGraph(g, LINE_COLOR[l], "");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    if (l == 0) {
      g->SetTitle(";10-90% rise [ns];Core position / true energy");
      g->Draw("AP");
      g->GetXaxis()->SetLimits(255, 345);
      g->GetYaxis()->SetRangeUser(y_lo, y_hi);
    } else {
      g->Draw("P SAME");
    }
    leg->AddEntry(g, LINE_LABEL[l], "pe");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

// Each line's rise distribution, its core events (the heights over the rise
// gain, so the core sits at the same place at every rise).
void PlotRiseDensity(TH2D *am, TH2D *ba) {
  const Double_t e[N_LINES] = {ReferenceLines::AM241, ReferenceLines::BA133_81,
                               EscapeCoreEnergy()};
  TH2D *source[N_LINES] = {am, ba, ba};
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TLegend *leg = PlottingUtils::AddLegend(0.5, 0.92, 0.66, 0.9);
  for (Int_t l = 0; l < N_LINES; l++) {
    TH1D *d = source[l]->ProjectionX(
        Form("rise_density_%d", l), source[l]->GetYaxis()->FindBin(e[l] - 0.6),
        source[l]->GetYaxis()->FindBin(e[l] + 0.4));
    d->SetDirectory(nullptr);
    d->Rebin(4);
    d->Scale(1.0 / d->Integral(d->FindBin(250.5), d->FindBin(329.5)));
    d->SetStats(0);
    d->SetLineColor(LINE_COLOR[l]);
    d->SetLineWidth(2);
    d->SetTitle(";10-90% rise [ns];Fraction / 1 ns");
    d->GetXaxis()->SetRangeUser(250, 330);
    d->Draw(l == 0 ? "hist" : "hist same");
    leg->AddEntry(d, LINE_LABEL[l], "l");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "rise_density_by_line", "",
                            PlotSaveOptions::kLINEAR);
}

void DrawRiseEnergyPlane(const std::vector<Event> &events, Bool_t corrected,
                         const TString &name) {
  TH2D *h = new TH2D("h2_" + name, "", 400, 60, 100, 90, 150, 600);
  for (size_t i = 0; i < events.size(); i++) {
    Double_t rise_ns = SAMPLE_NS * events[i].rise;
    Double_t e = corrected ? events[i].e / RiseGain(rise_ns) : events[i].e;
    h->Fill(e, rise_ns);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      Form(";%s [keV];Rise time 10-90%% [ns];Counts / 100 eV / 5 ns",
           corrected ? "Energy over the rise gain" : "Energy"));
  if (corrected)
    for (Int_t s = 1; s <= N_RISE_SLICES; s++) {
      TLine *edge = new TLine(60, RISE_SLICE_LO[s], 100, RISE_SLICE_LO[s]);
      edge->SetLineColor(kBlack);
      edge->SetLineStyle(2);
      edge->Draw();
    }
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

// The slices the result uses, as measured against over the rise gain, with
// 77.4 keV marked.
void PlotGhostCheck(const std::vector<TString> &runs, const TString &name,
                    const TString &label) {
  g_rise_gain_correction = kFALSE;
  TH1D *measured = FillRise(runs, name + "_measured", 0, RESULT_MAX_RISE_NS);
  g_rise_gain_correction = kTRUE;
  TH1D *corrected = FillRise(runs, name + "_corrected", 0, RESULT_MAX_RISE_NS);
  measured->Rebin(2);
  corrected->Rebin(2);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  canvas->SetLogy();
  corrected->SetTitle(";Energy [keV];Counts / 100 eV");
  corrected->GetXaxis()->SetRangeUser(62, 89);
  corrected->SetStats(0);
  corrected->SetLineColor(kBlue + 1);
  corrected->SetLineWidth(2);
  measured->SetLineColor(kGray + 2);
  measured->SetLineWidth(2);
  Double_t y_lo = 0.8 * corrected->GetMinimum(1);
  Double_t y_hi = 6.0 * corrected->GetMaximum();
  corrected->SetMinimum(y_lo);
  corrected->SetMaximum(y_hi);
  corrected->Draw("hist");
  measured->Draw("hist same");
  TLine *mark = new TLine(77.4, y_lo, 77.4, y_hi);
  mark->SetLineColor(kRed + 1);
  mark->SetLineStyle(2);
  mark->Draw();
  TLegend *leg = PlottingUtils::AddLegend(0.5, 0.92, 0.68, 0.9);
  leg->SetHeader(label + Form(", rise < %.0f ns", RESULT_MAX_RISE_NS));
  leg->AddEntry(corrected, "Over the rise gain", "l");
  leg->AddEntry(measured, "As measured", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "ghost_check_" + name, "",
                            PlotSaveOptions::kLOG);
  delete measured;
  delete corrected;
}

} // namespace CdTeRiseTime

void RiseTime() {
  using namespace CdTeRiseTime;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/rise_gain",
                                project_root + "/root_files");
  TGaxis::SetMaxDigits(3);

  std::vector<TString> am_two_days;
  am_two_days.push_back("PSI_Calibration_CdTe_Am_Si_Am_NoGe_20251010_1014_2");
  am_two_days.push_back("PSI_Calibration_CdTe_Am_Si_Am_20251013_0652_1");
  std::vector<TString> ge_1010;
  ge_1010.push_back("PSI_Beam_GeTarget_20251010_1525_2");
  std::vector<TString> noge_1010;
  noge_1010.push_back("PSI_Beam_NOGe_20251010_1106_3");
  PlotRiseVsEnergy(am_two_days, "source_am");
  PlotRiseVsEnergy(ge_1010, "beam_ge_1010");
  PlotRiseVsEnergy(noge_1010, "beam_noge_1010");

  TH2D *am_measured = FillRiseHeight(AmSourceRuns(), "am_measured", kFALSE);
  TH2D *ba_measured = FillRiseHeight(BaSourceRuns(), "ba_measured", kFALSE);
  TH2D *am_corrected = FillRiseHeight(AmSourceRuns(), "am_corrected", kTRUE);
  TH2D *ba_corrected = FillRiseHeight(BaSourceRuns(), "ba_corrected", kTRUE);
  PlotGainByLine(am_measured, ba_measured, "rise_gain_by_line",
                 "Heights as measured", 0.972, 1.006);
  PlotGainByLine(am_corrected, ba_corrected, "rise_gain_by_line_corrected",
                 "Heights over the Ba 81.00 keV rise gain", 0.994, 1.006);
  PlotRiseDensity(am_corrected, ba_corrected);

  TH1D *rise_ge = FillRiseDistribution(GeRuns(), "rise_ge", 62.0, 89.0);
  std::vector<Double_t> edges =
      EqualCountEdges(rise_ge, N_RISE_SLICES, RESULT_MAX_RISE_NS);
  std::cout << "RISE SLICES of equal counts below " << std::fixed
            << std::setprecision(0) << RESULT_MAX_RISE_NS
            << " ns (Ge runs, 62-89 keV):";
  for (size_t k = 0; k < edges.size(); k++)
    std::cout << " " << std::setprecision(2) << edges[k];
  std::cout << "; CdTeAnalysis.hpp has";
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    std::cout << " " << RISE_SLICE_LO[s];
  std::cout << " " << RESULT_MAX_RISE_NS << std::endl;
  for (Int_t s = 0; s < N_RISE_SLICES; s++)
    std::cout << "RISE SLICE " << std::setprecision(2) << RISE_SLICE_LO[s]
              << "-" << RISE_SLICE_HI[s] << " ns: " << std::setprecision(0)
              << rise_ge->Integral(rise_ge->FindBin(RISE_SLICE_LO[s] + 1e-6),
                                   rise_ge->FindBin(RISE_SLICE_HI[s] - 1e-6))
              << " events" << std::endl;
  delete rise_ge;

  std::vector<Event> beam = LoadEvents(GeRuns());
  DrawRiseEnergyPlane(beam, kFALSE, "rise_vs_energy_beam_ge");
  DrawRiseEnergyPlane(beam, kTRUE, "rise_vs_corrected_energy_beam_ge");

  PlotGhostCheck(GeRuns(), "ge_runs", "Ge beam runs");
  PlotGhostCheck(NoGePbRuns(), "noge_runs", "No-Ge beam runs");
}
