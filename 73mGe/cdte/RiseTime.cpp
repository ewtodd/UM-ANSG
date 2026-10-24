// The CdTe's rise-time (charge-loss) figures. Hole trapping lowers the pulse
// height of slow events; the beam-run lines are about three times broader
// than the source-run lines for that reason, which is why the Ge budget
// calibrates each rise-time slice on its own lines.
//   - Rise time (10-90% of the deconvolved step) against energy, source runs
//     against beam runs, on the Am-only scale.
//   - The relative pulse height against rise time from the Ba-133 81.00 keV
//     and Am-241 59.54 keV peaks of the source runs, per 10 ns slice relative
//     to the fast events, and the Ge beam runs' rise-energy plane before and
//     after dividing by it (the lines should stand up straight).
//
// Usage: root -b -q 'RiseTime.cpp+'
#include "CdTeAnalysis.hpp"

#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLegend.h>

namespace CdTeRiseTime {

using namespace CdTe;

struct Event {
  Float_t e, rise; // keV on the Am-only scale; samples
};

std::vector<Event> LoadEvents(const std::vector<TString> &runs) {
  TChain chain("trap");
  for (size_t i = 0; i < runs.size(); i++)
    chain.Add(TrapFile(runs[i]));
  Float_t energy[4], rise = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("energy", 1);
  chain.SetBranchStatus("rise", 1);
  chain.SetBranchAddress("energy", energy);
  chain.SetBranchAddress("rise", &rise);
  std::vector<Event> out;
  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    chain.GetEntry(i);
    if (energy[SHAPE] <= 0 || rise <= 0)
      continue;
    Event ev;
    ev.e = energy[SHAPE] * KEV_PER_UNIT;
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
    h->Fill(events[i].e, 10.0 * events[i].rise);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      Form(";Energy [keV];Rise time 10-90%% [ns];Counts / %.0f eV / 5 ns",
           1000.0 * bin_kev));
  PlottingUtils::SaveFigure(canvas, "rise_vs_energy_" + name, "",
                            PlotSaveOptions::kLINEAR);
}

// The slices of one sample (10 ns) the ridge is measured in, and the ones
// that define the fast events.
const Int_t RIDGE_RISE_LO = 25, RIDGE_RISE_HI = 40;
const Int_t FAST_LO = 26, FAST_HI = 29;

// A line's peak position in each rise slice, the highest bin between 84% and
// 103% of its energy refined by a Gaussian over +-0.4 keV, relative to the
// fast events'.
TGraphErrors *Ridge(const std::vector<Event> &events, Double_t e_line,
                    const TString &tag) {
  std::vector<TH1D *> slices;
  for (Int_t r = RIDGE_RISE_LO; r < RIDGE_RISE_HI; r++)
    slices.push_back(new TH1D(Form("ridge_%s_%d", tag.Data(), r), "", 300,
                              0.84 * e_line, 1.03 * e_line));
  for (size_t i = 0; i < events.size(); i++) {
    Int_t r = (Int_t)events[i].rise;
    if (r >= RIDGE_RISE_LO && r < RIDGE_RISE_HI)
      slices[r - RIDGE_RISE_LO]->Fill(events[i].e);
  }
  TGraphErrors *g = new TGraphErrors();
  for (Int_t r = RIDGE_RISE_LO; r < RIDGE_RISE_HI; r++) {
    TH1D *h = slices[r - RIDGE_RISE_LO];
    if (h->GetMaximum() < 40)
      continue;
    Double_t top = h->GetBinCenter(h->GetMaximumBin());
    TF1 f("ridge_gaus", "gaus", top - 0.4, top + 0.4);
    f.SetParameters(h->GetMaximum(), top, 0.3);
    h->Fit(&f, "QRN0L");
    g->SetPoint(g->GetN(), r + 0.5, f.GetParameter(1));
    g->SetPointError(g->GetN() - 1, 0, f.GetParError(1));
  }
  for (size_t i = 0; i < slices.size(); i++)
    delete slices[i];
  Double_t sw = 0, swx = 0;
  for (Int_t i = 0; i < g->GetN(); i++) {
    if (g->GetX()[i] < FAST_LO || g->GetX()[i] > FAST_HI + 1)
      continue;
    Double_t w = 1.0 / TMath::Power(g->GetEY()[i], 2);
    sw += w;
    swx += w * g->GetY()[i];
  }
  Double_t fast = swx / sw;
  for (Int_t i = 0; i < g->GetN(); i++) {
    g->SetPoint(i, 10.0 * g->GetX()[i], g->GetY()[i] / fast);
    g->SetPointError(i, 0, g->GetEY()[i] / fast);
    std::cout << "RIDGE " << tag << " rise " << g->GetX()[i]
              << " ns: relative height " << std::fixed << std::setprecision(4)
              << g->GetY()[i] << " +/- " << g->GetEY()[i] << std::endl;
  }
  return g;
}

// The relative height at a rise time (samples), flat outside the measured
// slices.
Double_t RelativeHeight(TGraphErrors *ridge, Double_t rise_samples) {
  Double_t x = 10.0 * rise_samples;
  if (x <= ridge->GetX()[0])
    return ridge->GetY()[0];
  if (x >= ridge->GetX()[ridge->GetN() - 1])
    return ridge->GetY()[ridge->GetN() - 1];
  return ridge->Eval(x);
}

void DrawRidges(TGraphErrors *ba, TGraphErrors *am) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(ba, kBlue + 1,
                                ";Rise time 10-90% [ns];Relative pulse height");
  PlottingUtils::ConfigureGraph(am, kRed + 1);
  ba->SetMarkerStyle(20);
  am->SetMarkerStyle(21);
  ba->Draw("AP");
  am->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.55, 0.9, 0.72, 0.9);
  leg->AddEntry(ba, "Ba-133 81.00 keV", "p");
  leg->AddEntry(am, "Am-241 59.54 keV", "p");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "rise_ridge_source_lines", "",
                            PlotSaveOptions::kLINEAR);
}

void DrawRiseEnergyPlane(const std::vector<Event> &events, TGraphErrors *ridge,
                         Bool_t corrected, const TString &name) {
  TH2D *h = new TH2D("h2_" + name, "", 400, 60, 100, 90, 150, 600);
  for (size_t i = 0; i < events.size(); i++) {
    Double_t e = corrected ? events[i].e / RelativeHeight(ridge, events[i].rise)
                           : events[i].e;
    h->Fill(e, 10.0 * events[i].rise);
  }
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      Form(";%s [keV];Rise time 10-90%% [ns];Counts / 100 eV / 5 ns",
           corrected ? "Corrected energy" : "Energy"));
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

} // namespace CdTeRiseTime

void RiseTime() {
  using namespace CdTeRiseTime;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG, project_root + "/plots/studies/psi_cdte/rise_time",
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

  // The four Ge beam runs from 10/10 and 10/11.
  const std::vector<TString> ge_runs = GeRuns();
  std::vector<TString> beam_runs(ge_runs.begin() + 1, ge_runs.end());
  std::vector<Event> ba = LoadEvents(BaSourceRuns());
  std::vector<Event> am = LoadEvents(AmSourceRuns());
  std::vector<Event> beam = LoadEvents(beam_runs);
  TGraphErrors *ridge_ba = Ridge(ba, 81.0, "Ba81");
  TGraphErrors *ridge_am = Ridge(am, 59.54, "Am59");
  DrawRidges(ridge_ba, ridge_am);
  DrawRiseEnergyPlane(beam, ridge_ba, kFALSE, "rise_vs_energy_beam_ge_zoom");
  DrawRiseEnergyPlane(beam, ridge_ba, kTRUE,
                      "rise_vs_corrected_energy_beam_ge");
}
