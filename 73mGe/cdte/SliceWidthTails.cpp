// Whether the beam lines' low-side tailing comes from the gain still varying
// inside a rise-time slice. The Ge runs below 320 ns are cut into 5, 8, 12 and
// 16 slices of equal counts (EqualCountEdges, CdTeAnalysis.hpp); in each, the
// 95.902 keV capture line is fitted with its tails free, then the in-run fit
// runs with its tails free and with the capture line's (GeBudget's nominal).
// If the tailing came from the slices' spread in rise, the free low tail would
// shrink, and the two fits would converge, as the slices narrow. It does not:
// the tail switches on above about 282 ns however narrow the slice, so the
// plots show it against each slice's mean rise (its RMS is printed).
//
// Usage: root -b -q 'SliceWidthTails.cpp+' (the escape tie at GeBudget's
// nominal; pass another Te Ka/Cd Kb as the argument). Plots go to
// plots/cdte/slice_width.
#include "CdTeAnalysis.hpp"

#include <TGaxis.h>
#include <TGraph.h>
#include <TLegend.h>

namespace CdTeSliceWidth {

using namespace CdTe;

// The study spans the slices up to 320 ns, beyond the result's 281 ns.
const Double_t STUDY_MAX_RISE_NS = 320.0;
const Int_t N_CUTS = 4;
const Int_t SLICE_COUNTS[N_CUTS] = {5, 8, 12, 16};
const Int_t CUT_COLORS[N_CUTS] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
const Int_t CUT_MARKERS[N_CUTS] = {20, 21, 22, 23};

struct Event {
  Float_t e, rise_ns; // keV over the rise gain
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
    Double_t rise_ns = SAMPLE_NS * rise;
    if (height <= 0 || rise <= 0 || rise_ns >= STUDY_MAX_RISE_NS)
      continue;
    Event ev;
    ev.e = height * KEV_PER_UNIT / RiseGain(rise_ns);
    ev.rise_ns = rise_ns;
    if (ev.e < 120.0)
      out.push_back(ev);
  }
  return out;
}

// The slice's spectrum as FillRise makes it, and the mean and RMS of its rise
// times over 62-89 keV.
TH1D *FillSlice(const std::vector<Event> &events, Double_t lo, Double_t hi,
                const TString &name, Double_t &rise_mean, Double_t &rise_rms) {
  TH1D *h = new TH1D(name, Form(";Energy [keV];Counts / %.0f eV", 1000.0 * BIN),
                     (Int_t)(120.0 / BIN + 0.5), 0, 120);
  h->SetDirectory(nullptr);
  Double_t n = 0, sum = 0, sum2 = 0;
  for (size_t i = 0; i < events.size(); i++) {
    if (events[i].rise_ns < lo || events[i].rise_ns >= hi)
      continue;
    h->Fill(events[i].e);
    if (events[i].e >= 62.0 && events[i].e < 89.0) {
      n += 1;
      sum += events[i].rise_ns;
      sum2 += events[i].rise_ns * events[i].rise_ns;
    }
  }
  rise_mean = sum / n;
  rise_rms = TMath::Sqrt(TMath::Max(sum2 / n - rise_mean * rise_mean, 0.0));
  return h;
}

// The low tail's counts over the core's: amplitude x decay / sqrt(2 pi).
Double_t TailShare(Double_t amplitude, Double_t decay) {
  return amplitude * decay / TMath::Sqrt(2.0 * TMath::Pi());
}

void DrawAgainstRise(TGraph *graphs[N_CUTS], const TString &y_title,
                     const TString &name, Double_t y_lo, Double_t y_hi) {
  std::vector<TString> labels;
  for (Int_t c = 0; c < N_CUTS; c++)
    labels.push_back(Form("%d slices of equal counts", SLICE_COUNTS[c]));
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  for (Int_t c = 0; c < N_CUTS; c++) {
    PlottingUtils::ConfigureGraph(graphs[c], CUT_COLORS[c],
                                  ";Mean rise in the slice [ns];" + y_title);
    graphs[c]->SetMarkerStyle(CUT_MARKERS[c]);
    graphs[c]->SetMinimum(y_lo);
    graphs[c]->SetMaximum(y_hi);
    graphs[c]->Draw(c == 0 ? "AP" : "P SAME");
    leg->AddEntry(graphs[c], labels[c], "p");
  }
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

} // namespace CdTeSliceWidth

void SliceWidthTails(Double_t te_cdkb = CdTe::TE_CDKB) {
  using namespace CdTeSliceWidth;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/slice_width",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;
  g_kappa = KAPPA_SOURCE;
  ResetKbGroups();
  SetEscapeTie(KB_KA_NOMINAL, te_cdkb);

  std::vector<Event> events = LoadEvents(GeRuns());
  TH1D *rise = new TH1D("rise_all", "", 32000, 0, STUDY_MAX_RISE_NS);
  rise->SetDirectory(nullptr);
  for (size_t i = 0; i < events.size(); i++)
    if (events[i].e >= 62.0 && events[i].e < 89.0)
      rise->Fill(events[i].rise_ns);

  TGraph *free_share[N_CUTS], *capture_share[N_CUTS];
  TGraph *free_chi2[N_CUTS], *capture_chi2[N_CUTS];
  for (Int_t c = 0; c < N_CUTS; c++) {
    free_share[c] = new TGraph();
    capture_share[c] = new TGraph();
    free_chi2[c] = new TGraph();
    capture_chi2[c] = new TGraph();
    std::vector<Double_t> edges =
        EqualCountEdges(rise, SLICE_COUNTS[c], STUDY_MAX_RISE_NS);
    for (Int_t s = 0; s + 1 < (Int_t)edges.size(); s++) {
      Double_t mean = 0, rms = 0;
      TString tag = Form("n%d_s%d", SLICE_COUNTS[c], s);
      TH1D *h =
          FillSlice(events, edges[s], edges[s + 1], "h_" + tag, mean, rms);
      g_tails = TAILS_CAPTURE;
      PeakFit capture = FitCaptureLine(h, tag);
      g_capture_meas = capture.mu;

      g_tails = TAILS_FREE;
      TF1 *f_free = FitInRun(h, "free_" + tag);
      g_tails = TAILS_CAPTURE;
      TF1 *f_capture = FitInRun(h, "capture_" + tag);
      Double_t chi2_free = f_free->GetChisquare() / f_free->GetNDF();
      Double_t chi2_capture = f_capture->GetChisquare() / f_capture->GetNDF();

      std::cout << "WIDTH " << SLICE_COUNTS[c] << " slices, " << std::fixed
                << std::setprecision(2) << edges[s] << "-" << edges[s + 1]
                << " ns (rise mean " << mean << ", RMS " << rms
                << " ns): capture line low exp " << std::setprecision(3)
                << capture.tail[0] << " x " << capture.tail[1]
                << " sigma | free fit low exp "
                << f_free->GetParameter(Q_LOW_EXP) << " x "
                << f_free->GetParameter(Q_LOW_EXP_DECAY) << " sigma, sigma "
                << f_free->GetParameter(Q_SIG0) << ", chi2/ndf "
                << std::setprecision(2) << chi2_free << ", Ge "
                << std::setprecision(4) << f_free->GetParameter(Q_GE_E)
                << " +- " << std::setprecision(1)
                << 1000 * f_free->GetParError(Q_GE_E)
                << " eV | capture tails: chi2/ndf " << std::setprecision(2)
                << chi2_capture << ", Ge " << std::setprecision(4)
                << f_capture->GetParameter(Q_GE_E) << " +- "
                << std::setprecision(1) << 1000 * f_capture->GetParError(Q_GE_E)
                << " eV" << std::endl;
      free_share[c]->SetPoint(s, mean,
                              TailShare(f_free->GetParameter(Q_LOW_EXP),
                                        f_free->GetParameter(Q_LOW_EXP_DECAY)));
      capture_share[c]->SetPoint(s, mean,
                                 TailShare(capture.tail[0], capture.tail[1]));
      free_chi2[c]->SetPoint(s, mean, chi2_free);
      capture_chi2[c]->SetPoint(s, mean, chi2_capture);
      delete f_free;
      delete f_capture;
      delete h;
    }
  }
  DrawAgainstRise(free_share, "Low tail over core, in-run fit with free tails",
                  "slice_width_tail_share_free", 0.0, 1.8);
  DrawAgainstRise(capture_share, "Low tail over core, capture line",
                  "slice_width_tail_share_capture", 0.0, 1.8);
  DrawAgainstRise(free_chi2, "#chi^{2}/ndf, in-run fit with free tails",
                  "slice_width_chi2_free", 0.8, 3.0);
  DrawAgainstRise(capture_chi2,
                  "#chi^{2}/ndf, in-run fit with the capture line's tails",
                  "slice_width_chi2_capture", 0.8, 3.0);
}
