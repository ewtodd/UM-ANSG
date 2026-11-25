// How the CdTe Ge result depends on the upper rise-time cut. For each cut, the
// Ge runs below it are split into four slices of equal counts
// (EqualCountEdges), each fitted in-run as GeBudget.cpp's nominal (tails
// free) and with the 95.902 keV capture line's tails; the slices' weighted
// mean is printed per cut. A cut that only removes events Ge and the
// calibration lines share would leave the mean where it is. Ge is read on
// the in-run scale: the source-run calibration holds the source runs' shape
// of each slice, which is measured only for the result's slices.
//
// Usage: root -b -q 'RiseCutScan.cpp+'. Plots go to plots/cdte/rise_cut.
#include "CdTeAnalysis.hpp"

#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>

namespace CdTeRiseCut {

using namespace CdTe;

const Int_t N_CUTS = 5;
const Double_t CUTS_NS[N_CUTS] = {281, 285, 290, 300, 320};
// GeBudget.cpp's nominal escape tie.
const Double_t KB_KA = 0.275, TE_CDKB = 0.352;
const Int_t N_TAILS = 2;
const Int_t TAIL_MODES[N_TAILS] = {TAILS_FREE, TAILS_CAPTURE};
const char *const TAIL_NAMES[N_TAILS] = {"tails free", "capture line's tails"};

} // namespace CdTeRiseCut

void RiseCutScan() {
  using namespace CdTeRiseCut;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/rise_cut",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;
  g_kappa = KAPPA_SOURCE;
  ResetKbGroups();
  SetEscapeTie(KB_KA, TE_CDKB);

  TH1D *rise = FillRiseDistribution(GeRuns(), "rise_ge", 62.0, 89.0);
  TGraphErrors *means[N_TAILS];
  for (Int_t m = 0; m < N_TAILS; m++)
    means[m] = new TGraphErrors();
  for (Int_t c = 0; c < N_CUTS; c++) {
    std::vector<Double_t> edges =
        EqualCountEdges(rise, N_RISE_SLICES, CUTS_NS[c]);
    Double_t sw[N_TAILS] = {0, 0}, swx[N_TAILS] = {0, 0};
    for (Int_t s = 0; s < N_RISE_SLICES; s++) {
      TH1D *h =
          FillRise(GeRuns(), Form("h_%d_%d", c, s), edges[s], edges[s + 1]);
      g_tails = TAILS_FREE;
      PeakFit capture = FitCaptureLine(h, Form("capture_%d_%d", c, s));
      g_capture_meas = capture.mu;
      for (Int_t m = 0; m < N_TAILS; m++) {
        g_tails = TAIL_MODES[m];
        TF1 *f = FitInRun(h, Form("fit_%d_%d_%d", c, s, m));
        Double_t ge = f->GetParameter(Q_GE_E), err = f->GetParError(Q_GE_E);
        sw[m] += 1.0 / (err * err);
        swx[m] += ge / (err * err);
        std::cout << "RISE-CUT " << std::fixed << std::setprecision(0)
                  << CUTS_NS[c] << " ns, slice " << std::setprecision(2)
                  << edges[s] << "-" << edges[s + 1] << " ns, " << TAIL_NAMES[m]
                  << ": Ge " << std::setprecision(4) << ge << " +- "
                  << std::setprecision(1) << 1000 * err << " eV, chi2/ndf "
                  << std::setprecision(2) << f->GetChisquare() / f->GetNDF()
                  << ", low exp " << std::setprecision(3)
                  << f->GetParameter(Q_LOW_EXP) << " x "
                  << f->GetParameter(Q_LOW_EXP_DECAY) << " sigma" << std::endl;
        delete f;
      }
      delete h;
    }
    for (Int_t m = 0; m < N_TAILS; m++) {
      Double_t mean = swx[m] / sw[m], err = 1.0 / TMath::Sqrt(sw[m]);
      std::cout << "RISE-CUT " << std::fixed << std::setprecision(0)
                << CUTS_NS[c] << " ns, " << TAIL_NAMES[m] << ": weighted mean "
                << std::setprecision(4) << mean << " +- "
                << std::setprecision(1) << 1000 * err << " eV" << std::endl;
      means[m]->SetPoint(c, CUTS_NS[c] + (m == 0 ? -0.5 : 0.5), mean);
      means[m]->SetPointError(c, 0, err);
    }
  }

  std::vector<TString> labels;
  for (Int_t m = 0; m < N_TAILS; m++)
    labels.push_back(TAIL_NAMES[m]);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  const Int_t colors[N_TAILS] = {kBlue + 1, kRed + 1};
  for (Int_t m = 0; m < N_TAILS; m++) {
    PlottingUtils::ConfigureGraph(
        means[m], colors[m],
        ";Upper rise-time cut [ns];Ge, four equal-count slices, in-run scale [keV]");
    means[m]->SetMarkerStyle(20 + m);
    means[m]->Draw(m == 0 ? "AP" : "P SAME");
    leg->AddEntry(means[m], labels[m], "pe");
  }
  PlottingUtils::SaveFigure(canvas, "rise_cut_scan", "",
                            PlotSaveOptions::kLINEAR);
}
