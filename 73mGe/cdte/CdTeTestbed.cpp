// A testbed for the CdTe Ge fit: whether the nominal in-run fit
// (CdTeAnalysis.hpp, GeBudget.cpp's nominal) reads Ge without bias. In each
// cell, a day's Ge runs in one rise slice, the data are fitted under a
// "truth" model; Poisson pseudo-spectra are drawn from that fit over the fit
// window, and each is fitted with the nominal model. Ge is read on the day's
// source-run scale (GainTransfer), the source lines held at their fits to
// the data; the bias is the fitted Ge less the truth's.
// Truths:
//   nominal        the nominal fit itself (closure: the fit's own bias);
//   step free      the step free, as the data prefer (0.02-0.05);
//   capture tails  both tails at the 95.902 keV capture line's;
//   source tail    the low tail at the source runs';
//   kappa free     the scale's curvature free (the beam fit's own).
// The capture line's position stays at its fit to the data, and the escape
// tie at GeBudget.cpp's nominal.
//
// Usage: root -b -q 'CdTeTestbed.cpp+(<truth 0-4>, <trials>)'. Plots go to
// plots/testbed/cdte.
#include "CdTeAnalysis.hpp"

#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom3.h>

namespace CdTeTestbedFit {

using namespace CdTe;

enum {
  T_NOMINAL,
  T_STEP_FREE,
  T_CAPTURE_TAILS,
  T_SOURCE_TAIL,
  T_KAPPA_FREE,
  N_TRUTHS
};
const char *const TRUTH_NAMES[N_TRUTHS] = {
    "nominal", "step free", "capture tails", "source tail", "kappa free"};
const char *const TRUTH_TAGS[N_TRUTHS] = {
    "nominal", "step_free", "capture_tails", "source_tail", "kappa_free"};

void SetNominal() {
  g_tails = TAILS_FREE;
  g_free_step = kFALSE;
  g_free_kappa = kFALSE;
  g_kappa = KAPPA_SOURCE;
  g_ge_own_width = kFALSE;
  g_ge_own_shape = kFALSE;
  g_flat_continuum = kFALSE;
  g_slope_shift = 0;
}

void SetTruth(Int_t truth) {
  SetNominal();
  if (truth == T_STEP_FREE)
    g_free_step = kTRUE;
  if (truth == T_CAPTURE_TAILS)
    g_tails = TAILS_CAPTURE;
  if (truth == T_SOURCE_TAIL)
    g_tails = TAILS_SOURCE;
  if (truth == T_KAPPA_FREE)
    g_free_kappa = kTRUE;
}

// A pseudo-spectrum: the data outside the fit window, Poisson draws of the
// truth's expectation inside it.
TH1D *Pseudo(TH1D *data, TF1 *truth, TRandom3 &random, const TString &name) {
  TH1D *h = (TH1D *)data->Clone(name);
  h->SetDirectory(nullptr);
  for (Int_t b = h->FindBin(g_lo); b <= h->FindBin(g_hi); b++) {
    Double_t expected = TMath::Max(truth->Eval(h->GetBinCenter(b)), 0.0);
    Double_t counts = random.Poisson(expected);
    h->SetBinContent(b, counts);
    h->SetBinError(b, TMath::Sqrt(TMath::Max(counts, 1.0)));
  }
  return h;
}

} // namespace CdTeTestbedFit

void CdTeTestbed(Int_t truth = 0, Int_t trials = 40) {
  using namespace CdTeTestbedFit;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/testbed/cdte",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;
  ResetKbGroups();
  SetEscapeTie(KB_KA_NOMINAL, TE_CDKB);
  TRandom3 random(4242 + truth);
  SetNominal();
  SourceLine lines[N_GE_DAYS][N_RISE_SLICES];
  FitSourceLines(lines);

  TGraphErrors *bias[N_GE_DAYS];
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    bias[d] = new TGraphErrors();
  Double_t sw = 0, swb = 0, swb2 = 0;
  for (Int_t d = 0; d < N_GE_DAYS; d++) {
    for (Int_t s = 0; s < N_RISE_SLICES; s++) {
      const SourceLine &line = lines[d][s];
      TH1D *h = FillRise(OfDay(GeRuns(), GE_DAYS[d]), Form("data_%d_%d", d, s),
                         RISE_SLICE_LO[s], RISE_SLICE_HI[s]);
      SetSourceLowTail(s, s);
      SetNominal();
      PeakFit capture = FitCaptureLine(h, Form("capture_%d_%d", d, s));
      g_capture_meas = capture.mu;

      SetNominal();
      TF1 *nominal = FitInRun(h, Form("nominal_%d_%d", d, s));
      Transferred nominal_ge = GainTransfer(line, nominal);
      Double_t weight = 1.0 / (nominal_ge.err * nominal_ge.err);
      SetTruth(truth);
      TF1 *truth_fit = FitInRun(h, Form("truth_%d_%d", d, s));
      const Double_t ge_truth = GainTransfer(line, truth_fit).e;

      SetNominal();
      std::vector<Double_t> deltas;
      for (Int_t t = 0; t < trials; t++) {
        TH1D *pseudo =
            Pseudo(h, truth_fit, random, Form("pseudo_%d_%d_%d", d, s, t));
        TF1 *f = FitInRun(pseudo, Form("fit_%d_%d_%d", d, s, t));
        deltas.push_back(GainTransfer(line, f).e - ge_truth);
        delete f;
        delete pseudo;
      }
      Double_t mean = 0, rms = 0;
      for (size_t k = 0; k < deltas.size(); k++)
        mean += deltas[k];
      mean /= deltas.size();
      for (size_t k = 0; k < deltas.size(); k++)
        rms += (deltas[k] - mean) * (deltas[k] - mean);
      rms = TMath::Sqrt(rms / (deltas.size() - 1));
      Double_t mean_err = rms / TMath::Sqrt((Double_t)deltas.size());
      std::cout << "TESTBED [" << TRUTH_NAMES[truth] << "] 10/" << GE_DAYS[d]
                << ", slice " << std::fixed << std::setprecision(1)
                << RISE_SLICE_LO[s] << "-" << RISE_SLICE_HI[s]
                << " ns: truth Ge " << std::setprecision(4) << ge_truth
                << " (nominal fit to the data " << nominal_ge.e
                << "), truth chi2/ndf " << std::setprecision(2)
                << truth_fit->GetChisquare() / truth_fit->GetNDF()
                << " | nominal fit to pseudo-data: bias " << std::showpos
                << std::setprecision(1) << 1000 * mean << std::noshowpos
                << " +- " << 1000 * mean_err << " eV, scatter " << 1000 * rms
                << " eV (fit error " << 1000 * nominal_ge.err << " eV), "
                << trials << " trials" << std::endl;
      Int_t point = bias[d]->GetN();
      bias[d]->SetPoint(
          point,
          0.5 * (TMath::Max(RISE_SLICE_LO[s], 250.0) + RISE_SLICE_HI[s]) +
              0.4 * (d - 0.5),
          1000 * mean);
      bias[d]->SetPointError(point, 0, 1000 * mean_err);
      sw += weight;
      swb += weight * mean;
      swb2 += weight * weight * mean_err * mean_err;
      delete nominal;
      delete truth_fit;
      delete h;
    }
  }
  std::cout << "TESTBED [" << TRUTH_NAMES[truth]
            << "] bias of the cells' weighted mean (nominal weights): "
            << std::showpos << std::fixed << std::setprecision(1)
            << 1000 * swb / sw << std::noshowpos << " +- "
            << 1000 * TMath::Sqrt(swb2) / sw << " eV" << std::endl;

  std::vector<TString> labels;
  for (Int_t d = 0; d < N_GE_DAYS; d++)
    labels.push_back(Form("10/%02d", GE_DAYS[d]));
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  const Int_t colors[N_GE_DAYS] = {kBlue + 1, kRed + 1};
  for (Int_t d = 0; d < N_GE_DAYS; d++) {
    PlottingUtils::ConfigureGraph(
        bias[d], colors[d],
        Form(";Rise time 10-90%% [ns];Ge fitted - truth [eV] (truth: %s)",
             TRUTH_NAMES[truth]));
    bias[d]->SetMarkerStyle(20 + d);
    bias[d]->GetXaxis()->SetLimits(250, 285);
    bias[d]->Draw(d == 0 ? "AP" : "P SAME");
    leg->AddEntry(bias[d], labels[d], "pe");
  }
  TLine *zero = new TLine(250, 0, 285, 0);
  zero->SetLineStyle(2);
  zero->Draw();
  PlottingUtils::SaveFigure(canvas, Form("cdte_testbed_%s", TRUTH_TAGS[truth]),
                            "", PlotSaveOptions::kLINEAR);
}
