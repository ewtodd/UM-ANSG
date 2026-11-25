// The production fits, the in-situ fits of the 01/13 signal runs and the 01/15
// rate-subtraction fit, redrawn for the paper and the dissertation from the
// objects CalibrateAndFit.cpp writes beside its own fit plots
// (PlottingUtils::SetFitFigureObjectsSaved). Each line is drawn
// as its whole shape, Gaussian core and low-side tail, on the continuum.
//
// Writes plots/czt/figures: each figure on a linear y axis, and the 01/13
// ones also, as log_*, on a log one.
//
// Usage: root -b -q 'FitFigure.cpp+'
#include "FigureStyle.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TPad.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TString.h>

#include <iostream>
#include <vector>

namespace {

const char *const PRODUCTION_TAG = "tau20_sigpb_esc_pxday";
const Double_t SHOW_LO_KEV = 64.0, SHOW_HI_KEV = 81.0;

// How RooFitUtils colours a fit's curves (RooFitUtils.cpp, PlotFit*Peak and
// AppendPeakGraphs), each curve being one component plus the background: the
// background green, each peak's Gaussian core black and its other parts in
// other colours, all of one peak in the same line style, and the fixed lines
// together in cyan.
const Color_t BACKGROUND_COLOR = kGreen;
const Color_t CORE_COLOR = kBlack;
const Color_t FIXED_LINES_COLOR = kCyan + 2;

struct FigureSpec {
  TString fit_name;    // CalibrateAndFit's plot name, without "log_"
  TString output_name; // this figure
  TString data_label;
};

struct FitObjects {
  TH1 *data = nullptr;
  TGraph *total = nullptr;
  TGraph *background = nullptr;
  TGraph *fixed_lines = nullptr;
  std::vector<TGraph *> peak_parts;
  Double_t fit_lo = 0, fit_hi = 0;
  TString chi2_label;
};

Bool_t ReadFitObjects(const TString &path, FitObjects &fit) {
  TFile file(path, "READ");
  if (file.IsZombie()) {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return kFALSE;
  }
  TH1 *data = dynamic_cast<TH1 *>(file.Get("data"));
  TGraph *total = dynamic_cast<TGraph *>(file.Get("total"));
  TParameter<Double_t> *lo =
      dynamic_cast<TParameter<Double_t> *>(file.Get("fit_range_low"));
  TParameter<Double_t> *hi =
      dynamic_cast<TParameter<Double_t> *>(file.Get("fit_range_high"));
  TNamed *label = dynamic_cast<TNamed *>(file.Get("label"));
  if (!data || !total || !lo || !hi || !label) {
    std::cerr << "ERROR: " << path << " lacks the fit figure's objects"
              << std::endl;
    return kFALSE;
  }
  data->SetDirectory(nullptr);
  fit.data = data;
  fit.total = total;
  fit.fit_lo = lo->GetVal();
  fit.fit_hi = hi->GetVal();
  fit.chi2_label = label->GetTitle();
  for (Int_t i = 0;; i++) {
    TGraph *g = dynamic_cast<TGraph *>(file.Get(Form("component_%d", i)));
    if (!g)
      break;
    if (g->GetLineColor() == BACKGROUND_COLOR)
      fit.background = g;
    else if (g->GetLineColor() == FIXED_LINES_COLOR)
      fit.fixed_lines = g;
    else
      fit.peak_parts.push_back(g);
  }
  if (!fit.background) {
    std::cerr << "ERROR: " << path << " has no background curve" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

// InitUtils turns on gROOT->ForceStyle(), which restyles every object read
// from a file and so erases the colours that tell the curves apart.
Bool_t LoadFit(const TString &path, FitObjects &fit) {
  Bool_t style_was_forced = gROOT->GetForceStyle();
  gROOT->ForceStyle(kFALSE);
  Bool_t loaded = ReadFitObjects(path, fit);
  gROOT->ForceStyle(style_was_forced);
  return loaded;
}

// One peak's whole shape on the background: every part in its line style,
// each with the background taken off once, then the background added back.
TGraph *PeakOnBackground(const FitObjects &fit, Style_t style) {
  Int_t n = fit.background->GetN();
  TGraph *g = new TGraph(n);
  for (Int_t k = 0; k < n; k++) {
    Double_t bkg = fit.background->GetY()[k];
    Double_t y = bkg;
    for (size_t p = 0; p < fit.peak_parts.size(); p++)
      if (fit.peak_parts[p]->GetLineStyle() == style)
        y += fit.peak_parts[p]->GetY()[k] - bkg;
    g->SetPoint(k, fit.background->GetX()[k], y);
  }
  return g;
}

Double_t PeakPositionKeV(const TGraph *core) {
  Double_t x_max = 0, y_max = -1e30;
  for (Int_t k = 0; k < core->GetN(); k++)
    if (core->GetY()[k] > y_max) {
      y_max = core->GetY()[k];
      x_max = core->GetX()[k];
    }
  return x_max;
}

void StyleCurve(TGraph *g, Color_t color, Style_t style) {
  PlottingUtils::ConfigureGraph(g, color);
  g->SetLineStyle(style);
}

void DrawFit(const TString &fits_dir, const FigureSpec &figure) {
  FitObjects fit;
  if (!LoadFit(fits_dir + "/" + figure.fit_name + ".root", fit))
    return;

  // The Ge line is the core peaking near 68.75 keV; every other peak is a Pb
  // K-alpha line.
  TGraph *ge = nullptr;
  std::vector<TGraph *> lead;
  for (size_t p = 0; p < fit.peak_parts.size(); p++) {
    if (fit.peak_parts[p]->GetLineColor() != CORE_COLOR)
      continue;
    TGraph *shape = PeakOnBackground(fit, fit.peak_parts[p]->GetLineStyle());
    Double_t position = PeakPositionKeV(fit.peak_parts[p]);
    std::cout << figure.output_name << ": peak at " << position << " keV"
              << std::endl;
    if (TMath::Abs(position - 68.75) < 1.0)
      ge = shape;
    else
      lead.push_back(shape);
  }
  if (!ge) {
    std::cerr << "ERROR: " << figure.fit_name << " has no Ge peak" << std::endl;
    return;
  }

  TH1 *data = fit.data;
  Double_t bin_ev = 1000.0 * data->GetBinWidth(1);
  data->GetXaxis()->SetRangeUser(SHOW_LO_KEV, SHOW_HI_KEV);
  Double_t data_max = data->GetMaximum();
  Double_t bkg_min =
      TMath::MinElement(fit.background->GetN(), fit.background->GetY());

  TGraph *pulls = new TGraph();
  for (Int_t b = 1; b <= data->GetNbinsX(); b++) {
    Double_t x = data->GetBinCenter(b);
    Double_t error = data->GetBinError(b);
    if (x < fit.fit_lo || x > fit.fit_hi || error <= 0)
      continue;
    pulls->SetPoint(pulls->GetN(), x,
                    (data->GetBinContent(b) - fit.total->Eval(x)) / error);
  }

  TString chi2_value = fit.chi2_label;
  chi2_value.Remove(0, chi2_value.Index("=") + 1);
  chi2_value = chi2_value.Strip(TString::kBoth);
  std::vector<TString> labels;
  labels.push_back(figure.data_label);
  labels.push_back("Fit, #chi^{2}/ndf = " + chi2_value);
  labels.push_back("^{73m}Ge, 68.75 keV");
  labels.push_back("Pb K#alpha_{2}, K#alpha_{1}");
  labels.push_back("^{113}Cd(n,#gamma) 95.9 keV escapes");

  // The legend sits beside the plot, centred on plot and panel together, so
  // the spectrum can fill the frame's height.
  const Int_t plot_h = 800, panel_h = 300;
  const Double_t data_share = 0.93;
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels),
      plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(plot, panel_h, upper, lower);
  upper->cd();
  PlottingUtils::ConfigureHistogram(
      data, kBlack, Form(";Energy [keV];Counts / %.0f eV", bin_ev));
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  // Counts in thousands, as on the figures with larger spectra.
  data->GetYaxis()->SetMaxDigits(3);
  const Double_t y_min = 0.85 * bkg_min;
  data->SetMinimum(y_min);
  data->SetMaximum(y_min + (data_max - y_min) / data_share);
  data->GetXaxis()->SetLabelSize(0);
  data->GetXaxis()->SetTitleSize(0);
  data->Draw("E");
  // The continuum is not drawn on its own: every line sits on it, and its
  // legend entry would push the legend onto the spectrum.
  if (fit.fixed_lines) {
    StyleCurve(fit.fixed_lines, kMagenta + 1, 1);
    fit.fixed_lines->Draw("L SAME");
  }
  for (size_t p = 0; p < lead.size(); p++) {
    StyleCurve(lead[p], kGreen + 2, 1);
    lead[p]->Draw("L SAME");
  }
  StyleCurve(ge, kAzure + 2, 1);
  ge->Draw("L SAME");
  StyleCurve(fit.total, kRed + 1, 1);
  fit.total->Draw("L SAME");
  data->Draw("E SAME");

  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  leg->AddEntry(data, labels[0], "lep");
  leg->AddEntry(fit.total, labels[1], "l");
  leg->AddEntry(ge, labels[2], "l");
  if (!lead.empty())
    leg->AddEntry(lead[0], labels[3], "l");
  if (fit.fixed_lines)
    leg->AddEntry(fit.fixed_lines, labels[4], "l");

  lower->cd();
  TH1D *pull_frame = FigureStyle::DrawPullPanel(pulls, SHOW_LO_KEV, SHOW_HI_KEV,
                                                "Energy [keV]");
  PlottingUtils::ScaleFigure(plot);
  // The pull panel's title beside its own short labels, not aligned with the
  // spectrum's.
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  PlottingUtils::SaveFigure(canvas, figure.output_name, "",
                            PlotSaveOptions::kLINEAR);

  // The log copy: labels at 1 and 5 x 10^n chosen from the data, its title
  // beside them.
  std::vector<TString> log_labels =
      FigureStyle::ConfigureLogYAxis(data, upper, y_min, data_max, data_share);
  FigureStyle::FitTitleToLabels(upper, data, log_labels);
  PlottingUtils::SaveFigure(canvas, figure.output_name, "",
                            PlotSaveOptions::kLOG);
}

// Free parameters of the 01/15 rate-subtraction fit (GainTransferFit.hpp,
// linear continuum): the Ge line's yield, centroid and width and the
// continuum's level and slope; its tail and the escapes are fixed.
const Int_t RATESUB_FREE_PARAMETERS = 5;

// The 01/15 Method 2 fit: the live-time residual of the signal and active
// background runs in the first of the binning phases whose fits are
// averaged, with the Ge line and the 95.902 keV line's escapes on the linear
// continuum. The residual goes below zero, so the figure is linear only.
void DrawRateSubFit(const TString &fits_dir, const TString &fit_name,
                    const TString &output_name, const TString &data_label) {
  TString path = fits_dir + "/" + fit_name + ".root";
  Bool_t style_was_forced = gROOT->GetForceStyle();
  gROOT->ForceStyle(kFALSE);
  TFile file(path, "READ");
  TH1 *data = dynamic_cast<TH1 *>(file.Get("data"));
  TGraph *total = dynamic_cast<TGraph *>(file.Get("total"));
  TGraph *continuum = dynamic_cast<TGraph *>(file.Get("component_0"));
  TGraph *escapes = dynamic_cast<TGraph *>(file.Get("component_1"));
  TParameter<Double_t> *lo_par =
      dynamic_cast<TParameter<Double_t> *>(file.Get("fit_range_low"));
  TParameter<Double_t> *hi_par =
      dynamic_cast<TParameter<Double_t> *>(file.Get("fit_range_high"));
  gROOT->ForceStyle(style_was_forced);
  if (!data || !total || !continuum || !escapes || !lo_par || !hi_par) {
    std::cerr << "ERROR: " << path
              << " lacks the rate-subtraction fit's objects" << std::endl;
    return;
  }
  data->SetDirectory(nullptr);
  const Double_t fit_lo = lo_par->GetVal(), fit_hi = hi_par->GetVal();

  // The Ge line on the continuum: the total less the escapes above it.
  TGraph *ge = new TGraph(total->GetN());
  for (Int_t k = 0; k < total->GetN(); k++)
    ge->SetPoint(k, total->GetX()[k],
                 total->GetY()[k] -
                     (escapes->GetY()[k] - continuum->GetY()[k]));

  Double_t chi2 = 0, y_lo = 1e30, y_hi = -1e30;
  Int_t n_bins = 0;
  TGraph *pulls = new TGraph();
  for (Int_t b = 1; b <= data->GetNbinsX(); b++) {
    Double_t x = data->GetBinCenter(b);
    Double_t y = data->GetBinContent(b), error = data->GetBinError(b);
    if (x < fit_lo || x > fit_hi || error <= 0)
      continue;
    Double_t pull = (y - total->Eval(x)) / error;
    chi2 += pull * pull;
    n_bins++;
    pulls->SetPoint(pulls->GetN(), x, pull);
    y_lo = TMath::Min(y_lo, y - error);
    y_hi = TMath::Max(y_hi, y + error);
  }
  const Int_t ndf = n_bins - RATESUB_FREE_PARAMETERS;
  std::cout << output_name << ": chi2 " << chi2 << " / " << ndf << std::endl;

  std::vector<TString> labels;
  labels.push_back(data_label);
  labels.push_back(Form("Fit, #chi^{2}/ndf = %.3f", chi2 / ndf));
  labels.push_back("^{73m}Ge, 68.75 keV");
  labels.push_back("^{113}Cd(n,#gamma) 95.9 keV escapes");

  const Int_t plot_h = 800, panel_h = 300;
  const Double_t data_share = 0.93;
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels),
      plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(plot, panel_h, upper, lower);
  upper->cd();
  PlottingUtils::ConfigureHistogram(
      data, kBlack,
      Form(";Energy [keV];Counts/s / %.0f eV", 1000.0 * data->GetBinWidth(1)));
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  data->GetXaxis()->SetRangeUser(fit_lo, fit_hi);
  const Double_t frame_lo = y_lo - 0.05 * (y_hi - y_lo);
  data->SetMinimum(frame_lo);
  data->SetMaximum(frame_lo + (y_hi - frame_lo) / data_share);
  data->GetXaxis()->SetLabelSize(0);
  data->GetXaxis()->SetTitleSize(0);
  data->Draw("E");
  StyleCurve(escapes, kMagenta + 1, 1);
  escapes->Draw("L SAME");
  StyleCurve(ge, kAzure + 2, 1);
  ge->Draw("L SAME");
  StyleCurve(total, kRed + 1, 1);
  total->Draw("L SAME");
  data->Draw("E SAME");

  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  leg->AddEntry(data, labels[0], "lep");
  leg->AddEntry(total, labels[1], "l");
  leg->AddEntry(ge, labels[2], "l");
  leg->AddEntry(escapes, labels[3], "l");

  lower->cd();
  TH1D *pull_frame =
      FigureStyle::DrawPullPanel(pulls, fit_lo, fit_hi, "Energy [keV]");
  PlottingUtils::ScaleFigure(plot);
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  PlottingUtils::SaveFigure(canvas, output_name, "", PlotSaveOptions::kLINEAR);
}

} // namespace

void FitFigure() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/czt/figures",
                                project_root + "/root_files");
  TString fits_dir =
      project_root + "/plots/czt/chain/" + PRODUCTION_TAG + "/fits";
  std::vector<FigureSpec> figures;
  FigureSpec cu10 = {
      "PostCal_CuShield10_sig_CuShieldSignal_10Percent_20260113_postcal_"
      "CuShield10",
      "czt_fit_cu10", "Signal run, Cu shield, 10% power"};
  FigureSpec cd10 = {
      "PostCal_CdShield10_sig_CdShieldSignal_10Percent_20260113_postcal_"
      "CdShield10",
      "czt_fit_cd10", "Signal run, Cd shield, 10% power"};
  FigureSpec cd25 = {
      "PostCal_CdShield25_sig_CdShieldSignal_25Percent_20260113_postcal_"
      "CdShield25",
      "czt_fit_cd25", "Signal run, Cd shield, 25% power"};
  figures.push_back(cu10);
  figures.push_back(cd10);
  figures.push_back(cd25);
  for (size_t i = 0; i < figures.size(); i++)
    DrawFit(fits_dir, figures[i]);
  DrawRateSubFit(fits_dir, "ratesub_fit_20260115", "czt_fit_0115",
                 "Background-subtracted signal run, no shield, 5% power");
}
