// The production fits, the in-situ fits of the 01/13 signal runs and the direct
// fit of the 01/15 no-shield run, redrawn for the paper and the dissertation
// from the objects CalibrateAndFit.cpp writes beside its own fit plots
// (PlottingUtils::SetFitFigureObjectsSaved). Each line is drawn
// as its whole shape, Gaussian core and low-side tail, on the continuum.
//
// Writes plots/czt/figures: each figure on a linear y axis and, as log_*, on
// a log one, at the letter's width and, as wide_*, wider for the
// dissertation's full-width panels (two to a page).
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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const char *const PRODUCTION_TAG = "tau20_sigpb_esc_pxday";
// Each run's Ge energy is quoted, as CombineGeResult.cpp and Table 5.1 quote
// it, at the midpoint of the gain-transfer fits (PRODUCTION_TAG, drawn here)
// and the offset-transfer ones.
const char *const OFFSET_TRANSFER_TAG = "tau20_sigpb_esc_oxfer_pxday";
// The figures show each fit's window and this much on either side.
const Double_t SHOW_MARGIN_KEV = 1.0;

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
  // The linear copy's axis starts just under the data rather than at 0.85
  // of the continuum, for a spectrum whose lines rise little above it.
  Bool_t tight_linear_axis;
  // The fit's row in CalibrateAndFit's results/ge_<tag>.result files, whose
  // Ge centroids' midpoint (OFFSET_TRANSFER_TAG) the legend quotes.
  TString result_label;
};

// The Ge centroid, keV, of the row labelled label in a ge_<tag>.result file
// (columns: method, ge_mu, ..., label); 0 if there is none.
Double_t ReadGeCentroid(const TString &path, const TString &label) {
  std::ifstream in(path.Data());
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token)
      tokens.push_back(token);
    if (tokens.size() >= 3 && TString(tokens.back()) == label)
      return std::atof(tokens[1].c_str());
  }
  return 0;
}

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

// The fit's curves at twice the style's width, to stand out from the data.
void StyleCurve(TGraph *g, Color_t color, Style_t style) {
  PlottingUtils::ConfigureGraph(g, color);
  g->SetLineStyle(style);
  g->SetLineWidth(2 * PlottingUtils::GetLineWidth());
}

// Canvas widths, pixels: the letter's single column, and the dissertation's
// full text width with two panels to a page.
const Int_t LETTER_CANVAS_W = 1400;
const Int_t DISSERTATION_CANVAS_W = 2000;

// The legend sits inside the frame's top left, over headroom left above the
// data, so the figure reads at a single column's width. A canvas wider than
// the letter's keeps the letter's pixel margins and text, so the extra width
// all goes to the frame; its files are named with prefix.
void DrawFit(const TString &fits_dir, const TString &gain_result_path,
             const TString &offset_result_path, const FigureSpec &figure,
             Int_t canvas_w, const TString &prefix) {
  FitObjects fit;
  if (!LoadFit(fits_dir + "/" + figure.fit_name + ".root", fit))
    return;
  const Double_t gain_kev =
      ReadGeCentroid(gain_result_path, figure.result_label);
  const Double_t offset_kev =
      ReadGeCentroid(offset_result_path, figure.result_label);
  if (gain_kev <= 0 || offset_kev <= 0) {
    std::cerr << "ERROR: no row " << figure.result_label << " in "
              << gain_result_path << " or " << offset_result_path << std::endl;
    return;
  }
  const Double_t ge_kev = 0.5 * (gain_kev + offset_kev);

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
  const Double_t show_lo = fit.fit_lo - SHOW_MARGIN_KEV;
  const Double_t show_hi = fit.fit_hi + SHOW_MARGIN_KEV;
  data->GetXaxis()->SetRangeUser(show_lo, show_hi);
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
  labels.push_back(TString::Format("^{73}Ge, %.4f keV", ge_kev));
  const TString lead_label = "Pb K#alpha_{2}, K#alpha_{1}";
  const TString escape_label = "^{113}Cd(n,#gamma) 95.9 keV escapes";
  if (!lead.empty())
    labels.push_back(lead_label);
  if (fit.fixed_lines)
    labels.push_back(escape_label);

  // The legend takes the top two fifths of the frame.
  const Int_t plot_h = 800, panel_h = 300;
  const Double_t data_share = 0.6;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(canvas_w, plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(canvas, panel_h, upper, lower);
  const Bool_t wider = canvas_w != LETTER_CANVAS_W;
  if (wider) {
    const Double_t margin_share = Double_t(LETTER_CANVAS_W) / canvas_w;
    TPad *pads[2] = {upper, lower};
    for (Int_t p = 0; p < 2; p++) {
      pads[p]->SetLeftMargin(pads[p]->GetLeftMargin() * margin_share);
      pads[p]->SetRightMargin(pads[p]->GetRightMargin() * margin_share);
    }
  }
  upper->cd();
  PlottingUtils::ConfigureHistogram(
      data, kBlack, Form(";Energy [keV];Counts / %.0f eV", bin_ev));
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  // Counts in thousands, as on the figures with larger spectra. ROOT writes
  // a narrow range of six-digit counts as 0.145 x 10^6 unless it may keep a
  // fourth digit.
  data->GetYaxis()->SetMaxDigits(figure.tight_linear_axis ? 4 : 3);
  // The log copy starts at the same minimum; ConfigureLogYAxis lowers it to a
  // labelled power of ten where that is near enough.
  const Double_t y_min = figure.tight_linear_axis
                             ? bkg_min - 0.15 * (data_max - bkg_min)
                             : 0.85 * bkg_min;
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

  TLegend *leg = FigureStyle::AddFrameLegend(upper, labels, kTRUE);
  leg->AddEntry(data, labels[0], "lep");
  leg->AddEntry(fit.total, labels[1], "l");
  leg->AddEntry(ge, labels[2], "l");
  if (!lead.empty())
    leg->AddEntry(lead[0], lead_label, "l");
  if (fit.fixed_lines)
    leg->AddEntry(fit.fixed_lines, escape_label, "l");

  lower->cd();
  TH1D *pull_frame =
      FigureStyle::DrawPullPanel(pulls, show_lo, show_hi, "Energy [keV]");
  pulls->SetMarkerSize(1.1);
  PlottingUtils::ScaleFigure(canvas, canvas_w, plot_h);
  // The pull panel's title beside its own short labels, not aligned with the
  // spectrum's.
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  // On a wider canvas ROOT's own title offset grows with the pad's width and
  // would leave the narrower margin; set it beside labels of up to three
  // digits instead.
  if (wider) {
    std::vector<TString> linear_labels;
    linear_labels.push_back("000");
    FigureStyle::FitTitleToLabels(upper, data, linear_labels);
  }
  const TString output_name = prefix + figure.output_name;
  PlottingUtils::SaveFigure(canvas, output_name, "", PlotSaveOptions::kLINEAR);
  // The log copy: labels at 1 and 5 x 10^n chosen from the data, its title
  // beside them.
  std::vector<TString> log_labels =
      FigureStyle::ConfigureLogYAxis(data, upper, y_min, data_max, data_share);
  FigureStyle::FitTitleToLabels(upper, data, log_labels);
  PlottingUtils::SaveFigure(canvas, output_name, "", PlotSaveOptions::kLOG);
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
  const TString gain_result_path =
      project_root + "/results/ge_" + PRODUCTION_TAG + ".result";
  const TString offset_result_path =
      project_root + "/results/ge_" + OFFSET_TRANSFER_TAG + ".result";
  FigureSpec cu10 = {
      "PostCal_CuShield10_sig_CuShieldSignal_10Percent_20260113_postcal_"
      "CuShield10",
      "czt_fit_cu10", "Signal run, Cu shield, 10% power", kFALSE,
      "Cu_Shield_Signal_10%_(01/13)"};
  FigureSpec cd10 = {
      "PostCal_CdShield10_sig_CdShieldSignal_10Percent_20260113_postcal_"
      "CdShield10",
      "czt_fit_cd10", "Signal run, Cd shield, 10% power", kFALSE,
      "Cd_Shield_Signal_10%_(01/13)"};
  FigureSpec cd25 = {
      "PostCal_CdShield25_sig_CdShieldSignal_25Percent_20260113_postcal_"
      "CdShield25",
      "czt_fit_cd25", "Signal run, Cd shield, 25% power", kFALSE,
      "Cd_Shield_Signal_25%_(01/13)"};
  // Method 2's direct fit of the no-shield run, whose Ge line rises about 10%
  // above the continuum.
  FigureSpec noshield = {
      "Direct_20260115_direct_NoShieldSignal_5Percent_20260115", "czt_fit_0115",
      "Signal run, no shield, 5% power", kTRUE,
      "No_Shield_Signal_5%_(01/15)_[direct]"};
  figures.push_back(cu10);
  figures.push_back(cd10);
  figures.push_back(cd25);
  figures.push_back(noshield);
  for (size_t i = 0; i < figures.size(); i++) {
    DrawFit(fits_dir, gain_result_path, offset_result_path, figures[i],
            LETTER_CANVAS_W, "");
    DrawFit(fits_dir, gain_result_path, offset_result_path, figures[i],
            DISSERTATION_CANVAS_W, "wide_");
  }
}
