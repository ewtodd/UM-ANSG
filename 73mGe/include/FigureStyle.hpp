// Pieces of the figures the 73mGe macros draw outside AU's own fit plots,
// kept in AU's style.
#pragma once

#include "PlottingUtils.hpp"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <vector>

namespace FigureStyle {

// A legend inside a pad's frame in one of its corners, sized for
// AddLegend()'s font in 50 px rows in the pad's own NDC. Works on a canvas, a
// pad of SplitPadForPanel(), or its panel. With symbol_column_px set, the
// symbols get a column that wide instead of TLegend's quarter of the box,
// which leaves a short line in a wide strip beside long labels.
inline TLegend *AddFrameLegend(TVirtualPad *pad,
                               const std::vector<TString> &labels,
                               Bool_t at_left, Bool_t at_top = kTRUE,
                               Double_t symbol_column_px = 0) {
  pad->cd();
  TCanvas *canvas = pad->GetCanvas();
  Double_t pad_w_px = pad->GetWNDC() * canvas->GetWw();
  Double_t pad_h_px = pad->GetHNDC() * canvas->GetWh();
  const Double_t padding_px = 15.0;
  Double_t width_px = PlottingUtils::LegendWidthPx(labels);
  if (symbol_column_px > 0) {
    Double_t text_px = 0;
    for (size_t i = 0; i < labels.size(); i++)
      text_px = TMath::Max(text_px, PlottingUtils::TextWidthPx(labels[i]));
    width_px = symbol_column_px + text_px + 2.0 * padding_px;
  }
  Double_t width = width_px / pad_w_px;
  Double_t height = 50.0 * labels.size() / pad_h_px;
  Double_t top = 1.0 - pad->GetTopMargin() - 0.02;
  if (!at_top)
    top = pad->GetBottomMargin() + 0.02 + height;
  Double_t left = pad->GetLeftMargin() + 0.01;
  if (!at_left)
    left = 1.0 - pad->GetRightMargin() - 0.01 - width;
  TLegend *legend =
      PlottingUtils::AddLegend(left, left + width, top - height, top);
  if (symbol_column_px > 0)
    legend->SetMargin(symbol_column_px / width_px);
  return legend;
}

// Draws, on the current pad, the pull panel of AU's
// PlottingUtils::PlotFitWithResiduals(): azure points, a centred delta/sigma
// axis over +-5.5, a dashed line at zero and dotted lines at +-3. Returns
// the panel's frame, whose y labels are PullPanelLabels().
inline TH1D *DrawPullPanel(TGraph *pulls, Double_t x_lo, Double_t x_hi,
                           const TString &x_title) {
  static Int_t frame_count = 0;
  TH1D *frame =
      new TH1D(Form("pull_frame_%d", frame_count++), "", 100, x_lo, x_hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(frame, kBlack,
                                    ";" + x_title + ";#delta/#sigma");
  frame->SetMinimum(-5.5);
  frame->SetMaximum(5.5);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->CenterTitle(kTRUE);
  frame->Draw("AXIS");
  TLine *zero = new TLine(x_lo, 0, x_hi, 0);
  zero->SetLineStyle(2);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();
  for (Int_t sign = -1; sign <= 1; sign += 2) {
    TLine *three_sigma = new TLine(x_lo, 3 * sign, x_hi, 3 * sign);
    three_sigma->SetLineColor(kGray + 2);
    three_sigma->SetLineStyle(3);
    three_sigma->SetLineWidth(PlottingUtils::GetLineWidth());
    three_sigma->Draw();
  }
  pulls->SetMarkerStyle(20);
  pulls->SetMarkerSize(0.8);
  pulls->SetMarkerColor(kAzure);
  pulls->SetLineColor(kAzure);
  pulls->Draw("P SAME");
  return frame;
}

// The y labels DrawPullPanel()'s axis shows.
inline std::vector<TString> PullPanelLabels() {
  std::vector<TString> labels;
  labels.push_back("5");
  labels.push_back("0");
  labels.push_back("#minus5");
  return labels;
}

// Values k x 10^n with k in multiples, inside [lo, hi], ascending.
inline std::vector<Double_t>
LogLabelValues(Double_t lo, Double_t hi, const std::vector<Int_t> &multiples) {
  std::vector<Double_t> values;
  for (Int_t decade = Int_t(TMath::Floor(TMath::Log10(lo))) - 1;
       decade <= Int_t(TMath::Ceil(TMath::Log10(hi))); decade++) {
    for (size_t m = 0; m < multiples.size(); m++) {
      Double_t value = multiples[m] * TMath::Power(10.0, decade);
      if (value >= lo && value <= hi)
        values.push_back(value);
    }
  }
  std::sort(values.begin(), values.end());
  return values;
}

// The largest k x 10^n below value, k in multiples.
inline Double_t LogLabelBelow(Double_t value,
                              const std::vector<Int_t> &multiples) {
  Double_t below = 0;
  Int_t top_decade = Int_t(TMath::Floor(TMath::Log10(value)));
  for (Int_t decade = top_decade - 2; decade <= top_decade; decade++) {
    for (size_t m = 0; m < multiples.size(); m++) {
      Double_t candidate = multiples[m] * TMath::Power(10.0, decade);
      if (candidate < value && candidate > below)
        below = candidate;
    }
  }
  return below;
}

// The maximum that puts data_max at data_share of a log frame from lo.
inline Double_t LogMaxForShare(Double_t lo, Double_t data_max,
                               Double_t data_share) {
  return lo * TMath::Power(data_max / lo, 1.0 / data_share);
}

// ROOT's text for a log-axis label at k x 10^n.
inline TString LogLabelText(Double_t value) {
  Int_t decade = Int_t(TMath::Floor(TMath::Log10(value) + 1e-9));
  Int_t multiple = TMath::Nint(value / TMath::Power(10.0, decade));
  if (multiple == 1)
    return Form("10^{%d}", decade);
  return Form("%d#times10^{%d}", multiple, decade);
}

// Sets a log y axis's range and labels from the data. The range runs from
// natural_min with data_max at data_share of the frame; the labels are at 1
// and 5 x 10^n. When that range holds fewer than two, the minimum drops to
// the next such value if that lowers it by at most 20%, leaving the label
// room inside the frame; if there are still fewer than two, the labels
// follow 1-2-5; failing that, the one label there is stays. Returns the
// labels drawn, for FitTitleToLabels().
inline std::vector<TString> ConfigureLogYAxis(TH1 *frame, TVirtualPad *pad,
                                              Double_t natural_min,
                                              Double_t data_max,
                                              Double_t data_share) {
  std::vector<Int_t> one_five;
  one_five.push_back(1);
  one_five.push_back(5);
  std::vector<Int_t> one_two_five;
  one_two_five.push_back(1);
  one_two_five.push_back(2);
  one_two_five.push_back(5);

  Double_t lo = natural_min;
  Double_t hi = LogMaxForShare(lo, data_max, data_share);
  std::vector<Double_t> shown = LogLabelValues(lo, hi, one_five);
  if (shown.size() < 2) {
    Double_t below = LogLabelBelow(lo, one_five);
    if (below >= 0.8 * natural_min) {
      // Half a label's height between the lowered label and the frame's foot.
      const Double_t pad_h_px = pad->GetAbsHNDC() * pad->GetCanvas()->GetWh();
      const Double_t pad_w_px = pad->GetAbsWNDC() * pad->GetCanvas()->GetWw();
      const Double_t frame_px =
          pad_h_px * (1.0 - pad->GetTopMargin() - pad->GetBottomMargin());
      const Double_t label_px =
          frame->GetYaxis()->GetLabelSize() * TMath::Min(pad_h_px, pad_w_px);
      lo = below;
      for (Int_t iteration = 0; iteration < 5; iteration++) {
        hi = LogMaxForShare(lo, data_max, data_share);
        Double_t decades = TMath::Log10(hi / lo);
        lo = below * TMath::Power(10.0, -0.5 * label_px * decades / frame_px);
      }
      hi = LogMaxForShare(lo, data_max, data_share);
      shown = LogLabelValues(lo, hi, one_five);
    }
    if (shown.size() < 2) {
      std::vector<Double_t> series = LogLabelValues(lo, hi, one_two_five);
      if (series.size() > shown.size())
        shown = series;
    }
  }
  frame->SetMinimum(lo);
  frame->SetMaximum(hi);
  TAxis *axis = frame->GetYaxis();
  axis->SetMoreLogLabels(kTRUE);
  std::vector<Int_t> all_multiples;
  for (Int_t k = 1; k <= 9; k++)
    all_multiples.push_back(k);
  std::vector<Double_t> candidates = LogLabelValues(lo, hi, all_multiples);
  std::vector<TString> texts;
  for (size_t c = 0; c < candidates.size(); c++) {
    Bool_t keep = kFALSE;
    for (size_t s = 0; s < shown.size(); s++) {
      if (TMath::Abs(candidates[c] / shown[s] - 1.0) < 1e-6)
        keep = kTRUE;
    }
    if (keep)
      texts.push_back(LogLabelText(candidates[c]));
    else
      axis->ChangeLabelByValue(candidates[c], -1, 0);
  }
  return texts;
}

// Places a y title from the labels its axis draws: the label offset, the
// widest label, a gap of half the title's height, then the title, whose
// centre ROOT puts 1.6 x title size x offset x pad width from the axis (for
// fonts sized relative to the pad). Call it after PlottingUtils::ScaleFigure,
// which resets title offsets.
inline void FitTitleToLabels(TVirtualPad *pad, TH1 *frame,
                             const std::vector<TString> &labels) {
  TAxis *axis = frame->GetYaxis();
  const Double_t pad_w_px = pad->GetAbsWNDC() * pad->GetCanvas()->GetWw();
  const Double_t pad_h_px = pad->GetAbsHNDC() * pad->GetCanvas()->GetWh();
  const Double_t min_px = TMath::Min(pad_w_px, pad_h_px);
  const Double_t label_px = axis->GetLabelSize() * min_px;
  const Double_t title_px = axis->GetTitleSize() * min_px;
  Double_t widest_px = 0;
  for (size_t t = 0; t < labels.size(); t++)
    widest_px = TMath::Max(widest_px,
                           PlottingUtils::TextWidthPx(labels[t], label_px, 43));
  const Double_t frame_w_px =
      pad_w_px * (1.0 - pad->GetLeftMargin() - pad->GetRightMargin());
  const Double_t centre_px = axis->GetLabelOffset() * frame_w_px + widest_px +
                             0.5 * title_px + 0.5 * title_px;
  axis->SetTitleOffset(centre_px / (1.6 * axis->GetTitleSize() * pad_w_px));
}

} // namespace FigureStyle
