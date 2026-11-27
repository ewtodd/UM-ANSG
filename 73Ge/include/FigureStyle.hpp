// Pieces of the figures the 73Ge macros draw outside AU's own fit plots,
// kept in AU's style.
#pragma once

#include "PlottingUtils.hpp"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
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

// The round values inside [lo, hi] with the largest step of 1, 2 or 5 x 10^k
// that gives at least two; empty if none does.
inline std::vector<Double_t> RoundLabelValues(Double_t lo, Double_t hi) {
  std::vector<Double_t> values;
  const Int_t steps_per_decade = 3;
  const Double_t step_multiples[steps_per_decade] = {5.0, 2.0, 1.0};
  Int_t top_decade = Int_t(TMath::Floor(TMath::Log10(hi - lo)));
  for (Int_t decade = top_decade; decade >= top_decade - 3; decade--) {
    for (Int_t s = 0; s < steps_per_decade; s++) {
      Double_t step = step_multiples[s] * TMath::Power(10.0, decade);
      values.clear();
      for (Double_t v = TMath::Ceil(lo / step - 1e-9) * step; v <= hi;
           v += step)
        values.push_back(v);
      if (values.size() >= 2)
        return values;
    }
  }
  values.clear();
  return values;
}

// The text of a round label as its multiple of 10^decade, with as many
// decimals as the step between labels (a power of ten times 1, 2 or 5) asks
// for; the 10^decade is written once at the top of the axis.
inline TString RoundLabelText(Double_t value, Double_t step, Int_t decade) {
  Int_t step_decade = Int_t(TMath::Floor(TMath::Log10(step) + 1e-9));
  Int_t decimals = TMath::Max(0, decade - step_decade);
  return Form("%.*f", decimals, value / TMath::Power(10.0, decade));
}

// A log y axis inside one decade that holds no k x 10^n: ROOT draws neither
// ticks nor labels there, so the labels at the round values and their ticks
// (minor ones halfway between) are drawn here, on both sides of the frame.
inline std::vector<TString> DrawRoundLogLabels(TH1 *frame, TVirtualPad *pad,
                                               Double_t lo, Double_t hi) {
  std::vector<TString> texts;
  std::vector<Double_t> values = RoundLabelValues(lo, hi);
  if (values.size() < 2)
    return texts;
  const Double_t step = values[1] - values[0];
  const Int_t decade = Int_t(TMath::Floor(TMath::Log10(values.back()) + 1e-9));
  TAxis *axis = frame->GetYaxis();
  const Double_t left = pad->GetLeftMargin();
  const Double_t right = 1.0 - pad->GetRightMargin();
  const Double_t bottom = pad->GetBottomMargin();
  const Double_t top = 1.0 - pad->GetTopMargin();
  const Double_t major_tick = axis->GetTickLength() * (right - left);
  const Double_t minor_tick = 0.5 * major_tick;
  const Double_t label_x = left - axis->GetLabelOffset();
  pad->cd();
  // Every half step from half a step below the first label to half a step
  // above the last: labels at the odd positions, minor ticks between.
  const Int_t n_positions = 2 * (Int_t)values.size() + 1;
  for (Int_t i = 0; i < n_positions; i++) {
    Double_t v = values[0] + (i - 1) * 0.5 * step;
    if (v < lo || v > hi)
      continue;
    Double_t y =
        bottom + (top - bottom) * TMath::Log10(v / lo) / TMath::Log10(hi / lo);
    Bool_t is_label = (i % 2 == 1);
    Double_t length = is_label ? major_tick : minor_tick;
    TLine *left_tick = new TLine();
    left_tick->SetBit(kCanDelete);
    left_tick->DrawLineNDC(left, y, left + length, y);
    TLine *right_tick = new TLine();
    right_tick->SetBit(kCanDelete);
    right_tick->DrawLineNDC(right - length, y, right, y);
    if (!is_label)
      continue;
    TString text = RoundLabelText(v, step, decade);
    TLatex *label = new TLatex();
    label->SetBit(kCanDelete);
    label->SetTextFont(axis->GetLabelFont());
    label->SetTextSize(axis->GetLabelSize());
    label->SetTextAlign(32);
    label->DrawLatexNDC(label_x, y, text);
    texts.push_back(text);
  }
  // The common power of ten above the axis, where ROOT puts a linear axis's.
  TLatex *exponent = new TLatex();
  exponent->SetBit(kCanDelete);
  exponent->SetTextFont(axis->GetLabelFont());
  exponent->SetTextSize(axis->GetLabelSize());
  exponent->SetTextAlign(11);
  exponent->DrawLatexNDC(left, top + 0.01, Form("#times10^{%d}", decade));
  return texts;
}

// Sets a log y axis's range and labels from the data. The range runs from
// natural_min with data_max at data_share of the frame; the labels are at 1
// and 5 x 10^n. When that range holds fewer than two, the minimum drops to
// the next such value if that lowers it by at most 20%, leaving the label
// room inside the frame; if there are still fewer than two, the labels
// follow 1-2-5; failing that, the one label there is stays. A range that
// holds no such value at all (one inside a single decade) gets labels at
// round values instead (DrawRoundLogLabels).
// Returns the labels drawn, for FitTitleToLabels().
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
  // No label at all: the range lies inside one decade, and lowering never
  // reached a label (lowering always leaves one).
  Bool_t round_labels = shown.empty();
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
  if (round_labels)
    texts = DrawRoundLogLabels(frame, pad, lo, hi);
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
