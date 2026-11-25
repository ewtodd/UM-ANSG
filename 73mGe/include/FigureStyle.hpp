// Pieces of the figures the 73mGe macros draw outside AU's own fit plots,
// kept in AU's style.
#pragma once

#include "PlottingUtils.hpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TVirtualPad.h>

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
// axis over +-5.5, a dashed line at zero and dotted lines at +-3.
inline void DrawPullPanel(TGraph *pulls, Double_t x_lo, Double_t x_hi,
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
}

} // namespace FigureStyle
