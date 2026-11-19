// The recommended energy of the 68.75 keV gamma ray from the measurements
// made without a Ge detector (STATUS.md, Plan): this work's CZT, the
// collaborator's CdTe and Si, and Jones and Kraner's Si(Li) values.
//
// Jones and Kraner calibrated on 241Am at 59.537 keV and Pb K-alpha at
// 72.8042 and 74.9694 keV. On DDEP's and NIST SRD 128's values their
// energies shift by an amount that depends on how they built the
// calibration, which their papers do not give. The plots put their lines on
// a fake ADC, exactly linear in energy on their scale, and recalibrate it on
// today's values every way their three lines allow.
#include "FigureStyle.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

const Int_t N_REF = 3;
const char *const REF_NAMES[N_REF] = {"Am-241", "Pb Ka2", "Pb Ka1"};
// The energies Jones and Kraner calibrated on, keV (Phys. Rev. C 4, 125;
// Phys. Rev. A 11, 1347).
const Double_t JK_REF_KEV[N_REF] = {59.537, 72.8042, 74.9694};
const Double_t CURRENT_REF_KEV[N_REF] = {
    ReferenceLines::AM241, ReferenceLines::PB_KA2, ReferenceLines::PB_KA1};
// Their 1975 value, read on their calibration.
const Double_t GE_JK_KEV = 68.755;

// The fake ADC: on Jones and Kraner's calibration every line reads back the
// energy it was put in at.
const Double_t ADC_CHANNELS_PER_KEV = 40.0;
Double_t Channel(Double_t jk_kev) { return ADC_CHANNELS_PER_KEV * jk_kev; }
Double_t JKEnergyKeV(Double_t channel) {
  return channel / ADC_CHANNELS_PER_KEV;
}

// The ADC calibrated on today's reference energies: a weighted least-squares
// line through the three lines' channels (two nonzero weights give the line
// through those two), or, with all weights zero, the quadratic through all
// three.
Double_t CalibratedKeV(const Double_t weights[N_REF], Double_t channel) {
  if (weights[0] == 0 && weights[1] == 0 && weights[2] == 0) {
    Double_t e = 0;
    for (Int_t i = 0; i < N_REF; i++) {
      Double_t lagrange = 1.0;
      for (Int_t j = 0; j < N_REF; j++)
        if (j != i)
          lagrange *= (channel - Channel(JK_REF_KEV[j])) /
                      (Channel(JK_REF_KEV[i]) - Channel(JK_REF_KEV[j]));
      e += lagrange * CURRENT_REF_KEV[i];
    }
    return e;
  }
  Double_t sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (Int_t i = 0; i < N_REF; i++) {
    Double_t w = weights[i];
    Double_t x = Channel(JK_REF_KEV[i]);
    Double_t y = CURRENT_REF_KEV[i];
    sw += w;
    swx += w * x;
    swy += w * y;
    swxx += w * x * x;
    swxy += w * x * y;
  }
  Double_t slope = (sw * swxy - swx * swy) / (sw * swxx - swx * swx);
  Double_t intercept = (swy - slope * swx) / sw;
  return intercept + slope * channel;
}

// How far a calibration on today's values reads a channel above Jones and
// Kraner's, in eV.
Double_t ShiftEV(const Double_t weights[N_REF], Double_t channel) {
  return 1000.0 * (CalibratedKeV(weights, channel) - JKEnergyKeV(channel));
}

// The lowest and highest shift at a channel over weighted least-squares
// lines, the weights spanning a factor 10^4 each way.
void LinearEnvelopeEV(Double_t channel, Double_t &low, Double_t &high) {
  low = 1e30;
  high = -1e30;
  for (Int_t a = -16; a <= 16; a++)
    for (Int_t b = -16; b <= 16; b++) {
      Double_t weights[N_REF] = {std::pow(10.0, a / 4.0),
                                 std::pow(10.0, b / 4.0), 1.0};
      Double_t s = ShiftEV(weights, channel);
      low = std::min(low, s);
      high = std::max(high, s);
    }
}

struct Calibration {
  TString name;
  Double_t weights[N_REF];
  Int_t color;
  Int_t marker;
};

struct Input {
  TString name;
  Double_t mu, err; // keV
  Bool_t jones_kraner;
};

TGraph *ShiftCurve(const Calibration &cal, Double_t ch_lo, Double_t ch_hi) {
  const Int_t n_points = 400;
  TGraph *g = new TGraph(n_points);
  for (Int_t k = 0; k < n_points; k++) {
    Double_t ch = ch_lo + (ch_hi - ch_lo) * k / (n_points - 1);
    g->SetPoint(k, ch, ShiftEV(cal.weights, ch));
  }
  g->SetLineColor(cal.color);
  g->SetLineWidth(3);
  return g;
}

TGraph *EnvelopeBand(Double_t ch_lo, Double_t ch_hi) {
  const Int_t n_points = 200;
  TGraph *band = new TGraph(2 * n_points);
  for (Int_t k = 0; k < n_points; k++) {
    Double_t ch = ch_lo + (ch_hi - ch_lo) * k / (n_points - 1);
    Double_t low = 0, high = 0;
    LinearEnvelopeEV(ch, low, high);
    band->SetPoint(k, ch, high);
    band->SetPoint(2 * n_points - 1 - k, ch, low);
  }
  band->SetFillColorAlpha(kGray, 0.5);
  band->SetLineWidth(0);
  return band;
}

// Today's calibrations minus Jones and Kraner's over a channel range; with
// mark_ge, each calibration's reading of the 68.75 keV line marked on its
// channel. Hands back the drawn objects for the legend.
void DrawShiftPanel(const std::vector<Calibration> &cals, Double_t ch_lo,
                    Double_t ch_hi, Double_t y_lo, Double_t y_hi,
                    const TString &frame_name, Bool_t mark_ge,
                    std::vector<TGraph *> &curves, TGraph *&band,
                    TGraph *&ref_changes, std::vector<TGraph *> &ge_markers) {
  TH1D *frame = new TH1D(frame_name, ";ADC channel;Today's #minus J&K's [eV]",
                         100, ch_lo, ch_hi);
  frame->SetMinimum(y_lo);
  frame->SetMaximum(y_hi);
  frame->SetStats(kFALSE);
  frame->GetXaxis()->SetNdivisions(505);
  frame->Draw("AXIS");
  band = EnvelopeBand(ch_lo, ch_hi);
  band->Draw("F SAME");
  curves.clear();
  for (size_t c = 0; c < cals.size(); c++) {
    curves.push_back(ShiftCurve(cals[c], ch_lo, ch_hi));
    curves.back()->Draw("L SAME");
  }
  ref_changes = new TGraph(N_REF);
  for (Int_t i = 0; i < N_REF; i++)
    ref_changes->SetPoint(i, Channel(JK_REF_KEV[i]),
                          1000.0 * (CURRENT_REF_KEV[i] - JK_REF_KEV[i]));
  ref_changes->SetMarkerStyle(20);
  ref_changes->SetMarkerSize(2.0);
  ref_changes->Draw("P SAME");
  Double_t ch_ge = Channel(GE_JK_KEV);
  TLine *ge_line = new TLine(ch_ge, y_lo, ch_ge, y_hi);
  ge_line->SetLineStyle(2);
  ge_line->SetLineColor(kGray + 2);
  ge_line->Draw();
  ge_markers.clear();
  if (!mark_ge)
    return;
  for (size_t c = 0; c < cals.size(); c++) {
    TGraph *m = new TGraph(1);
    m->SetPoint(0, ch_ge, ShiftEV(cals[c].weights, ch_ge));
    m->SetMarkerStyle(cals[c].marker);
    m->SetMarkerColor(cals[c].color);
    m->SetLineColor(cals[c].color);
    m->SetLineWidth(3);
    m->SetMarkerSize(2.2);
    m->Draw("P SAME");
    ge_markers.push_back(m);
  }
}

void PlotRecalibration(const std::vector<Calibration> &cals) {
  std::vector<TString> labels;
  Double_t ch_ge = Channel(GE_JK_KEV);
  for (size_t c = 0; c < cals.size(); c++)
    labels.push_back(cals[c].name +
                     Form(": %+.2f eV", ShiftEV(cals[c].weights, ch_ge)));
  Double_t env_low = 0, env_high = 0;
  LinearEnvelopeEV(ch_ge, env_low, env_high);
  TString band_label =
      Form("Linear, any weights: %+.2f to %+.2f eV", env_low, env_high);
  TString refs_label = "Calibration lines";
  TString jk_label = "J&K's calibration";
  TString ge_label = Form("68.75 keV line, %.3f keV on J&K's", GE_JK_KEV);

  // Full range: the ADC against energy above, today's calibrations below.
  const Double_t ch_lo = Channel(57.5), ch_hi = Channel(77.5);
  const Int_t canvas_w = 1400, plot_h = 800, panel_h = 700;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(canvas_w, plot_h + panel_h);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(canvas, panel_h, upper, lower);
  upper->cd();
  TH1D *frame =
      new TH1D("adc_frame", ";ADC channel;Energy [keV]", 100, ch_lo, ch_hi);
  frame->SetMinimum(57.5);
  frame->SetMaximum(77.5);
  frame->SetStats(kFALSE);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw("AXIS");
  TGraph *jk_line = new TGraph(2);
  jk_line->SetPoint(0, ch_lo, JKEnergyKeV(ch_lo));
  jk_line->SetPoint(1, ch_hi, JKEnergyKeV(ch_hi));
  jk_line->SetLineColor(kBlack);
  jk_line->SetLineWidth(3);
  jk_line->Draw("L SAME");
  TGraph *refs = new TGraph(N_REF);
  for (Int_t i = 0; i < N_REF; i++)
    refs->SetPoint(i, Channel(JK_REF_KEV[i]), CURRENT_REF_KEV[i]);
  refs->SetMarkerStyle(20);
  refs->SetMarkerSize(2.0);
  refs->Draw("P SAME");
  TGraph *ge_point = new TGraph(1);
  ge_point->SetPoint(0, ch_ge, GE_JK_KEV);
  ge_point->SetMarkerStyle(29);
  ge_point->SetMarkerSize(3.0);
  ge_point->SetMarkerColor(kOrange + 7);
  ge_point->Draw("P SAME");
  TLine *ge_line = new TLine(ch_ge, 57.5, ch_ge, 77.5);
  ge_line->SetLineStyle(2);
  ge_line->SetLineColor(kGray + 2);
  ge_line->Draw();
  std::vector<TString> upper_labels;
  upper_labels.push_back(jk_label);
  upper_labels.push_back(refs_label + ", today's energies");
  upper_labels.push_back(ge_label);
  TLegend *upper_leg = FigureStyle::AddFrameLegend(upper, upper_labels, kTRUE);
  upper_leg->AddEntry(jk_line, upper_labels[0], "l");
  upper_leg->AddEntry(refs, upper_labels[1], "p");
  upper_leg->AddEntry(ge_point, upper_labels[2], "p");

  lower->cd();
  std::vector<TGraph *> curves, ge_markers;
  TGraph *band = nullptr, *ref_changes = nullptr;
  DrawShiftPanel(cals, ch_lo, ch_hi, -0.2, 8.0, "shift_frame", kFALSE, curves,
                 band, ref_changes, ge_markers);
  std::vector<TString> lower_labels = labels;
  lower_labels.push_back(band_label);
  lower_labels.push_back(refs_label);
  TLegend *lower_leg = FigureStyle::AddFrameLegend(lower, lower_labels, kFALSE);
  for (size_t c = 0; c < curves.size(); c++)
    lower_leg->AddEntry(curves[c], labels[c], "l");
  lower_leg->AddEntry(band, band_label, "f");
  lower_leg->AddEntry(ref_changes, refs_label, "p");
  PlottingUtils::ScaleFigure(canvas, canvas_w, plot_h);
  PlottingUtils::SaveFigure(canvas, "jones_kraner_recalibration", "",
                            PlotSaveOptions::kLINEAR);

  // Around the line: where each calibration reads it.
  TCanvas *zcanvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  zcanvas->SetCanvasSize(canvas_w, plot_h);
  std::vector<TGraph *> zcurves, zmarkers;
  TGraph *zband = nullptr, *zrefs = nullptr;
  DrawShiftPanel(cals, ch_ge - 50.0, ch_ge + 50.0, 1.50, 2.30, "zoom_frame",
                 kTRUE, zcurves, zband, zrefs, zmarkers);
  std::vector<TString> zoom_labels = labels;
  zoom_labels.push_back(band_label);
  TLegend *zleg =
      FigureStyle::AddFrameLegend(zcanvas, zoom_labels, kTRUE, kFALSE);
  for (size_t c = 0; c < zcurves.size(); c++)
    zleg->AddEntry(zmarkers[c], labels[c], "pl");
  zleg->AddEntry(zband, band_label, "f");
  PlottingUtils::SaveFigure(zcanvas, "jones_kraner_recalibration_zoom", "",
                            PlotSaveOptions::kLINEAR);
}

void PrintAverage(const std::vector<Input> &inputs, Double_t jk_shift_kev,
                  Double_t jk_shift_err_kev, const TString &title,
                  Double_t &mean_kev, Double_t &mean_err_kev) {
  std::vector<Double_t> mus, errs;
  for (size_t i = 0; i < inputs.size(); i++) {
    Double_t mu = inputs[i].mu;
    Double_t err = inputs[i].err;
    if (inputs[i].jones_kraner) {
      mu += jk_shift_kev;
      err = std::sqrt(err * err + jk_shift_err_kev * jk_shift_err_kev);
    }
    mus.push_back(mu);
    errs.push_back(err);
  }
  Double_t sw = 0, swx = 0;
  for (size_t i = 0; i < mus.size(); i++) {
    sw += 1.0 / (errs[i] * errs[i]);
    swx += mus[i] / (errs[i] * errs[i]);
  }
  Double_t mean = swx / sw;
  Double_t chi2 = 0;
  for (size_t i = 0; i < mus.size(); i++)
    chi2 += std::pow((mus[i] - mean) / errs[i], 2);
  std::cout << title << ": " << std::fixed << std::setprecision(4) << mean
            << " +/- " << std::sqrt(1.0 / sw) << " keV  (chi2/ndf "
            << std::setprecision(1) << chi2 << "/" << mus.size() - 1 << ")"
            << std::endl;
  for (size_t i = 0; i < mus.size(); i++)
    std::cout << "  " << std::left << std::setw(34) << inputs[i].name
              << std::right << std::setprecision(4) << mus[i] << " +/- "
              << errs[i] << "  weight " << std::setprecision(3)
              << (1.0 / (errs[i] * errs[i])) / sw << std::endl;
  mean_kev = mean;
  mean_err_kev = std::sqrt(1.0 / sw);
}

struct Row {
  TString label;
  Double_t mu, err; // keV, as plotted
  Int_t marker;
  Bool_t in_mean;
};

// Every measurement against the recommended value: the rows that enter the
// mean first, then the ones left out, top to bottom.
void PlotComparison(const std::vector<Row> &rows, Double_t mean_kev,
                    Double_t mean_err_kev) {
  Int_t n = (Int_t)rows.size();
  // Empty rows above the data hold the legend.
  const Int_t legend_rows = 4;
  const Double_t e_lo = 68.700, e_hi = 68.825;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(1400, 800);
  canvas->SetLeftMargin(0.40);
  canvas->SetGridy(kFALSE);
  TH2D *frame = new TH2D("comparison_frame", "", 100, e_lo, e_hi,
                         n + legend_rows, 0.5, n + legend_rows + 0.5);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(frame, kBlack, ";E_{#gamma} [keV];");
  for (Int_t i = 0; i < n; i++)
    frame->GetYaxis()->SetBinLabel(n - i, rows[i].label);
  for (Int_t k = 1; k <= legend_rows; k++)
    frame->GetYaxis()->SetBinLabel(n + k, "");
  frame->GetYaxis()->SetTickLength(0);
  frame->GetXaxis()->SetNdivisions(505);
  frame->Draw("AXIS");

  TBox *band =
      new TBox(mean_kev - mean_err_kev, 0.5, mean_kev + mean_err_kev, n + 0.5);
  band->SetFillColorAlpha(kAzure + 1, 0.35);
  band->SetLineWidth(0);
  band->Draw();
  TLine *mean_line = new TLine(mean_kev, 0.5, mean_kev, n + 0.5);
  mean_line->SetLineColor(kAzure + 2);
  mean_line->SetLineWidth(PlottingUtils::GetLineWidth());
  mean_line->Draw();

  Int_t n_in_mean = 0;
  for (Int_t i = 0; i < n; i++)
    if (rows[i].in_mean)
      n_in_mean++;
  TLine *divider =
      new TLine(e_lo, n - n_in_mean + 0.5, e_hi, n - n_in_mean + 0.5);
  divider->SetLineStyle(2);
  divider->SetLineColor(kGray + 2);
  divider->Draw();

  TGraphErrors *in_mean = new TGraphErrors();
  TGraphErrors *left_out = new TGraphErrors();
  TGraphErrors *evaluation = new TGraphErrors();
  for (Int_t i = 0; i < n; i++) {
    TGraphErrors *g = left_out;
    if (rows[i].in_mean)
      g = in_mean;
    else if (rows[i].marker == 25)
      g = evaluation;
    Int_t k = g->GetN();
    g->SetPoint(k, rows[i].mu, n - i);
    g->SetPointError(k, rows[i].err, 0);
  }
  PlottingUtils::ConfigureGraph(in_mean, kBlack);
  PlottingUtils::ConfigureGraph(left_out, kBlack);
  PlottingUtils::ConfigureGraph(evaluation, kGray + 2);
  in_mean->SetMarkerStyle(20);
  left_out->SetMarkerStyle(24);
  evaluation->SetMarkerStyle(25);
  in_mean->SetMarkerSize(1.8);
  left_out->SetMarkerSize(1.8);
  evaluation->SetMarkerSize(1.8);
  in_mean->Draw("P SAME");
  left_out->Draw("P SAME");
  evaluation->Draw("P SAME");
  frame->Draw("AXIS SAME");

  std::vector<TString> labels;
  labels.push_back("No Ge detector: in the mean");
  labels.push_back("Ge or unstated detector");
  labels.push_back(Form("Mean: %.4f #pm %.4f keV", mean_kev, mean_err_kev));
  TLegend *leg = FigureStyle::AddFrameLegend(canvas, labels, kFALSE);
  leg->AddEntry(in_mean, labels[0], "pe");
  leg->AddEntry(left_out, labels[1], "pe");
  leg->AddEntry(band, labels[2], "f");
  PlottingUtils::SaveFigure(canvas, "recommended_value_comparison", "",
                            PlotSaveOptions::kLINEAR);
}

} // namespace

void RecommendedValue() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/evaluation",
                                project_root + "/root_files");

  std::vector<Calibration> cals;
  Calibration ka_only = {
      "K#alpha_{2} and K#alpha_{1} only", {0, 1, 1}, kBlue, 21};
  Calibration am_ka2 = {"Am and K#alpha_{2}", {1, 1, 0}, kRed + 1, 22};
  Calibration am_ka1 = {"Am and K#alpha_{1}", {1, 0, 1}, kGreen + 2, 23};
  Calibration all_linear = {"All three, linear", {1, 1, 1}, kBlack, 20};
  Calibration all_quadratic = {
      "All three, quadratic", {0, 0, 0}, kMagenta + 1, 33};
  cals.push_back(ka_only);
  cals.push_back(am_ka2);
  cals.push_back(am_ka1);
  cals.push_back(all_linear);
  cals.push_back(all_quadratic);

  Double_t ch_ge = Channel(GE_JK_KEV);
  std::cout << "Reference lines, today's - Jones and Kraner's:" << std::endl;
  for (Int_t i = 0; i < N_REF; i++)
    std::cout << "  " << REF_NAMES[i] << ": " << std::fixed
              << std::setprecision(2)
              << 1000.0 * (CURRENT_REF_KEV[i] - JK_REF_KEV[i]) << " eV"
              << std::endl;
  std::cout << "The " << std::setprecision(3) << GE_JK_KEV
            << " keV line on today's calibrations:" << std::endl;
  Double_t low = 1e30, high = -1e30;
  for (size_t c = 0; c < cals.size(); c++) {
    Double_t s = ShiftEV(cals[c].weights, ch_ge);
    low = std::min(low, s);
    high = std::max(high, s);
    std::cout << "  " << cals[c].name << ": " << std::setprecision(2) << s
              << " eV" << std::endl;
  }
  Double_t env_low = 0, env_high = 0;
  LinearEnvelopeEV(ch_ge, env_low, env_high);
  low = std::min(low, env_low);
  high = std::max(high, env_high);
  std::cout << "  linear, any weights: " << env_low << " to " << env_high
            << " eV" << std::endl;

  PlotRecalibration(cals);

  // The adjustment is the middle of every calibration's shift, its error the
  // half-range.
  Double_t jk_shift_kev = 0.5 * (low + high) / 1000.0;
  Double_t jk_shift_err_kev = 0.5 * (high - low) / 1000.0;
  std::cout << "Jones and Kraner adjustment: " << std::setprecision(2)
            << jk_shift_kev * 1000.0 << " +/- " << jk_shift_err_kev * 1000.0
            << " eV" << std::endl;

  // This work's CZT (CombineGeResult.cpp); the collaborator's CdTe and Si
  // (their draft's table); Jones and Kraner's 1971 (p,p'gamma) group, taken
  // with a Si(Li) (their Fig. 2), and their 1975 (n,gamma) value, taken with
  // a Si(Li).
  std::vector<Input> inputs;
  Input czt = {"CZT (this work)", 68.7567, 0.0191, kFALSE};
  Input cdte = {"CdTe (collaborator)", 68.7484, 0.0189, kFALSE};
  Input si1 = {"Si-1 (collaborator)", 68.7529, 0.0114, kFALSE};
  Input si2 = {"Si-2 (collaborator)", 68.7555, 0.0128, kFALSE};
  Input jk71 = {"Jones and Kraner 1971 (p,p'g)", 68.7465, 0.0090, kTRUE};
  Input jk75 = {"Jones and Kraner 1975 (n,g)", 68.7550, 0.0050, kTRUE};
  inputs.push_back(czt);
  inputs.push_back(cdte);
  inputs.push_back(si1);
  inputs.push_back(si2);
  inputs.push_back(jk71);
  inputs.push_back(jk75);

  std::cout << std::endl;
  Double_t published_mean = 0, published_err = 0;
  PrintAverage(inputs, 0.0, 0.0, "As published", published_mean, published_err);
  Double_t mean = 0, mean_err = 0;
  PrintAverage(inputs, jk_shift_kev, jk_shift_err_kev,
               "Jones and Kraner on DDEP and SRD 128", mean, mean_err);

  // The figure: the inputs as averaged, then the values left out. Jones and
  // Kraner's other 1971 values do not say which detector took them (Phys.
  // Rev. C 4, 125); they share the inputs' calibration, so they take the same
  // adjustment. Collar et al. (Phys. Rev. D 103, 122003) used a Ge detector.
  // ENSDF's adopted energy is shown for comparison.
  std::vector<Row> rows;
  const char *const labels[6] = {"OSU CZT",
                                 "PSI CdTe",
                                 "PSI Si, run 1",
                                 "PSI Si, run 2",
                                 "Jones & Kraner 1971, (p,p'#gamma)",
                                 "Jones & Kraner 1975, (n,#gamma)"};
  for (size_t i = 0; i < inputs.size(); i++) {
    Row r = {labels[i], inputs[i].mu, inputs[i].err, 20, kTRUE};
    if (inputs[i].jones_kraner) {
      r.mu += jk_shift_kev;
      r.err = std::sqrt(r.err * r.err + jk_shift_err_kev * jk_shift_err_kev);
    }
    rows.push_back(r);
  }
  Double_t jk_alpha_err =
      std::sqrt(0.0127 * 0.0127 + jk_shift_err_kev * jk_shift_err_kev);
  Double_t jk_capture_err =
      std::sqrt(0.0180 * 0.0180 + jk_shift_err_kev * jk_shift_err_kev);
  Row jk71_alpha = {"Jones & Kraner 1971, (#alpha,#alpha'#gamma)",
                    68.7660 + jk_shift_kev, jk_alpha_err, 24, kFALSE};
  Row jk71_capture = {"Jones & Kraner 1971, (n,#gamma)", 68.7460 + jk_shift_kev,
                      jk_capture_err, 24, kFALSE};
  Row collar = {"Collar et al. 2021, Ge", 68.734, 0.020, 24, kFALSE};
  Row ensdf = {"ENSDF", ReferenceLines::GE73M_ENSDF, 0.007, 25, kFALSE};
  rows.push_back(jk71_alpha);
  rows.push_back(jk71_capture);
  rows.push_back(collar);
  rows.push_back(ensdf);
  PlotComparison(rows, mean, mean_err);
}
