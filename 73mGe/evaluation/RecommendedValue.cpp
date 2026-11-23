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
#include <fstream>
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

// The weighted least-squares line y = intercept + slope * x through three
// points; two nonzero weights give the line through those two.
void WeightedLine(const Double_t x[N_REF], const Double_t y[N_REF],
                  const Double_t weights[N_REF], Double_t &intercept,
                  Double_t &slope) {
  Double_t sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (Int_t i = 0; i < N_REF; i++) {
    Double_t w = weights[i];
    sw += w;
    swx += w * x[i];
    swy += w * y[i];
    swxx += w * x[i] * x[i];
    swxy += w * x[i] * y[i];
  }
  slope = (sw * swxy - swx * swy) / (sw * swxx - swx * swx);
  intercept = (swy - slope * swx) / sw;
}

// The ADC calibrated on today's reference energies: a weighted least-squares
// line through the three lines' channels, or, with all weights zero, the
// quadratic through all three.
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
  Double_t channels[N_REF];
  for (Int_t i = 0; i < N_REF; i++)
    channels[i] = Channel(JK_REF_KEV[i]);
  Double_t intercept = 0, slope = 0;
  WeightedLine(channels, CURRENT_REF_KEV, weights, intercept, slope);
  return intercept + slope * channel;
}

// The collaborator's calibration energies: Am-241 at DDEP's 59.5409 keV and
// Bearden's Pb K-alpha2 and K-alpha1, as stated in Phys. Rev. D 110, 083043
// (2024). arXiv:2609.20848 does not state them for the PSI values, but the
// collaborator's fit notebook for a PSI CdTe run (Fit_ThursdayNight.ipynb)
// uses the same Pb energies.
const Double_t BEARDEN_PB_REF_KEV[N_REF] = {59.5409, 72.8042, 74.9694};

// How far today's Pb K-alpha energies move an energy read on a weighted
// least-squares line through BEARDEN_PB_REF_KEV, in eV.
Double_t BeardenPbShiftEV(const Double_t weights[N_REF], Double_t energy_kev) {
  Double_t intercept = 0, slope = 0;
  WeightedLine(BEARDEN_PB_REF_KEV, CURRENT_REF_KEV, weights, intercept, slope);
  return 1000.0 * (intercept + slope * energy_kev - energy_kev);
}

// The range of that shift over lines in which Am-241 weighs at least as much
// as each Pb line: both papers calibrate on it, and it is the best-measured
// line (long source runs; in the PSI Si data its position is twice as precise
// as Pb K-alpha1's).
void BeardenPbEnvelopeEV(Double_t energy_kev, Double_t &low, Double_t &high) {
  low = 1e30;
  high = -1e30;
  for (Int_t a = -16; a <= 0; a++)
    for (Int_t b = -16; b <= 0; b++) {
      Double_t weights[N_REF] = {1.0, std::pow(10.0, a / 4.0),
                                 std::pow(10.0, b / 4.0)};
      Double_t s = BeardenPbShiftEV(weights, energy_kev);
      low = std::min(low, s);
      high = std::max(high, s);
    }
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

// The reference energies a value was calibrated on.
enum ReferenceSet { kTodaysReferences, kJonesKraner, kBeardenPb };

struct Input {
  TString name;     // printed
  TString label;    // in the figure
  Double_t mu, err; // keV
  ReferenceSet references;
};

// The shifts that put each reference set on today's energies, keV.
struct Adjustment {
  Double_t jk = 0, jk_err = 0;
  Double_t bearden_pb = 0, bearden_pb_err = 0;
};

void Adjust(ReferenceSet references, const Adjustment &adjustment, Double_t &mu,
            Double_t &err) {
  Double_t shift = 0, shift_err = 0;
  if (references == kJonesKraner) {
    shift = adjustment.jk;
    shift_err = adjustment.jk_err;
  } else if (references == kBeardenPb) {
    shift = adjustment.bearden_pb;
    shift_err = adjustment.bearden_pb_err;
  }
  mu += shift;
  err = std::sqrt(err * err + shift_err * shift_err);
}

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

// An error-weighted mean with every value independent; Jones and Kraner's two
// were taken four years apart. What the values do share, the reference
// energies and each group's adjustment, is under 0.3 eV and nearly equal
// for all of them: a BLUE with it off the diagonal gives the same mean to
// 0.01 eV.
void PrintAverage(const std::vector<Input> &inputs,
                  const Adjustment &adjustment, const TString &title,
                  Double_t &mean_kev, Double_t &mean_err_kev) {
  std::vector<Double_t> mus, errs;
  for (size_t i = 0; i < inputs.size(); i++) {
    Double_t mu = inputs[i].mu;
    Double_t err = inputs[i].err;
    Adjust(inputs[i].references, adjustment, mu, err);
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
// mean first, then the ones left out, top to bottom. The rows fill the frame,
// so the legend sits beside it.
void PlotComparison(const std::vector<Row> &rows, Double_t mean_kev,
                    Double_t mean_err_kev) {
  Int_t n = (Int_t)rows.size();
  const Double_t e_lo = 68.700, e_hi = 68.900;
  std::vector<TString> labels;
  labels.push_back("No Ge detector: in the mean");
  labels.push_back("Ge or unstated detector");
  labels.push_back(Form("Mean: %.4f #pm %.4f keV", mean_kev, mean_err_kev));
  Double_t legend_w = PlottingUtils::LegendWidthPx(labels);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)legend_w);
  plot->cd();
  plot->SetLeftMargin(0.42);
  plot->SetGridy(kFALSE);
  TH2D *frame =
      new TH2D("comparison_frame", "", 100, e_lo, e_hi, n, 0.5, n + 0.5);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(frame, kBlack, ";E_{#gamma} [keV];");
  for (Int_t i = 0; i < n; i++)
    frame->GetYaxis()->SetBinLabel(n - i, rows[i].label);
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

  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
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

  // Each adjustment is the middle of its range of shifts, its error the
  // half-range.
  Adjustment adjustment;
  adjustment.jk = 0.5 * (low + high) / 1000.0;
  adjustment.jk_err = 0.5 * (high - low) / 1000.0;
  std::cout << "Jones and Kraner adjustment: " << std::setprecision(2)
            << adjustment.jk * 1000.0 << " +/- " << adjustment.jk_err * 1000.0
            << " eV" << std::endl;
  Double_t pb_low = 0, pb_high = 0;
  BeardenPbEnvelopeEV(GE_JK_KEV, pb_low, pb_high);
  adjustment.bearden_pb = 0.5 * (pb_low + pb_high) / 1000.0;
  adjustment.bearden_pb_err = 0.5 * (pb_high - pb_low) / 1000.0;
  const Double_t equal_weights[N_REF] = {1, 1, 1};
  std::cout << "Bearden Pb K-alpha adjustment (Am weighted at least as each Pb "
               "line): "
            << pb_low << " to " << pb_high << " eV, equal weights "
            << BeardenPbShiftEV(equal_weights, GE_JK_KEV) << " eV -> "
            << adjustment.bearden_pb * 1000.0 << " +/- "
            << adjustment.bearden_pb_err * 1000.0 << " eV" << std::endl;

  // This work's CZT (CombineGeResult.cpp); the collaborator's CdTe and Si
  // (arXiv:2609.20848, Table III); Jones and Kraner's 1971 (p,p'gamma) group,
  // taken with a Si(Li) (their Fig. 2), and their 1975 (n,gamma) value, taken
  // with a Si(Li).
  std::vector<Input> inputs;
  // Named by year and detector, as in the paper.
  Input jk71 = {"1971 Jones and Kraner (p,p'g)",
                "1971 Jones & Kraner (p,p'#gamma)", 68.7465, 0.0090,
                kJonesKraner};
  Input jk75 = {"1975 Jones and Kraner (n,g)", "1975 Jones & Kraner (n,#gamma)",
                68.7550, 0.0050, kJonesKraner};
  Input cdte = {"2025 CdTe", "2025 CdTe", 68.7484, 0.0189, kBeardenPb};
  Input si1 = {"2025 Si 1", "2025 Si 1", 68.7529, 0.0114, kBeardenPb};
  Input si2 = {"2025 Si 2", "2025 Si 2", 68.7555, 0.0128, kBeardenPb};
  Input czt = {"2026 CZT", "2026 CZT", 68.7567, 0.0191, kTodaysReferences};
  inputs.push_back(jk71);
  inputs.push_back(jk75);
  inputs.push_back(cdte);
  inputs.push_back(si1);
  inputs.push_back(si2);
  inputs.push_back(czt);

  std::cout << std::endl;
  Adjustment none;
  Double_t published_mean = 0, published_err = 0;
  PrintAverage(inputs, none, "As published", published_mean, published_err);
  Adjustment jk_only = adjustment;
  jk_only.bearden_pb = 0;
  jk_only.bearden_pb_err = 0;
  Double_t jk_only_mean = 0, jk_only_err = 0;
  PrintAverage(inputs, jk_only, "Jones and Kraner adjusted only", jk_only_mean,
               jk_only_err);
  Double_t mean = 0, mean_err = 0;
  PrintAverage(inputs, adjustment, "All on DDEP and SRD 128", mean, mean_err);

  // The ionization from the recoil: the summed gamma-plus-recoil peak of
  // Phys. Rev. D 110, 083043 (68.8158 keV, calibrated on Bearden's Pb
  // K-alpha, so it takes the same adjustment) minus the recommended value.
  // Its other uncertainty, 4.6 eV, is that paper's total on the summed peak.
  // RecoilEnergy.cpp turns it into a quenching factor.
  const Double_t SUM_PEAK_KEV = 68.8158, SUM_PEAK_ERR_KEV = 0.0046;
  Double_t sum_low = 0, sum_high = 0;
  BeardenPbEnvelopeEV(SUM_PEAK_KEV, sum_low, sum_high);
  Double_t sum_shift_kev = 0.5 * (sum_low + sum_high) / 1000.0;
  Double_t sum_shift_err_kev = 0.5 * (sum_high - sum_low) / 1000.0;
  Double_t sum_kev = SUM_PEAK_KEV + sum_shift_kev;
  Double_t sum_err_kev = std::sqrt(SUM_PEAK_ERR_KEV * SUM_PEAK_ERR_KEV +
                                   sum_shift_err_kev * sum_shift_err_kev);
  Double_t ionization_ev = 1000.0 * (sum_kev - mean);
  Double_t ionization_err_ev =
      1000.0 * std::sqrt(sum_err_kev * sum_err_kev + mean_err * mean_err);
  std::cout << std::endl
            << "Summed peak on today's references: " << std::setprecision(4)
            << sum_kev << " keV (" << std::showpos << std::setprecision(2)
            << 1000.0 * sum_shift_kev << std::noshowpos << " +/- "
            << 1000.0 * sum_shift_err_kev << " eV)" << std::endl;
  std::cout << "Ionization: " << std::setprecision(1) << ionization_ev
            << " +/- " << ionization_err_ev << " eV" << std::endl;

  std::ofstream result(
      (project_root + "/results/recommended_value.result").Data());
  result << "# RecommendedValue.cpp: the 68.75 keV gamma ray from measurements "
            "without a Ge detector, and the ionization of the 2024 summed "
            "peak; keV unless noted"
         << std::endl;
  result << std::fixed << std::setprecision(6);
  result << "mean " << mean << std::endl;
  result << "mean_err " << mean_err << std::endl;
  result << "mean_unadjusted " << published_mean << std::endl;
  result << "jk_shift " << adjustment.jk << std::endl;
  result << "jk_shift_err " << adjustment.jk_err << std::endl;
  result << "bearden_pb_shift " << adjustment.bearden_pb << std::endl;
  result << "bearden_pb_shift_err " << adjustment.bearden_pb_err << std::endl;
  result << "sum_peak " << sum_kev << std::endl;
  result << "sum_peak_err " << sum_err_kev << std::endl;
  result << std::setprecision(3);
  result << "ionization_ev " << ionization_ev << std::endl;
  result << "ionization_err_ev " << ionization_err_ev << std::endl;

  // The figure: the inputs as averaged, then the values left out. Jones and
  // Kraner's other 1971 values do not say which detector took them (Phys.
  // Rev. C 4, 125); they share the inputs' calibration, so they take the same
  // adjustment. Collar et al. (Phys. Rev. D 103, 122003) used a Ge detector,
  // Weishaupt and Rabenstein (Z. Phys. A 251, 105) Ge(Li) detectors.
  // ENSDF's adopted energy is shown for comparison.
  std::vector<Row> rows;
  for (size_t i = 0; i < inputs.size(); i++) {
    Row r = {inputs[i].label, inputs[i].mu, inputs[i].err, 20, kTRUE};
    Adjust(inputs[i].references, adjustment, r.mu, r.err);
    rows.push_back(r);
  }
  Row jk71_alpha = {"1971 Jones & Kraner (#alpha,#alpha'#gamma)", 68.7660,
                    0.0127, 24, kFALSE};
  Adjust(kJonesKraner, adjustment, jk71_alpha.mu, jk71_alpha.err);
  Row jk71_capture = {"1971 Jones & Kraner (n,#gamma)", 68.7460, 0.0180, 24,
                      kFALSE};
  Adjust(kJonesKraner, adjustment, jk71_capture.mu, jk71_capture.err);
  Row weishaupt = {"1972 Weishaupt & Rabenstein", 68.84, 0.05, 24, kFALSE};
  Row collar = {"2021 Collar et al.", 68.734, 0.020, 24, kFALSE};
  Row ensdf = {"2019 ENSDF evaluation", ReferenceLines::GE73M_ENSDF, 0.007, 25,
               kFALSE};
  rows.push_back(jk71_alpha);
  rows.push_back(jk71_capture);
  rows.push_back(weishaupt);
  rows.push_back(collar);
  rows.push_back(ensdf);
  PlotComparison(rows, mean, mean_err);
}
