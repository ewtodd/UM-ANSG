// The CZT event selection, drawn for the dissertation: Filter.cpp's
// FilterFile cuts applied event by event to a run's full event tree
// (filtered/<run>.root's bef_tree), which is read and never written.
//   czt_filter_x_position      every interaction's x in the 01/13 Cd-shield
//                              10% signal run, on a log axis, with a dashed
//                              line at each pixel centre;
//   czt_filter_depth_energy    the same run's interactions in depth and
//                              energy, below a dashed line the depths the
//                              selection removes;
//   czt_filter_spectra_signal  that run's accepted and rejected events;
//   czt_filter_spectra_source  the same for the 01/16 Am-241 and Ba-133
//                              source run;
//   czt_filter_xy_position_cu10, _zoom  every interaction's y against x in
//                              the 01/13 Cu-shield 10% signal run, whose map
//                              has every pixel, with the pixel outlines:
//                              the whole detector, and a crystal's corner.
// An energy is the event's energy times its pixel's 558 keV gain where the
// pixel has one, as FilterFile matches them, before any calibration. It also
// prints, per run, the share of interactions the M400 places at a pixel
// centre. Writes plots/czt/figures (PNG and PDF).
//
// Usage: root -b -q 'FilterFigure.cpp+'
#include "Constants.hpp"
#include "FigureStyle.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"

#include <TArrayI.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <vector>

namespace {

// Spectrum bins in the M400's on-board keV, and the window, whose upper edge
// is on their grid.
const Double_t BIN_KEV = 0.261;
const Double_t E_LO_KEV = 40.0;
const Int_t N_E_BINS = (Int_t)((110.0 - E_LO_KEV) / BIN_KEV);
const Double_t E_HI_KEV = E_LO_KEV + N_E_BINS * BIN_KEV;

// Position histograms. The map of y against x has coarser bins, which keep
// the pixel centres visible at a figure's size.
const Int_t N_X_BINS = 1000;
const Int_t N_XY_BINS = 500;
const Double_t X_RANGE_MM = 22.0;
// The map's close-up: one corner of a crystal, past its outermost pixels,
// starting just short of -22 mm so the two axes' first labels do not meet.
const Double_t ZOOM_LO_MM = -21.9, ZOOM_HI_MM = -15.0;
const Int_t N_Z_BINS = 500;
const Double_t Z_MAX_MM = 10.0;

const Double_t GE_LINE_KEV = 68.752;

std::vector<std::vector<Float_t>> g_pixel_gain;

// Filter.cpp's pixel gain table, the 558 keV gains of PixelCalibration.cpp.
Bool_t LoadPixelGains() {
  TFile *file = IO::OpenForReading("pixel_calibration.root");
  TH2D *gain = (file && !file->IsZombie())
                   ? static_cast<TH2D *>(file->Get("pixel_gain"))
                   : nullptr;
  const Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  const Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();
  if (!gain || gain->GetNbinsX() != nx || gain->GetNbinsY() != ny) {
    std::cerr << "ERROR: no usable pixel_gain in pixel_calibration.root"
              << std::endl;
    delete file;
    return kFALSE;
  }
  g_pixel_gain.assign(nx, std::vector<Float_t>(ny, 0.0f));
  for (Int_t ix = 0; ix < nx; ix++)
    for (Int_t iy = 0; iy < ny; iy++)
      g_pixel_gain[ix][iy] = (Float_t)gain->GetBinContent(ix + 1, iy + 1);
  file->Close();
  delete file;
  return kTRUE;
}

// Filter.cpp's PixelGain, GetPixelIndex and GetCrystalIndex.
Float_t PixelGain(Int_t ix, Int_t iy) {
  if (ix < 0 || iy < 0 || ix >= (Int_t)g_pixel_gain.size() ||
      iy >= (Int_t)g_pixel_gain[0].size())
    return 0.0f;
  return g_pixel_gain[ix][iy];
}

Int_t GetPixelIndex(Float_t pos, const std::vector<Float_t> &centers) {
  for (Int_t i = 0; i < (Int_t)centers.size(); i++) {
    if (TMath::Abs(pos - centers[i]) <= Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      return i;
  }
  return -1;
}

// The index of the pixel centre nearest pos: the pixel whose area holds it,
// for interactions the M400 places between centres.
Int_t NearestCentre(Float_t pos, const std::vector<Float_t> &centres) {
  Int_t best = 0;
  for (Int_t i = 1; i < (Int_t)centres.size(); i++)
    if (TMath::Abs(pos - centres[i]) < TMath::Abs(pos - centres[best]))
      best = i;
  return best;
}

Int_t GetCrystalIndex(Float_t x, Float_t y) {
  if (x < 0 && y < 0)
    return 0;
  if (x > 0 && y < 0)
    return 1;
  if (x < 0 && y > 0)
    return 2;
  if (x > 0 && y > 0)
    return 3;
  return -1;
}

// How many interactions the M400 places within the selection's band about a
// pixel centre: in x, in y and in both.
struct CentreCounts {
  Long64_t total = 0, in_x = 0, in_y = 0, in_both = 0;
  void Add(Bool_t at_x, Bool_t at_y) {
    total++;
    if (at_x)
      in_x++;
    if (at_y)
      in_y++;
    if (at_x && at_y)
      in_both++;
  }
};

struct SelectionHists {
  TH1D *x = nullptr;
  TH2D *xy = nullptr;
  TH2D *depth_energy = nullptr;
  TH1D *accepted = nullptr, *rejected = nullptr;
  // The interactions that pass every cut but the pixel centre's, in a pixel
  // with a gain (the nearest centre's for those between centres), and of
  // them the ones at a pixel centre in x and y, which the selection keeps:
  // what leaving out charge sharing does.
  TH1D *all_events = nullptr, *single_pixel = nullptr;
  // Every interaction, and the single interactions within 2 mm of the
  // cathode.
  CentreCounts all, single_near_cathode;
};

// One run's interactions through FilterFile's cuts, in its order: the
// crystal, one interaction, the depth, a pixel centre with a gain, and
// pile-up judged against the previous event in the stream.
SelectionHists FillSelection(const TString &run, const TString &tag) {
  SelectionHists h;
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  TTree *tree = (file && !file->IsZombie())
                    ? static_cast<TTree *>(file->Get("bef_tree"))
                    : nullptr;
  if (!tree) {
    std::cerr << "ERROR: no bef_tree in filtered/" << run << ".root"
              << std::endl;
    delete file;
    return h;
  }
  Float_t energy = 0, x = 0, y = 0, z = 0;
  UInt_t eventTime = 0;
  Int_t liveTime = 0, nInteractions = 0, interaction = 0;
  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("zmm", &z);
  tree->SetBranchAddress("eventTime", &eventTime);
  tree->SetBranchAddress("liveTime", &liveTime);
  tree->SetBranchAddress("nInteractions", &nInteractions);
  tree->SetBranchAddress("interaction", &interaction);

  const TString spectrum_title = TString::Format(
      ";Uncorrected Energy [keV];Counts / %.0f eV", 1000.0 * BIN_KEV);
  h.x = new TH1D("x_" + tag,
                 TString::Format(";X Position [mm];Counts / %.0f #mum",
                                 1000.0 * 2.0 * X_RANGE_MM / N_X_BINS),
                 N_X_BINS, -X_RANGE_MM, X_RANGE_MM);
  h.xy = new TH2D("xy_" + tag, ";X Position [mm];Y Position [mm];Counts",
                  N_XY_BINS, -X_RANGE_MM, X_RANGE_MM, N_XY_BINS, -X_RANGE_MM,
                  X_RANGE_MM);
  h.depth_energy = new TH2D(
      "depth_energy_" + tag, ";Uncorrected Energy [keV];Z Position [mm];Counts",
      N_E_BINS, E_LO_KEV, E_HI_KEV, N_Z_BINS, 0.0, Z_MAX_MM);
  h.accepted =
      new TH1D("accepted_" + tag, spectrum_title, N_E_BINS, E_LO_KEV, E_HI_KEV);
  h.rejected =
      new TH1D("rejected_" + tag, spectrum_title, N_E_BINS, E_LO_KEV, E_HI_KEV);
  h.all_events = new TH1D("all_events_" + tag, spectrum_title, N_E_BINS,
                          E_LO_KEV, E_HI_KEV);
  h.single_pixel = new TH1D("single_pixel_" + tag, spectrum_title, N_E_BINS,
                            E_LO_KEV, E_HI_KEV);
  TH1 *hists[7] = {h.x,        h.xy,         h.depth_energy, h.accepted,
                   h.rejected, h.all_events, h.single_pixel};
  for (Int_t k = 0; k < 7; k++)
    hists[k]->SetDirectory(nullptr);

  Long64_t previous_event_time = -1;
  Bool_t same_ms_as_previous = kFALSE;
  const Long64_t n_entries = tree->GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    if (interaction == 0) {
      same_ms_as_previous = (previous_event_time == (Long64_t)eventTime);
      previous_event_time = eventTime;
    }
    const Int_t ix = GetPixelIndex(x, Constants::PIXEL_CENTERS_X_MM);
    const Int_t iy = GetPixelIndex(y, Constants::PIXEL_CENTERS_Y_MM);
    const Float_t gain = PixelGain(ix, iy);
    const Float_t live_us = liveTime * Constants::TENS_OF_NS_TO_S * 1e6;
    const Bool_t piled_up = live_us < Constants::PILEUP_LIVETIME_THRESHOLD_US &&
                            same_ms_as_previous;
    // Every cut but the pixel centre's.
    const Bool_t other_cuts = GetCrystalIndex(x, y) >= 0 &&
                              nInteractions == 1 &&
                              z >= Constants::FILTER_DEPTH_MM && !piled_up;
    const Bool_t accepted = other_cuts && ix >= 0 && iy >= 0 && gain > 0;
    const Float_t nearest_gain =
        PixelGain(NearestCentre(x, Constants::PIXEL_CENTERS_X_MM),
                  NearestCentre(y, Constants::PIXEL_CENTERS_Y_MM));
    if (other_cuts && nearest_gain > 0) {
      const Float_t matched = energy * nearest_gain;
      if (matched > E_LO_KEV && matched < E_HI_KEV) {
        h.all_events->Fill(matched);
        if (ix >= 0 && iy >= 0)
          h.single_pixel->Fill(matched);
      }
    }

    h.x->Fill(x);
    h.xy->Fill(x, y);
    h.all.Add(ix >= 0, iy >= 0);
    if (nInteractions == 1 && z >= Constants::FILTER_DEPTH_MM)
      h.single_near_cathode.Add(ix >= 0, iy >= 0);
    const Float_t shown_energy = (gain > 0) ? energy * gain : energy;
    if (shown_energy <= E_LO_KEV || shown_energy >= E_HI_KEV)
      continue;
    if (z > 0)
      h.depth_energy->Fill(shown_energy, z);
    if (accepted)
      h.accepted->Fill(shown_energy);
    else
      h.rejected->Fill(shown_energy);
  }
  file->Close();
  delete file;
  std::cout << "FILTER-FIGURE " << run << ": " << n_entries << " interactions, "
            << (Long64_t)h.accepted->GetEntries() << " accepted and "
            << (Long64_t)h.rejected->GetEntries() << " rejected in " << E_LO_KEV
            << "-" << E_HI_KEV << " keV" << std::endl;
  const CentreCounts *counts[2] = {&h.all, &h.single_near_cathode};
  const char *names[2] = {"all interactions",
                          "single interactions within 2 mm of the cathode"};
  for (Int_t k = 0; k < 2; k++)
    std::cout << "FILTER-FIGURE " << run << ", " << names[k] << " ("
              << counts[k]->total << "): at a pixel centre in x " << std::fixed
              << std::setprecision(1)
              << 100.0 * counts[k]->in_x / counts[k]->total << "%, in y "
              << 100.0 * counts[k]->in_y / counts[k]->total << "%, in both "
              << 100.0 * counts[k]->in_both / counts[k]->total << "%"
              << std::defaultfloat << std::setprecision(6) << std::endl;
  return h;
}

// The y title beside labels as wide as the axis maximum's.
// Above five digits ROOT labels the axis in thousands (or millions) and puts
// the power of ten above it, so the labels are those of the scaled maximum.
void FitTitleToMaximum(TVirtualPad *pad, TH1 *frame, Double_t maximum) {
  Double_t labelled = maximum;
  while (labelled >= 1e5)
    labelled /= 1e3;
  std::vector<TString> labels;
  labels.push_back(TString::Format("%.0f", labelled));
  FigureStyle::FitTitleToLabels(pad, frame, labels);
}

void DrawXPosition(const SelectionHists &h) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  PlottingUtils::ConfigureHistogram(h.x, kBlack);
  h.x->SetFillStyle(0);
  h.x->SetMinimum(5.0);
  h.x->SetMaximum(3.0 * h.x->GetMaximum());
  h.x->Draw("HIST");
  canvas->Update();
  // A line at each pixel centre, about which the selection keeps a band of
  // PIXEL_ACCEPT_HALFWIDTH_MM, narrower than a bin. Dashed and drawn over the
  // histogram, whose spike at each centre would otherwise hide it.
  const Double_t y_lo = h.x->GetMinimum(), y_hi = h.x->GetMaximum();
  for (size_t i = 0; i < Constants::PIXEL_CENTERS_X_MM.size(); i++) {
    const Double_t c = Constants::PIXEL_CENTERS_X_MM[i];
    TLine *centre = new TLine(c, y_lo, c, y_hi);
    centre->SetLineColor(kRed + 1);
    centre->SetLineStyle(2);
    centre->SetLineWidth(PlottingUtils::GetLineWidth());
    centre->Draw();
  }
  PlottingUtils::SaveFigure(canvas, "czt_filter_x_position", "",
                            PlotSaveOptions::kLOG);
}

// The colour axis in the x and y axes' sizes, with room for its labels and
// title.
void StyleColourAxis(TH2 *h, TCanvas *canvas) {
  TAxis *z_axis = h->GetZaxis();
  const TAxis *y_axis = h->GetYaxis();
  z_axis->SetLabelFont(y_axis->GetLabelFont());
  z_axis->SetLabelSize(y_axis->GetLabelSize());
  z_axis->SetTitleFont(y_axis->GetTitleFont());
  z_axis->SetTitleSize(y_axis->GetTitleSize());
  z_axis->SetTitleOffset(1.25);
  canvas->SetRightMargin(0.22);
}

// One crystal's pixel edges along an axis from its PIXELS_PER_CRYSTAL centres
// starting at first: halfway between neighbouring centres, and half a pitch
// beyond the outermost ones.
const Int_t PIXELS_PER_CRYSTAL = 11;

std::vector<Double_t> PixelEdges(const std::vector<Float_t> &centres,
                                 Int_t first) {
  const Int_t last = first + PIXELS_PER_CRYSTAL - 1;
  std::vector<Double_t> edges;
  edges.push_back(centres[first] - 0.5 * (centres[first + 1] - centres[first]));
  for (Int_t i = first; i < last; i++)
    edges.push_back(0.5 * (centres[i] + centres[i + 1]));
  edges.push_back(centres[last] + 0.5 * (centres[last] - centres[last - 1]));
  return edges;
}

// Whether a position lies between a crystal's outermost pixel centres along
// an axis, the only positions the M400 reports.
Bool_t InReportedRange(Double_t pos, const std::vector<Float_t> &centres) {
  for (Int_t first = 0; first < (Int_t)centres.size();
       first += PIXELS_PER_CRYSTAL)
    if (pos >= centres[first] && pos <= centres[first + PIXELS_PER_CRYSTAL - 1])
      return kTRUE;
  return kFALSE;
}

// A red line from (x1, y1) to (x2, y2), clipped to the square [lo, hi] the
// map shows, or nothing if it lies outside.
void DrawOutline(Double_t x1, Double_t y1, Double_t x2, Double_t y2,
                 Double_t lo, Double_t hi) {
  x1 = TMath::Max(x1, lo);
  y1 = TMath::Max(y1, lo);
  x2 = TMath::Min(x2, hi);
  y2 = TMath::Min(y2, hi);
  if (x1 > x2 || y1 > y2)
    return;
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(kRed);
  // Twice the standard width, to stay visible over the map in print.
  line->SetLineWidth(2 * PlottingUtils::GetLineWidth());
  line->Draw();
}

// Every pixel's outline, crystal by crystal, within the square [lo, hi].
void DrawPixelOutlines(Double_t lo, Double_t hi) {
  const std::vector<Float_t> &cx = Constants::PIXEL_CENTERS_X_MM;
  const std::vector<Float_t> &cy = Constants::PIXEL_CENTERS_Y_MM;
  for (Int_t fx = 0; fx < (Int_t)cx.size(); fx += PIXELS_PER_CRYSTAL)
    for (Int_t fy = 0; fy < (Int_t)cy.size(); fy += PIXELS_PER_CRYSTAL) {
      const std::vector<Double_t> ex = PixelEdges(cx, fx);
      const std::vector<Double_t> ey = PixelEdges(cy, fy);
      for (size_t i = 0; i < ex.size(); i++)
        if (ex[i] >= lo && ex[i] <= hi)
          DrawOutline(ex[i], ey.front(), ex[i], ey.back(), lo, hi);
      for (size_t i = 0; i < ey.size(); i++)
        if (ey[i] >= lo && ey[i] <= hi)
          DrawOutline(ex.front(), ey[i], ex.back(), ey[i], lo, hi);
    }
}

// Every interaction's y against its x over the square [lo, hi], on a log
// colour scale in Cividis, with each pixel's outline in red, so that the
// interactions placed at pixel centres show in the middle of the pixels. The
// canvas is shaped so that the frame is square, and the empty bins where
// the M400 reports positions take the scale's lowest colour.
const Double_t EMPTY_BIN_LEVEL = 0.5;

void DrawXYPosition(const SelectionHists &h, const TString &name, Double_t lo,
                    Double_t hi) {
  TH2D *shown = static_cast<TH2D *>(h.xy->Clone(name));
  shown->SetDirectory(nullptr);
  for (Int_t bx = 1; bx <= shown->GetNbinsX(); bx++)
    for (Int_t by = 1; by <= shown->GetNbinsY(); by++)
      if (shown->GetBinContent(bx, by) <= 0 &&
          InReportedRange(shown->GetXaxis()->GetBinCenter(bx),
                          Constants::PIXEL_CENTERS_X_MM) &&
          InReportedRange(shown->GetYaxis()->GetBinCenter(by),
                          Constants::PIXEL_CENTERS_Y_MM))
        shown->SetBinContent(bx, by, EMPTY_BIN_LEVEL);
  // The palette every other figure uses, put back afterwards.
  TArrayI default_palette = TColor::GetPalette();
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetCanvasSize(1100, 900);
  canvas->SetGrid(0, 0);
  PlottingUtils::Configure2DHistogram(shown, canvas);
  gStyle->SetPalette(kCividis);
  canvas->SetLogz(kTRUE);
  shown->SetMinimum(EMPTY_BIN_LEVEL);
  StyleColourAxis(shown, canvas);
  // Room for labels like -20 and the y title beside them.
  canvas->SetLeftMargin(0.18);
  shown->GetXaxis()->SetRangeUser(lo, hi);
  shown->GetYaxis()->SetRangeUser(lo, hi);
  shown->Draw("COLZ");
  std::vector<TString> widest_label;
  widest_label.push_back("#minus20");
  FigureStyle::FitTitleToLabels(canvas, shown, widest_label);
  canvas->Update();
  DrawPixelOutlines(lo, hi);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
  gStyle->SetPalette(default_palette.GetSize(), default_palette.GetArray());
}

void DrawDepthEnergy(const SelectionHists &h) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::Configure2DHistogram(h.depth_energy, canvas);
  canvas->SetLogz(kFALSE);
  StyleColourAxis(h.depth_energy, canvas);
  h.depth_energy->Draw("COLZ");
  canvas->Update();
  // The depths the selection removes, more than 2 mm from the cathode, lie
  // below a dashed line at the cut.
  TLine *cut = new TLine(E_LO_KEV, Constants::FILTER_DEPTH_MM, E_HI_KEV,
                         Constants::FILTER_DEPTH_MM);
  cut->SetLineColor(kWhite);
  cut->SetLineStyle(2);
  cut->SetLineWidth(PlottingUtils::GetLineWidth());
  cut->Draw();
  PlottingUtils::SaveFigure(canvas, "czt_filter_depth_energy", "",
                            PlotSaveOptions::kLINEAR);
}

// Two spectra of one run, red over black, with a legend in the frame's top
// right and, for a signal run, the Ge line marked.
void DrawSpectrumPair(TH1D *red, const TString &red_label, TH1D *black,
                      const TString &black_label, const TString &name,
                      Bool_t mark_ge) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  // Room for five-digit labels and the title beside them.
  canvas->SetLeftMargin(0.17);
  PlottingUtils::ConfigureHistogram(red, kRed + 1);
  PlottingUtils::ConfigureHistogram(black, kBlack);
  red->SetFillStyle(0);
  black->SetFillStyle(0);
  const Double_t y_max =
      1.15 * TMath::Max(red->GetMaximum(), black->GetMaximum());
  red->SetMinimum(0.0);
  red->SetMaximum(y_max);
  // No 0 on the y axis: it would meet the x axis's first label in the corner.
  red->GetYaxis()->ChangeLabel(1, -1, 0.0);
  red->Draw("HIST");
  black->Draw("HIST SAME");
  red->Draw("HIST SAME");
  if (mark_ge) {
    TLine *ge = new TLine(GE_LINE_KEV, 0.0, GE_LINE_KEV, y_max);
    ge->SetLineColor(kViolet);
    ge->SetLineStyle(2);
    ge->SetLineWidth(PlottingUtils::GetLineWidth());
    ge->Draw();
    TLatex *label = new TLatex(GE_LINE_KEV - 0.8, 0.93 * y_max,
                               "^{73}Ge 68.75 keV #gamma ray");
    label->SetTextFont(43);
    label->SetTextSize(38); // AddText()'s annotation size
    label->SetTextColor(kViolet);
    label->SetTextAlign(32);
    label->Draw();
  }
  std::vector<TString> labels;
  labels.push_back(red_label);
  labels.push_back(black_label);
  TLegend *legend = FigureStyle::AddFrameLegend(canvas, labels, kFALSE);
  legend->AddEntry(red, labels[0], "l");
  legend->AddEntry(black, labels[1], "l");
  legend->Draw();
  FitTitleToMaximum(canvas, red, y_max);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

void DrawSpectra(const SelectionHists &h, const TString &name, Bool_t mark_ge) {
  DrawSpectrumPair(h.accepted, "Accepted", h.rejected, "Rejected", name,
                   mark_ge);
}

// What leaving out charge sharing costs in the spectra's window: of the
// events that pass every other cut, the share in a single pixel.
void PrintSharingCost(const SelectionHists &h, const TString &run) {
  std::cout << "FILTER-FIGURE " << run << ", " << E_LO_KEV << "-" << E_HI_KEV
            << " keV, every cut but the pixel centre's: "
            << (Long64_t)h.all_events->GetEntries() << " events, of which "
            << (Long64_t)h.single_pixel->GetEntries() << " in a single pixel ("
            << std::fixed << std::setprecision(1)
            << 100.0 * h.single_pixel->GetEntries() / h.all_events->GetEntries()
            << "%)" << std::defaultfloat << std::setprecision(6) << std::endl;
}

// ------------------------------------------------ charge sharing at 558 keV

// The 558.456 keV line the pixel gains are measured on (PixelCalibration.cpp,
// the 01/16 graphite-castle run), in the interactions the M400 places at a
// pixel centre and in those it places between centres, whose charge it shares
// between pixels. The interactions are taken as the gain fits take them:
// every interaction, with none of the selection's other cuts, at a centre
// when within PIXEL_ACCEPT_HALFWIDTH_MM of it in x and y. Each is matched with
// its nearest pixel's gain; pixels without one are left out.
const Double_t CD558_KEV = 558.456;
const Double_t GAIN_LINE_LO_KEV = 545.0;
const Int_t N_GAIN_LINE_BINS = (Int_t)(25.0 / BIN_KEV);
const Double_t GAIN_LINE_HI_KEV = GAIN_LINE_LO_KEV + N_GAIN_LINE_BINS * BIN_KEV;
// Where each line is fitted, a Gaussian on a linear continuum.
const Double_t GAIN_FIT_LO_KEV = 545.0, GAIN_FIT_HI_KEV = 569.5;

struct LineFit {
  Bool_t ok = kFALSE;
  Double_t mu = 0, mu_err = 0, sigma = 0, peak_over_continuum = 0;
};

LineFit FitGainLine(TH1D *h) {
  LineFit out;
  TF1 f(TString(h->GetName()) + "_fit", "gaus(0) + pol1(3)", GAIN_FIT_LO_KEV,
        GAIN_FIT_HI_KEV);
  f.SetParameters(h->GetMaximum(), CD558_KEV, 2.0,
                  h->GetBinContent(h->FindBin(GAIN_FIT_LO_KEV)), 0.0);
  TFitResultPtr result = h->Fit(&f, "QRNS");
  out.ok = result.Get() && result->IsValid();
  out.mu = f.GetParameter(1);
  out.mu_err = f.GetParError(1);
  out.sigma = TMath::Abs(f.GetParameter(2));
  out.peak_over_continuum =
      f.GetParameter(0) / (f.GetParameter(3) + f.GetParameter(4) * out.mu);
  return out;
}

void DrawGainLineSharing() {
  const TString run =
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116;
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  TTree *tree = (file && !file->IsZombie())
                    ? static_cast<TTree *>(file->Get("bef_tree"))
                    : nullptr;
  if (!tree) {
    std::cerr << "ERROR: no bef_tree in filtered/" << run << ".root"
              << std::endl;
    delete file;
    return;
  }
  Float_t energy = 0, x = 0, y = 0;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("energykeV", 1);
  tree->SetBranchStatus("xmm", 1);
  tree->SetBranchStatus("ymm", 1);
  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  const TString title = TString::Format(
      ";Uncorrected Energy [keV];Counts / %.0f eV", 1000.0 * BIN_KEV);
  TH1D *single = new TH1D("gain_line_single", title, N_GAIN_LINE_BINS,
                          GAIN_LINE_LO_KEV, GAIN_LINE_HI_KEV);
  TH1D *shared = new TH1D("gain_line_shared", title, N_GAIN_LINE_BINS,
                          GAIN_LINE_LO_KEV, GAIN_LINE_HI_KEV);
  single->SetDirectory(nullptr);
  shared->SetDirectory(nullptr);
  const std::vector<Float_t> &cx = Constants::PIXEL_CENTERS_X_MM;
  const std::vector<Float_t> &cy = Constants::PIXEL_CENTERS_Y_MM;
  const Long64_t n_entries = tree->GetEntries();
  for (Long64_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    const Int_t ix = NearestCentre(x, cx), iy = NearestCentre(y, cy);
    const Float_t gain = PixelGain(ix, iy);
    if (!(gain > 0))
      continue;
    const Float_t matched = energy * gain;
    if (matched < GAIN_LINE_LO_KEV || matched >= GAIN_LINE_HI_KEV)
      continue;
    const Bool_t at_centre =
        TMath::Abs(x - cx[ix]) < Constants::PIXEL_ACCEPT_HALFWIDTH_MM &&
        TMath::Abs(y - cy[iy]) < Constants::PIXEL_ACCEPT_HALFWIDTH_MM;
    if (at_centre)
      single->Fill(matched);
    else
      shared->Fill(matched);
  }
  file->Close();
  delete file;

  // Every interaction: what the gain fits would see without the cut.
  TH1D *all = static_cast<TH1D *>(single->Clone("gain_line_all"));
  all->SetDirectory(nullptr);
  all->Add(shared);
  TH1D *hists[3] = {single, shared, all};
  const char *names[3] = {"single pixel", "shared", "all"};
  const Double_t total = all->GetEntries();
  for (Int_t k = 0; k < 3; k++) {
    const LineFit fit = FitGainLine(hists[k]);
    std::cout << "FILTER-FIGURE " << run << ", 558 keV line, " << names[k]
              << ": " << (Long64_t)hists[k]->GetEntries() << " interactions ("
              << std::fixed << std::setprecision(1)
              << 100.0 * hists[k]->GetEntries() / total << "% of "
              << GAIN_LINE_LO_KEV << "-" << GAIN_LINE_HI_KEV << " keV), peak "
              << std::setprecision(3) << fit.mu << " +/- " << fit.mu_err
              << " keV, sigma " << fit.sigma << " keV, peak/continuum "
              << std::setprecision(1) << fit.peak_over_continuum
              << (fit.ok ? "" : " [fit failed]") << std::defaultfloat
              << std::setprecision(6) << std::endl;
  }
  DrawSpectrumPair(single, "Single pixel", all, "All events",
                   "czt_gain_sharing_558", kFALSE);
}

} // namespace

void FilterFigure() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/czt/figures",
                                project_root + "/root_files");
  if (!LoadPixelGains())
    return;
  const SelectionHists signal =
      FillSelection(Constants::CDSHIELDSIGNAL_10PERCENT_20260113, "signal");
  const SelectionHists source =
      FillSelection(Constants::POSTREACTOR_AM241_BA133_20260116, "source");
  // The Cu-shield signal run of the same day, whose map has every pixel: in
  // the Cd-shield 10% runs, and partly the 25% signal run, some pixels report
  // no interactions at all.
  const SelectionHists cu_signal =
      FillSelection(Constants::CUSHIELDSIGNAL_10PERCENT_20260113, "cu_signal");
  if (!signal.x || !source.x || !cu_signal.x)
    return;
  DrawXPosition(signal);
  DrawXYPosition(cu_signal, "czt_filter_xy_position_cu10", -X_RANGE_MM,
                 X_RANGE_MM);
  DrawXYPosition(cu_signal, "czt_filter_xy_position_cu10_zoom", ZOOM_LO_MM,
                 ZOOM_HI_MM);
  DrawDepthEnergy(signal);
  DrawSpectra(signal, "czt_filter_spectra_signal", kTRUE);
  DrawSpectra(source, "czt_filter_spectra_source", kFALSE);
  PrintSharingCost(signal, Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  PrintSharingCost(cu_signal, Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
  PrintSharingCost(source, Constants::POSTREACTOR_AM241_BA133_20260116);
  DrawGainLineSharing();
}
