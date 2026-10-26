// Where in the CZT the 95.902 keV 113Cd(n,g) line is detected, per shielding
// configuration, to test the escape simulation's two pictures (geant4/):
// captures inside the crystal, within a neutron mean free path of the faces
// the neutrons enter, whose conversion deposits stay at the capture, against
// photons from captures in a Cd shield, which convert within a photon
// attenuation length (about 0.5 mm at 96 keV) of the faces they reach.
//
// Single-interaction events on a pixel centre, gain-corrected with the pixel
// table and past the pile-up cut, as Filter.cpp selects them. Each line is a
// window less the continuum of side bands either side, scaled by width. For
// each configuration (its signal and background runs together):
//   - the depth profile (z = 10 mm is the cathode, the face the sources and
//     the Ge target looked at) of the 95.9 keV line against Pb K-beta from
//     the shielding, photons entering from outside at a similar energy in the
//     same runs; and the Am-241 59.54 keV source line, photons entering the
//     cathode;
//   - the lateral profile: the line's rate per pixel against the pixel's
//     ring counted from its crystal's edge (0 = outermost), over the side
//     bands' rate in the same ring, so pixel acceptance cancels; normalised
//     to the central rings 3-5;
//   - the spectrum with the windows, to check them.
// Depth is reconstructed only above about 0.7 mm: the anode face itself is
// not seen.
//
// Usage: root -b -q 'Capture96Positions.cpp+'
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"

#include <TBox.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPad.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <vector>

namespace Capture96 {

// Energies are binned no finer than the ADC channel (Constants::BIN_WIDTH_KEV,
// 255 eV), on a grid from 40 keV; every window edge sits on the grid, so a
// window holds whole channel-width bins.
const Double_t GRID_START_KEV = 40.0;
const Int_t N_GRID_BINS = 274; // to 109.87 keV

Double_t OnGrid(Double_t e) {
  return GRID_START_KEV +
         TMath::Nint((e - GRID_START_KEV) / Constants::BIN_WIDTH_KEV) *
             Constants::BIN_WIDTH_KEV;
}

// A line's window and the side bands either side of it, keV.
struct LineWindow {
  const char *name;
  Double_t lo, hi, s1_lo, s1_hi, s2_lo, s2_hi;
  Double_t SideScale() const {
    return (hi - lo) / ((s1_hi - s1_lo) + (s2_hi - s2_lo));
  }
};

LineWindow OnGridWindow(const char *name, Double_t lo, Double_t hi,
                        Double_t s1_lo, Double_t s1_hi, Double_t s2_lo,
                        Double_t s2_hi) {
  LineWindow w = {name,          OnGrid(lo),    OnGrid(hi),   OnGrid(s1_lo),
                  OnGrid(s1_hi), OnGrid(s2_lo), OnGrid(s2_hi)};
  return w;
}

const LineWindow LINE96 =
    OnGridWindow("95.9 keV line", 94.9, 96.9, 91.5, 93.5, 98.3, 100.3);
const LineWindow PB_KBETA =
    OnGridWindow("Pb K#beta", 84.0, 88.0, 80.5, 82.5, 89.5, 91.5);
const LineWindow AM59 =
    OnGridWindow("Am-241 59.54 keV", 58.6, 60.5, 55.5, 57.4, 61.0, 62.9);

const Int_t N_Z_BINS = 40;
const Int_t PIXELS_PER_CRYSTAL = 11, N_RINGS = 6;

struct Configuration {
  const char *tag, *label;
  std::vector<TString> runs;
};

struct Profile {
  TH1D *z_line, *z_side, *ring_line, *ring_side;
};

std::vector<std::vector<Float_t>> g_pixel_gain;

// The pixel gain table of PixelCalibration.cpp, as Filter.cpp reads it.
Bool_t LoadPixelGain() {
  TFile *file = IO::OpenForReading("pixel_calibration.root");
  if (!file)
    return kFALSE;
  TH2D *gain = static_cast<TH2D *>(file->Get("pixel_gain"));
  if (!gain)
    return kFALSE;
  const Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  const Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();
  g_pixel_gain.assign(nx, std::vector<Float_t>(ny, 0.0f));
  for (Int_t ix = 0; ix < nx; ix++)
    for (Int_t iy = 0; iy < ny; iy++)
      g_pixel_gain[ix][iy] = gain->GetBinContent(ix + 1, iy + 1);
  file->Close();
  return kTRUE;
}

Int_t PixelIndex(Float_t pos, const std::vector<Float_t> &centers) {
  for (Int_t i = 0; i < (Int_t)centers.size(); i++)
    if (TMath::Abs(pos - centers[i]) <= Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      return i;
  return -1;
}

// The pixel's ring from its crystal's edge: 0 for the outermost pixels.
Int_t Ring(Int_t ix, Int_t iy) {
  Int_t lx = ix % PIXELS_PER_CRYSTAL, ly = iy % PIXELS_PER_CRYSTAL;
  Int_t last = PIXELS_PER_CRYSTAL - 1;
  return TMath::Min(TMath::Min(lx, last - lx), TMath::Min(ly, last - ly));
}

Bool_t InWindow(Double_t e, Double_t lo, Double_t hi) {
  return e >= lo && e < hi;
}

TH1D *NewHist(const TString &name, Int_t n, Double_t lo, Double_t hi) {
  TH1D *h = new TH1D(name, "", n, lo, hi);
  h->SetDirectory(nullptr);
  h->Sumw2();
  return h;
}

Profile NewProfile(const TString &tag) {
  Profile p;
  p.z_line = NewHist("z_line_" + tag, N_Z_BINS, 0, 10);
  p.z_side = NewHist("z_side_" + tag, N_Z_BINS, 0, 10);
  p.ring_line = NewHist("ring_line_" + tag, N_RINGS, -0.5, N_RINGS - 0.5);
  p.ring_side = NewHist("ring_side_" + tag, N_RINGS, -0.5, N_RINGS - 0.5);
  return p;
}

// One run's single-interaction events into the spectrum and each window's
// profile.
void FillRun(const TString &run, TH1D *spectrum,
             const std::vector<LineWindow> &windows,
             std::vector<Profile> &profiles) {
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  if (!file) {
    std::cout << "CAPTURE96 missing " << run << std::endl;
    return;
  }
  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  Float_t energy = 0, x = 0, y = 0, z = 0;
  UInt_t event_time = 0;
  Int_t live_time = 0, n_interactions = 0;
  tree->SetBranchStatus("*", 0);
  const char *branches[7] = {"energykeV", "xmm",      "ymm",          "zmm",
                             "eventTime", "liveTime", "nInteractions"};
  for (Int_t b = 0; b < 7; b++)
    tree->SetBranchStatus(branches[b], 1);
  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("zmm", &z);
  tree->SetBranchAddress("eventTime", &event_time);
  tree->SetBranchAddress("liveTime", &live_time);
  tree->SetBranchAddress("nInteractions", &n_interactions);
  const Long64_t n_entries = tree->GetEntries();
  Long64_t previous_time = -1;
  for (Long64_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    Long64_t delta_time =
        previous_time < 0 ? 1 : (Long64_t)event_time - previous_time;
    previous_time = event_time;
    if (n_interactions != 1 || x == 0 || y == 0)
      continue;
    Float_t live_time_us = live_time * Constants::TENS_OF_NS_TO_S * 1e6;
    if (live_time_us < Constants::PILEUP_LIVETIME_THRESHOLD_US &&
        delta_time == 0)
      continue;
    Int_t ix = PixelIndex(x, Constants::PIXEL_CENTERS_X_MM);
    Int_t iy = PixelIndex(y, Constants::PIXEL_CENTERS_Y_MM);
    if (ix < 0 || iy < 0 || g_pixel_gain[ix][iy] <= 0)
      continue;
    Double_t e = energy * g_pixel_gain[ix][iy];
    if (spectrum)
      spectrum->Fill(e);
    Int_t ring = Ring(ix, iy);
    for (size_t w = 0; w < windows.size(); w++) {
      const LineWindow &win = windows[w];
      if (InWindow(e, win.lo, win.hi)) {
        profiles[w].z_line->Fill(z);
        profiles[w].ring_line->Fill(ring);
      } else if (InWindow(e, win.s1_lo, win.s1_hi) ||
                 InWindow(e, win.s2_lo, win.s2_hi)) {
        profiles[w].z_side->Fill(z);
        profiles[w].ring_side->Fill(ring);
      }
    }
  }
  std::cout << "CAPTURE96 " << run << ": " << n_entries << " rows" << std::endl;
  file->Close();
}

// The window's depth profile less the side bands', normalised to one.
TH1D *NetDepth(const Profile &p, const LineWindow &win, const TString &name,
               Double_t &net) {
  TH1D *h = (TH1D *)p.z_line->Clone(name);
  h->SetDirectory(nullptr);
  h->Add(p.z_side, -win.SideScale());
  net = h->Integral();
  h->Scale(1.0 / net);
  return h;
}

// Per-pixel rate of the line over that of the side bands, per ring, with the
// central rings (3-5) at one.
TH1D *RingRatio(const Profile &p, const LineWindow &win, const TString &name) {
  TH1D *h = (TH1D *)p.ring_line->Clone(name);
  h->SetDirectory(nullptr);
  h->Add(p.ring_side, -win.SideScale());
  h->Divide(p.ring_side);
  Double_t sw = 0, swx = 0;
  for (Int_t b = 4; b <= N_RINGS; b++) {
    Double_t w = 1.0 / TMath::Power(h->GetBinError(b), 2);
    sw += w;
    swx += w * h->GetBinContent(b);
  }
  h->Scale(sw / swx);
  return h;
}

Double_t NearCathode(TH1D *h, Double_t depth_mm) {
  return h->Integral(h->FindBin(10.0 - depth_mm + 0.001), N_Z_BINS);
}

const Int_t COLORS[6] = {kBlack,    kRed + 1,   kOrange + 7,
                         kBlue + 1, kAzure + 7, kGreen + 2};

void DrawDepths(const std::vector<TH1D *> &hists,
                const std::vector<TString> &labels,
                const std::vector<Int_t> &colors,
                const std::vector<Int_t> &styles, const TString &name) {
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels), 800,
      kTRUE);
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  for (size_t k = 0; k < hists.size(); k++) {
    PlottingUtils::ConfigureHistogram(
        hists[k], colors[k], ";Depth z [mm] (10 = cathode);Fraction / 0.25 mm");
    hists[k]->SetStats(0);
    hists[k]->SetLineStyle(styles[k]);
    hists[k]->SetMinimum(1e-3);
    hists[k]->SetMaximum(1.0);
    hists[k]->GetXaxis()->SetRangeUser(4.0, 10.0);
    hists[k]->Draw(k == 0 ? "hist e" : "hist e same");
    leg->AddEntry(hists[k], labels[k], "l");
  }
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLOG);
}

// The spectrum around the windows on a linear axis, each line's window and
// side bands shaded in its colour.
void PlotWindows(TH1D *spectrum, const std::vector<LineWindow> &windows,
                 const TString &label, const TString &name) {
  const Int_t window_color[2] = {kRed + 1, kBlue + 1};
  const Int_t side_color[2] = {kRed - 9, kAzure - 9};
  std::vector<TString> labels;
  for (size_t w = 0; w < windows.size(); w++) {
    labels.push_back(Form("%s window", windows[w].name));
    labels.push_back(Form("%s side bands", windows[w].name));
  }
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels, label));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size() + 1);
  leg->SetHeader(label);
  PlottingUtils::ConfigureHistogram(spectrum, kBlack, "");
  spectrum->SetStats(0);
  spectrum->GetXaxis()->SetRangeUser(78, 105);
  Double_t y_lo = 0.9 * spectrum->GetMinimum();
  Double_t y_hi = 1.1 * spectrum->GetMaximum();
  spectrum->SetMinimum(y_lo);
  spectrum->SetMaximum(y_hi);
  spectrum->Draw("hist");
  for (size_t w = 0; w < windows.size() && w < 2; w++) {
    const Double_t lo[3] = {windows[w].lo, windows[w].s1_lo, windows[w].s2_lo};
    const Double_t hi[3] = {windows[w].hi, windows[w].s1_hi, windows[w].s2_hi};
    for (Int_t k = 0; k < 3; k++) {
      TBox *box = new TBox(lo[k], y_lo, hi[k], y_hi);
      box->SetFillColorAlpha(k == 0 ? window_color[w] : side_color[w], 0.4);
      box->SetLineWidth(0);
      box->Draw();
      if (k < 2)
        leg->AddEntry(box, labels[2 * w + k], "f");
    }
  }
  spectrum->Draw("hist same");
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
}

std::vector<Configuration> Configurations() {
  std::vector<Configuration> configs;
  Configuration cu = {"cu_0113", "Cu shield, 01/13", {}};
  cu.runs.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
  cu.runs.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);
  Configuration cd10 = {"cd10_0113", "Cd shield 10%, 01/13", {}};
  cd10.runs.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  cd10.runs.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  Configuration cd25 = {"cd25_0113", "Cd shield 25%, 01/13", {}};
  cd25.runs.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
  cd25.runs.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
  Configuration none = {"none_0115", "No shield, 01/15", {}};
  none.runs.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_20260115);
  none.runs.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115);
  Configuration castle = {"castle_0116", "Graphite castle, 01/16", {}};
  castle.runs.push_back(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116);
  castle.runs.push_back(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116);
  configs.push_back(cu);
  configs.push_back(cd10);
  configs.push_back(cd25);
  configs.push_back(none);
  configs.push_back(castle);
  return configs;
}

void SetUp(const TString &project_root) {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root +
                                    "/plots/czt/studies/capture96_positions",
                                project_root + "/root_files");
  PlottingUtils::SetResidualHistPlotsEnabled(kFALSE);
  TGaxis::SetMaxDigits(3);
}

} // namespace Capture96

void Capture96Positions() {
  using namespace Capture96;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  SetUp(project_root);
  if (!LoadPixelGain()) {
    std::cout << "CAPTURE96 no pixel gain table" << std::endl;
    return;
  }
  std::vector<Configuration> configs = Configurations();

  std::vector<LineWindow> windows;
  windows.push_back(LINE96);
  windows.push_back(PB_KBETA);
  std::vector<TH1D *> depth96, rings;
  std::vector<TString> labels;
  for (size_t c = 0; c < configs.size(); c++) {
    TString tag = configs[c].tag;
    TH1D *spectrum =
        NewHist("spectrum_" + tag, N_GRID_BINS, GRID_START_KEV,
                GRID_START_KEV + N_GRID_BINS * Constants::BIN_WIDTH_KEV);
    spectrum->SetTitle(
        Form(";Energy [keV];Counts / %d eV", Constants::BIN_WIDTH_EV));
    std::vector<Profile> profiles;
    for (size_t w = 0; w < windows.size(); w++)
      profiles.push_back(NewProfile(Form("%s_%d", tag.Data(), (Int_t)w)));
    for (size_t r = 0; r < configs[c].runs.size(); r++)
      FillRun(configs[c].runs[r], spectrum, windows, profiles);

    Double_t net96 = 0, net_pb = 0;
    TH1D *z96 = NetDepth(profiles[0], LINE96, "z96_" + tag, net96);
    TH1D *z_pb = NetDepth(profiles[1], PB_KBETA, "zpb_" + tag, net_pb);
    TH1D *ring96 = RingRatio(profiles[0], LINE96, "ring96_" + tag);
    std::cout << "CAPTURE96 " << configs[c].label << ": net 95.9 keV "
              << std::fixed << std::setprecision(0) << net96 << ", Pb Kb "
              << net_pb << " | within 1 mm of the cathode: 95.9 keV "
              << std::setprecision(3) << NearCathode(z96, 1.0) << ", Pb Kb "
              << NearCathode(z_pb, 1.0)
              << " | last 0.25 mm: " << NearCathode(z96, 0.25) << ", "
              << NearCathode(z_pb, 0.25) << std::endl;
    std::cout << "CAPTURE96 " << configs[c].label
              << ": 95.9 keV ring rate / central";
    for (Int_t b = 1; b <= N_RINGS; b++)
      std::cout << "  " << b - 1 << ": " << std::setprecision(2)
                << ring96->GetBinContent(b) << " +- " << ring96->GetBinError(b);
    std::cout << std::endl;
    depth96.push_back(z96);
    rings.push_back(ring96);
    labels.push_back(configs[c].label);

    std::vector<TH1D *> pair;
    pair.push_back(z96);
    pair.push_back(z_pb);
    std::vector<TString> pair_labels;
    pair_labels.push_back(TString(configs[c].label) + ", 95.9 keV");
    pair_labels.push_back(TString(configs[c].label) + ", Pb K#beta");
    std::vector<Int_t> pair_colors, pair_styles;
    pair_colors.push_back(kRed + 1);
    pair_colors.push_back(kBlue + 1);
    pair_styles.push_back(1);
    pair_styles.push_back(2);
    DrawDepths(pair, pair_labels, pair_colors, pair_styles,
               "capture96_depth_vs_pb_" + tag);

    PlotWindows(spectrum, windows, configs[c].label, "spectrum_" + tag);
  }

  // The Am-241 source line: photons entering the cathode.
  std::vector<LineWindow> am_window;
  am_window.push_back(AM59);
  std::vector<Profile> am;
  am.push_back(NewProfile("am"));
  FillRun(Constants::POSTREACTOR_AM241_20260113, nullptr, am_window, am);
  Double_t net_am = 0;
  TH1D *z_am = NetDepth(am[0], AM59, "z_am", net_am);
  std::cout << "CAPTURE96 Am-241 01/13: within 1 mm of the cathode "
            << std::setprecision(3) << NearCathode(z_am, 1.0)
            << ", last 0.25 mm " << NearCathode(z_am, 0.25) << std::endl;

  std::vector<TH1D *> all(depth96);
  all.push_back(z_am);
  std::vector<TString> all_labels(labels);
  all_labels.push_back("Am-241 59.54 keV source, 01/13");
  std::vector<Int_t> colors, styles;
  for (size_t c = 0; c < depth96.size(); c++) {
    colors.push_back(COLORS[c + 1]);
    styles.push_back(1);
  }
  colors.push_back(COLORS[0]);
  styles.push_back(2);
  DrawDepths(all, all_labels, colors, styles, "capture96_depth");

  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels));
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  for (size_t c = 0; c < rings.size(); c++) {
    PlottingUtils::ConfigureHistogram(
        rings[c], COLORS[c + 1],
        ";Pixel ring from the crystal edge;Signal / background, centre = 1");
    rings[c]->SetStats(0);
    rings[c]->SetMarkerStyle(20 + c);
    rings[c]->SetMinimum(0.5);
    rings[c]->SetMaximum(1.5);
    rings[c]->Draw(c == 0 ? "e1 p" : "e1 p same");
    leg->AddEntry(rings[c], labels[c], "pe");
  }
  PlottingUtils::SaveFigure(canvas, "capture96_rings", "",
                            PlotSaveOptions::kLINEAR);
}
