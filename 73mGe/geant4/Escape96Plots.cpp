// The 95.902 keV line's K X-ray escapes in the Geant4 runs of
// ~/Geant4/ANSG/Cd-Detector-Escape-Peaks (its run_escape96.sh writes them to
// root_files/geant4/<geometry>_<run>): the simulated spectra with every
// escape line marked at the capture line's energy minus the X-ray's (DDEP),
// as simulated (sharp lines) and smeared with the CdTe source runs' width;
// each sub-line's intensity and position; and the four families'
// intensities across the runs.
//
// The dissertation's figures (plots/czt/figures) are the CZT ones: the
// in-crystal run from 60 to 80 keV, smeared with the escape lines' measured
// width (results/escape_fractions.result's capture_sigma, so
// czt/CaptureEscapeFraction.cpp runs first), and the families of captures
// in the crystal against 95.902 keV photons from outside it.
//
// A sub-line's counts are those within 40 eV of it (the sub-lines of the
// K-beta1' groups, 34-51 eV apart, share one window) less the continuum of
// side bands 60-120 eV either side, over the full-energy line's counted the
// same way. A conversion escape leaves exactly 95.902 keV minus the DDEP
// X-ray energy the generator emits; a photon escape leaves 95.902 keV minus
// Geant4's own X-ray energy, which can differ by ~20 eV, so the window holds
// both.
//
// Usage: root -b -q 'Escape96Plots.cpp+'
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "ReferenceLines.hpp"

#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Escape96 {

const Double_t E_LINE = ReferenceLines::CD113_CAPTURE;
const Double_t HALF_WINDOW = 0.040, SIDE_NEAR = 0.060, SIDE_FAR = 0.120;
// The CdTe source runs' Gaussian width near 70 keV (SourceCurvature.cpp).
const Double_t SMEAR_SIGMA = 0.33;
const Double_t SHOW_LO = 60.0, SHOW_HI = 97.5;
// The dissertation's CZT figure: the escapes alone, smeared with the width
// the CZT fits give the escape lines (CaptureEscapeFraction.cpp's
// capture_sigma, the 01/15 capture line's).
const Double_t CZT_SHOW_LO = 60.0, CZT_SHOW_HI = 80.0;
const Double_t CZT_BIN_KEV = 0.05;

struct Run {
  const char *tag, *label;
};

const Int_t N_RUNS = 5;
// The CZT runs the dissertation shows: captures in the crystal, and photons
// from captures outside it (a Cd shield) onto its face.
const char *const IN_CRYSTAL_TAG = "czt96_surface";
const char *const OUTSIDE_TAG = "czt96_beam";
const Run RUNS[N_RUNS] = {{"cdte96_surface", "CdTe, captures near surface"},
                          {"cdte96_surface_noconv", "CdTe, no conversion"},
                          {"cdte96_beam", "CdTe, external beam"},
                          {"czt96_surface", "CZT, captures near surface"},
                          {"czt96_beam", "CZT, external beam"}};

// The escapes' sub-lines, DDEP X-ray energies (keV); a K-beta1' group's two
// components bound its one window.
struct SubLine {
  const char *name;
  Int_t family; // 0 Cd K-alpha, 1 Cd K-beta, 2 Te K-alpha, 3 Te K-beta
  Double_t x_lo, x_hi;
};

const Int_t N_SUB = 8;
const SubLine SUB[N_SUB] = {
    {"Cd Ka1", 0, 23.1739, 23.1739},     {"Cd Ka2", 0, 22.9843, 22.9843},
    {"Cd Kb1+Kb3", 1, 26.0615, 26.0958}, {"Cd Kb2", 1, 26.6440, 26.6440},
    {"Te Ka1", 2, 27.4726, 27.4726},     {"Te Ka2", 2, 27.2020, 27.2020},
    {"Te Kb1+Kb3", 3, 30.9446, 30.9960}, {"Te Kb2", 3, 31.7008, 31.7008}};

const Int_t N_FAMILIES = 4;
const char *const FAMILY_NAME[N_FAMILIES] = {"Cd K#alpha", "Cd K#beta",
                                             "Te K#alpha", "Te K#beta"};
const Int_t FAMILY_COLOR[N_FAMILIES] = {kRed + 1, kOrange + 7, kBlue + 1,
                                        kAzure + 7};
const char *const GE_LABEL = "^{73}Ge 68.752 keV #gamma ray";

struct Counted {
  Double_t counts = 0, err = 0, mean = 0, mode = 0;
};

// Counts in [lo, hi] (keV, deposited) less the side bands' continuum, the
// continuum-subtracted mean position and the fullest bin.
Counted CountWindow(TH1D *h, Double_t lo, Double_t hi) {
  Counted c;
  Int_t b_lo = h->FindBin(lo), b_hi = h->FindBin(hi);
  Int_t s1 = h->FindBin(lo - SIDE_FAR), s2 = h->FindBin(lo - SIDE_NEAR);
  Int_t s3 = h->FindBin(hi + SIDE_NEAR), s4 = h->FindBin(hi + SIDE_FAR);
  Double_t side = h->Integral(s1, s2) + h->Integral(s3, s4);
  Double_t side_per_bin = side / (s2 - s1 + 1 + s4 - s3 + 1);
  Double_t peak = h->Integral(b_lo, b_hi);
  Int_t n_bins = b_hi - b_lo + 1;
  c.counts = peak - side_per_bin * n_bins;
  c.err = TMath::Sqrt(peak + side_per_bin * n_bins * n_bins /
                                 (s2 - s1 + 1 + s4 - s3 + 1));
  Double_t sw = 0, swx = 0, best = -1;
  for (Int_t b = b_lo; b <= b_hi; b++) {
    Double_t w = h->GetBinContent(b) - side_per_bin;
    if (w > 0) {
      sw += w;
      swx += w * h->GetBinCenter(b);
    }
    if (h->GetBinContent(b) > best) {
      best = h->GetBinContent(b);
      c.mode = h->GetBinCenter(b);
    }
  }
  c.mean = sw > 0 ? swx / sw : 0;
  return c;
}

// A "key value" line's value from a results file; 0 if it is not there.
Double_t ReadResultValue(const TString &path, const TString &key) {
  std::ifstream in(path.Data());
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream fields(line);
    std::string name;
    Double_t value = 0;
    if (fields >> name >> value && TString(name) == key)
      return value;
  }
  return 0;
}

// Every crystal's deposits of one run, in 1 eV bins.
TH1D *LoadRun(const TString &dir, const TString &name) {
  TH1D *h = new TH1D(name, "", 50000, 50.0, 100.0);
  h->SetDirectory(nullptr);
  for (Int_t crystal = 0; crystal < 4; crystal++) {
    TChain chain(Form("CZT_%d", crystal));
    // The sim merges its per-thread files into output0.root; older runs
    // kept one file per thread.
    if (gSystem->AccessPathName(dir + "/output0.root"))
      chain.Add(dir + "/output0_t*.root");
    else
      chain.Add(dir + "/output0.root");
    if (chain.GetEntries() == 0)
      continue;
    Double_t edep = 0;
    chain.SetBranchAddress("fEDep", &edep);
    const Long64_t n_entries = chain.GetEntries();
    for (Long64_t i = 0; i < n_entries; i++) {
      chain.GetEntry(i);
      h->Fill(edep);
    }
  }
  return h;
}

// h folded with a Gaussian of width sigma (keV), into bins of width_kev
// over [lo, hi].
TH1D *Smear(TH1D *h, Double_t sigma, Double_t width_kev, const TString &name,
            Double_t lo = SHOW_LO, Double_t hi = SHOW_HI) {
  Int_t n_out = (Int_t)((hi - lo) / width_kev + 0.5);
  TH1D *out = new TH1D(name, "", n_out, lo, hi);
  out->SetDirectory(nullptr);
  for (Int_t b = 1; b <= h->GetNbinsX(); b++) {
    Double_t n = h->GetBinContent(b);
    if (n <= 0)
      continue;
    Double_t x = h->GetBinCenter(b);
    Int_t first = out->FindBin(x - 5.0 * sigma);
    Int_t last = out->FindBin(x + 5.0 * sigma);
    for (Int_t o = TMath::Max(first, 1); o <= TMath::Min(last, n_out); o++) {
      Double_t lo = (out->GetXaxis()->GetBinLowEdge(o) - x) / sigma;
      Double_t hi = (out->GetXaxis()->GetBinUpEdge(o) - x) / sigma;
      out->AddBinContent(o, n * (TMath::Freq(hi) - TMath::Freq(lo)));
    }
  }
  return out;
}

void MarkLines(Double_t y_lo, Double_t y_hi, TLegend *leg) {
  for (Int_t f = 0; f < N_FAMILIES; f++) {
    for (Int_t s = 0; s < N_SUB; s++) {
      if (SUB[s].family != f)
        continue;
      Double_t e = E_LINE - 0.5 * (SUB[s].x_lo + SUB[s].x_hi);
      TLine *line = new TLine(e, y_lo, e, y_hi);
      line->SetLineColor(FAMILY_COLOR[f]);
      line->SetLineStyle(2);
      line->Draw();
      if (s == 0 || SUB[s - 1].family != f)
        leg->AddEntry(line, Form("%s escape", FAMILY_NAME[f]), "l");
    }
  }
  TLine *ge = new TLine(ReferenceLines::GE73M_ENSDF, y_lo,
                        ReferenceLines::GE73M_ENSDF, y_hi);
  ge->SetLineColor(kGreen + 2);
  // Solid, to stand out from the dashed escapes: the Te K-alpha2 escape lies
  // 52 eV below it.
  ge->SetLineWidth(2 * PlottingUtils::GetLineWidth());
  ge->SetLineStyle(1);
  ge->Draw();
  leg->AddEntry(ge, GE_LABEL, "l");
}

const Int_t PLOT_HEIGHT_PX = 800;

void PlotSpectrum(TH1D *h, const TString &title, const TString &label,
                  const TString &note, const TString &name,
                  Double_t show_lo = SHOW_LO, Double_t show_hi = SHOW_HI) {
  std::vector<TString> labels;
  labels.push_back(note);
  for (Int_t f = 0; f < N_FAMILIES; f++)
    labels.push_back(Form("%s escape", FAMILY_NAME[f]));
  labels.push_back(GE_LABEL);
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels, label),
      PLOT_HEIGHT_PX, kTRUE);
  plot->cd();
  h->SetTitle(title);
  h->SetStats(0);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetRangeUser(show_lo, show_hi);
  Double_t y_lo = 0.5 * TMath::Max(h->GetMinimum(0), 1.0);
  Double_t data_max = h->GetMaximum();
  h->SetMinimum(y_lo);
  const Double_t y_hi = 5.0 * data_max;
  h->SetMaximum(y_hi);
  h->Draw("hist");
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size() + 1);
  leg->SetHeader(label);
  leg->AddEntry((TObject *)nullptr, note, "");
  // The marks span the whole frame.
  MarkLines(y_lo, y_hi, leg);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLOG);
}

// The four families' intensities, one bar per run.
void PlotFamilies(const std::vector<TH1D *> &hists,
                  const std::vector<TString> &labels,
                  const std::vector<Int_t> &colors, const TString &title,
                  const TString &name) {
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels),
      PLOT_HEIGHT_PX, kTRUE);
  TLegend *leg =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  plot->cd();
  gStyle->SetPaintTextFormat(".2f");
  const Double_t width = 0.8 / hists.size();
  for (size_t k = 0; k < hists.size(); k++) {
    TH1D *fam = hists[k];
    for (Int_t f = 0; f < N_FAMILIES; f++)
      fam->GetXaxis()->SetBinLabel(f + 1, FAMILY_NAME[f]);
    fam->SetStats(0);
    fam->SetTitle(title);
    fam->SetFillColor(colors[k]);
    fam->SetLineColor(colors[k]);
    fam->SetBarWidth(width);
    fam->SetBarOffset(0.1 + width * k);
    fam->SetMinimum(0.03);
    fam->SetMaximum(20.0);
    // The bar values' text size (marker size times the style's), legible
    // when the figure is printed at a page's width.
    fam->SetMarkerSize(2.0);
    fam->Draw(k == 0 ? "bar text0" : "bar text0 same");
    leg->AddEntry(fam, labels[k], "f");
  }
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLOG);
}

} // namespace Escape96

void Escape96Plots() {
  using namespace Escape96;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG, project_root + "/plots/czt/studies/geant4_escapes",
      project_root + "/root_files");
  std::vector<TH1D *> family_hists;
  std::vector<Int_t> shown;
  TH1D *in_crystal_smeared = nullptr;
  const Double_t czt_sigma = ReadResultValue(
      project_root + "/results/escape_fractions.result", "capture_sigma");
  if (czt_sigma <= 0)
    std::cout << "ESCAPE96 no capture_sigma in results/escape_fractions.result "
                 "(run czt/CaptureEscapeFraction.cpp); no CZT spectrum figure"
              << std::endl;
  for (Int_t r = 0; r < N_RUNS; r++) {
    TString dir = project_root + "/root_files/geant4/" + RUNS[r].tag;
    TH1D *h = LoadRun(dir, Form("run_%d", r));
    if (h->GetEntries() == 0) {
      std::cout << "ESCAPE96 " << RUNS[r].tag << ": no output in " << dir
                << std::endl;
      continue;
    }
    Counted full = CountWindow(h, E_LINE - HALF_WINDOW, E_LINE + HALF_WINDOW);
    std::cout << "ESCAPE96 " << RUNS[r].tag << " (" << RUNS[r].label
              << "): full-energy " << std::fixed << std::setprecision(0)
              << full.counts << " counts" << std::endl;
    Double_t family[N_FAMILIES] = {0, 0, 0, 0};
    Double_t sub_counts[N_SUB];
    for (Int_t s = 0; s < N_SUB; s++) {
      Double_t lo = E_LINE - SUB[s].x_hi - HALF_WINDOW;
      Double_t hi = E_LINE - SUB[s].x_lo + HALF_WINDOW;
      Counted c = CountWindow(h, lo, hi);
      sub_counts[s] = c.counts;
      family[SUB[s].family] += c.counts;
      Double_t expected = E_LINE - 0.5 * (SUB[s].x_lo + SUB[s].x_hi);
      std::cout << "ESCAPE96 " << RUNS[r].tag << "   " << std::left
                << std::setw(11) << SUB[s].name << std::right << " expected "
                << std::setprecision(3) << expected << " keV, simulated mean "
                << c.mean << " (fullest bin " << c.mode << ")  "
                << std::setprecision(3) << 100.0 * c.counts / full.counts
                << " +- " << 100.0 * c.err / full.counts << " % of full energy"
                << std::endl;
    }
    TH1D *fam = new TH1D(Form("families_%d", r), "", N_FAMILIES, 0, N_FAMILIES);
    fam->SetDirectory(nullptr);
    std::cout << "ESCAPE96 " << RUNS[r].tag << "   families [%]:";
    for (Int_t f = 0; f < N_FAMILIES; f++) {
      fam->SetBinContent(f + 1, 100.0 * family[f] / full.counts);
      std::cout << " " << FAMILY_NAME[f] << " " << std::setprecision(3)
                << 100.0 * family[f] / full.counts;
    }
    // The ratios GainTransferFit.hpp takes from czt96_surface
    // (GEANT4_INSIDE_*).
    std::cout << "  | Cd Ka/Kb " << std::setprecision(4)
              << family[0] / family[1] << ", Te Ka/Cd Ka "
              << family[2] / family[0] << ", Te Kb/Te Ka "
              << family[3] / family[2] << std::endl;
    // Each line's share of its family, the escape shares
    // include/ReferenceLines.hpp takes from czt96_surface and cdte96_surface.
    std::cout << "ESCAPE96 " << RUNS[r].tag << "   shares in family:";
    for (Int_t s = 0; s < N_SUB; s++)
      std::cout << " " << SUB[s].name << " " << std::setprecision(4)
                << sub_counts[s] / family[SUB[s].family];
    std::cout << std::endl;
    family_hists.push_back(fam);
    shown.push_back(r);

    TH1D *sharp = (TH1D *)h->Rebin(10, Form("sharp_%d", r));
    PlotSpectrum(sharp, ";Deposited energy [keV];Counts / 10 eV", RUNS[r].label,
                 "as simulated", TString("escape96_") + RUNS[r].tag);
    if (TString(RUNS[r].tag) == IN_CRYSTAL_TAG && czt_sigma > 0)
      in_crystal_smeared =
          Smear(h, czt_sigma, CZT_BIN_KEV, "czt_in_crystal_smeared",
                CZT_SHOW_LO, CZT_SHOW_HI);
    TH1D *smeared = Smear(h, SMEAR_SIGMA, 0.05, Form("smeared_%d", r));
    PlotSpectrum(smeared, ";Deposited energy [keV];Counts / 50 eV",
                 RUNS[r].label,
                 Form("smeared, #sigma = %.0f eV", 1000.0 * SMEAR_SIGMA),
                 TString("escape96_") + RUNS[r].tag + "_smeared");
    delete h;
  }

  const Int_t run_color[N_RUNS] = {kRed + 1, kMagenta + 1, kOrange + 7,
                                   kBlue + 1, kAzure + 7};
  std::vector<TString> labels;
  std::vector<Int_t> colors;
  for (size_t k = 0; k < shown.size(); k++) {
    labels.push_back(RUNS[shown[k]].label);
    colors.push_back(run_color[shown[k]]);
  }
  PlotFamilies(family_hists, labels, colors,
               ";Escape family;Escape / full-energy counts [%]",
               "escape96_families");

  // The dissertation's figures (plots/czt/figures): the in-crystal run whose
  // ratios the CZT fits take (GainTransferFit.hpp's GEANT4_INSIDE_*), and
  // the CZT runs' families, captures in the crystal against 95.902 keV
  // photons from outside it (ESCAPES_OUTSIDE, the Cd shield's).
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/czt/figures",
                                project_root + "/root_files");
  if (in_crystal_smeared)
    PlotSpectrum(
        in_crystal_smeared,
        Form(";Deposited Energy [keV];Counts / %.0f eV", 1000.0 * CZT_BIN_KEV),
        "Geant4, captures in the crystal",
        Form("smeared, #sigma = %.0f eV", 1000.0 * czt_sigma),
        "czt_escape96_simulated", CZT_SHOW_LO, CZT_SHOW_HI);
  std::vector<TH1D *> czt_hists;
  std::vector<TString> czt_labels;
  std::vector<Int_t> czt_colors;
  for (size_t k = 0; k < shown.size(); k++) {
    const TString tag = RUNS[shown[k]].tag;
    if (tag == IN_CRYSTAL_TAG) {
      czt_hists.push_back(family_hists[k]);
      czt_labels.push_back("Captures in the crystal");
      czt_colors.push_back(kBlue + 1);
    } else if (tag == OUTSIDE_TAG) {
      czt_hists.push_back(family_hists[k]);
      czt_labels.push_back("95.902 keV photons from outside");
      czt_colors.push_back(kAzure + 7);
    }
  }
  if (!czt_hists.empty())
    PlotFamilies(czt_hists, czt_labels, czt_colors,
                 ";Escaping X Ray;Escape / Full-Energy Counts [%]",
                 "czt_escape96_families");
}
