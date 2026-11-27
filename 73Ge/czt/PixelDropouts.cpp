// Pixels that report no interactions in a run although they work in others:
// in every CZT run, each interaction is counted in the pixel whose centre is
// nearest in x and in y, from the full event tree (filtered/<run>.root's
// bef_tree, read only). A pixel is silent if it holds less than
// SILENT_FRACTION of the run's median pixel; a crystal's edge pixels, which
// the M400 gives positions only on their inner side, still hold about half
// of it. A run whose median pixel holds fewer than MIN_MEDIAN_COUNTS
// interactions is not judged.
//
// Writes results/pixel_dropouts.result and the dissertation appendix's table
// of the runs with silent pixels.
//
// Usage: root -b -q 'PixelDropouts.cpp+'
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "ProjectPaths.hpp"

#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

const Double_t SILENT_FRACTION = 0.01;
const Double_t MIN_MEDIAN_COUNTS = 100.0;

const TString TABLE_PATH =
    TString(gSystem->HomeDirectory()) +
    "/Dissertation/appendices/appendix_czt_dropouts/czt_silent_pixels_table.tex";

// Each run as the dissertation names it; day 1 is 01/13 and day 2 01/15.
struct RunName {
  TString run;
  TString label;
};

std::vector<RunName> Runs() {
  std::vector<RunName> runs;
  runs.push_back(
      {Constants::PASSIVEBACKGROUND_20260112, "01/12 passive background"});
  runs.push_back({Constants::CALIBRATION_20260112, "01/12 calibration"});
  runs.push_back({Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113,
                  "01/13 active background test, 5\\%"});
  runs.push_back({Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113,
                  "01/13 active background test, 90\\%"});
  runs.push_back({Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                  "01/13 Cu-shield signal, 10\\%"});
  runs.push_back({Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                  "01/13 Cu-shield background, 10\\%"});
  runs.push_back({Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                  "01/13 Cd-shield signal, 10\\%"});
  runs.push_back({Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                  "01/13 Cd-shield background, 10\\%"});
  runs.push_back({Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                  "01/13 Cd-shield signal, 25\\%"});
  runs.push_back({Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                  "01/13 Cd-shield background, 25\\%"});
  runs.push_back(
      {Constants::POSTREACTOR_AM241_20260113, "01/13 \\iso{241}{Am} source"});
  runs.push_back({Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                  "01/14 Cu-shield signal, 10\\%"});
  runs.push_back({Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                  "01/14 Cu-shield background, 10\\%"});
  runs.push_back({Constants::CUSHIELDSIGNAL_90PERCENT_20260114,
                  "01/14 Cu-shield signal, 90\\%"});
  runs.push_back({Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                  "01/15 no-shield signal, 5\\%"});
  runs.push_back({Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                  "01/15 no-shield background, 5\\%"});
  runs.push_back(
      {Constants::POSTREACTOR_AM241_20260115, "01/15 \\iso{241}{Am} source"});
  runs.push_back(
      {Constants::POSTREACTOR_BA133_20260115, "01/15 \\iso{133}{Ba} source"});
  runs.push_back({Constants::SHUTTERCLOSED_20260115, "01/15 shutter closed"});
  runs.push_back({Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
                  "01/16 Ge on the CZT, 0.5\\%"});
  runs.push_back({Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
                  "01/16 active background, 0.5\\%"});
  runs.push_back({Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
                  "01/16 graphite castle signal, 10\\%"});
  runs.push_back(
      {Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
       "01/16 graphite castle background, 10\\%"});
  runs.push_back({Constants::POSTREACTOR_AM241_BA133_20260116,
                  "01/16 \\iso{241}{Am} and \\iso{133}{Ba} source"});
  return runs;
}

// The index of the pixel centre nearest pos: the pixel whose area holds it.
Int_t NearestCentre(Float_t pos, const std::vector<Float_t> &centres) {
  Int_t best = 0;
  for (Int_t i = 1; i < (Int_t)centres.size(); i++)
    if (TMath::Abs(pos - centres[i]) < TMath::Abs(pos - centres[best]))
      best = i;
  return best;
}

struct SilentPixel {
  Int_t ix = 0, iy = 0;
  Long64_t counts = 0;
};

struct RunPixels {
  Bool_t read = kFALSE;
  Long64_t interactions = 0;
  Double_t median = 0;
  Bool_t judged = kFALSE;
  std::vector<SilentPixel> silent;
};

RunPixels CountPixels(const TString &run) {
  RunPixels out;
  TFile *file = IO::OpenForReading("filtered/" + run + ".root");
  TTree *tree = (file && !file->IsZombie())
                    ? static_cast<TTree *>(file->Get("bef_tree"))
                    : nullptr;
  if (!tree) {
    std::cerr << "ERROR: no bef_tree in filtered/" << run << ".root"
              << std::endl;
    delete file;
    return out;
  }
  const std::vector<Float_t> &cx = Constants::PIXEL_CENTERS_X_MM;
  const std::vector<Float_t> &cy = Constants::PIXEL_CENTERS_Y_MM;
  std::vector<std::vector<Long64_t>> counts(
      cx.size(), std::vector<Long64_t>(cy.size(), 0));
  Float_t x = 0, y = 0;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("xmm", 1);
  tree->SetBranchStatus("ymm", 1);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  out.interactions = tree->GetEntries();
  for (Long64_t i = 0; i < out.interactions; i++) {
    tree->GetEntry(i);
    counts[NearestCentre(x, cx)][NearestCentre(y, cy)]++;
  }
  file->Close();
  delete file;

  std::vector<Long64_t> all;
  for (size_t ix = 0; ix < cx.size(); ix++)
    for (size_t iy = 0; iy < cy.size(); iy++)
      all.push_back(counts[ix][iy]);
  std::sort(all.begin(), all.end());
  out.median = 0.5 * (all[(all.size() - 1) / 2] + all[all.size() / 2]);
  out.read = kTRUE;
  out.judged = out.median >= MIN_MEDIAN_COUNTS;
  if (!out.judged)
    return out;
  for (size_t ix = 0; ix < cx.size(); ix++)
    for (size_t iy = 0; iy < cy.size(); iy++)
      if (counts[ix][iy] < SILENT_FRACTION * out.median) {
        SilentPixel p;
        p.ix = (Int_t)ix;
        p.iy = (Int_t)iy;
        p.counts = counts[ix][iy];
        out.silent.push_back(p);
      }
  return out;
}

// A pixel as its centre, (x, y) in mm.
TString CentreText(const SilentPixel &p) {
  return TString::Format("(%.1f, %.1f)", Constants::PIXEL_CENTERS_X_MM[p.ix],
                         Constants::PIXEL_CENTERS_Y_MM[p.iy]);
}

} // namespace

void PixelDropouts() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/czt/studies",
                                project_root + "/root_files");
  const std::vector<RunName> runs = Runs();
  std::vector<RunPixels> results;
  for (size_t r = 0; r < runs.size(); r++) {
    results.push_back(CountPixels(runs[r].run));
    const RunPixels &p = results.back();
    std::cout << "PIXELS " << runs[r].run << ": " << p.interactions
              << " interactions, median pixel " << p.median;
    if (!p.judged)
      std::cout << ", too few to judge";
    else
      std::cout << ", " << p.silent.size() << " silent";
    for (size_t k = 0; k < p.silent.size(); k++)
      std::cout << " " << CentreText(p.silent[k]) << ":" << p.silent[k].counts;
    std::cout << std::endl;
  }

  const TString result_path = project_root + "/results/pixel_dropouts.result";
  std::ofstream result(result_path.Data());
  result << "# PixelDropouts.cpp: per run, the interactions, the median "
            "pixel's count, and each silent pixel (below "
         << SILENT_FRACTION
         << " of the median) as its centre in mm and its count; runs whose "
            "median pixel holds fewer than "
         << MIN_MEDIAN_COUNTS << " are not judged" << std::endl;
  for (size_t r = 0; r < runs.size(); r++) {
    const RunPixels &p = results[r];
    result << runs[r].run << " interactions " << p.interactions << " median "
           << p.median << " silent "
           << (p.judged ? TString::Format("%d", (Int_t)p.silent.size())
                        : TString("unjudged"));
    for (size_t k = 0; k < p.silent.size(); k++)
      result << " " << Constants::PIXEL_CENTERS_X_MM[p.silent[k].ix] << ","
             << Constants::PIXEL_CENTERS_Y_MM[p.silent[k].iy] << ":"
             << p.silent[k].counts;
    result << std::endl;
  }
  result.close();
  std::cout << "Wrote " << result_path << std::endl;

  // The appendix's table: the runs with silent pixels, each pixel by its
  // centre and its count against the run's median pixel.
  gSystem->mkdir(gSystem->GetDirName(TABLE_PATH), kTRUE);
  std::ofstream table(TABLE_PATH.Data());
  table << "% Written by UM-ANSG 73Ge/czt/PixelDropouts.cpp; edit that, not "
           "this file."
        << std::endl;
  table << "\\small" << std::endl;
  table << "\\begin{tabular}{@{}lrrl@{}}" << std::endl;
  table << "\\toprule" << std::endl;
  table << "Run & Interactions & Median pixel & Silent pixels, $(x, y)$ "
           "[mm]: interactions \\\\"
        << std::endl;
  table << "\\midrule" << std::endl;
  for (size_t r = 0; r < runs.size(); r++) {
    const RunPixels &p = results[r];
    if (!p.judged || p.silent.empty())
      continue;
    for (size_t k = 0; k < p.silent.size(); k++) {
      if (k == 0)
        table << runs[r].label << " & "
              << TString::Format("%.1f M", 1e-6 * p.interactions) << " & "
              << TString::Format("%.0f", p.median);
      else
        table << " & & ";
      table << " & " << CentreText(p.silent[k]) << ": " << p.silent[k].counts
            << " \\\\" << std::endl;
    }
  }
  table << "\\bottomrule" << std::endl;
  table << "\\end{tabular}" << std::endl;
  table.close();
  std::cout << "Wrote " << TABLE_PATH << std::endl;
}
