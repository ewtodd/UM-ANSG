#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TTree.h>
#include <cstdlib>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

static std::mutex print_mutex;

static std::vector<std::vector<Float_t>> g_pixel_gain;

static void LoadPixelGainOrDie() {
  TFile *file = IO::OpenForReading("pixel_calibration.root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot open root_files/pixel_calibration.root. "
                 "Run PixelCalibration first."
              << std::endl;
    std::exit(1);
  }
  TH2D *gain = static_cast<TH2D *>(file->Get("pixel_gain"));
  if (!gain) {
    std::cerr << "ERROR: pixel_gain TH2D not found in "
                 "root_files/pixel_calibration.root."
              << std::endl;
    file->Close();
    delete file;
    std::exit(1);
  }
  Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();
  if (gain->GetNbinsX() != nx || gain->GetNbinsY() != ny) {
    std::cerr << "ERROR: pixel_gain dims " << gain->GetNbinsX() << "x"
              << gain->GetNbinsY() << " != expected " << nx << "x" << ny
              << std::endl;
    file->Close();
    delete file;
    std::exit(1);
  }
  g_pixel_gain.assign(nx, std::vector<Float_t>(ny, 0.0f));
  Int_t n_good = 0;
  for (Int_t ix = 0; ix < nx; ix++)
    for (Int_t iy = 0; iy < ny; iy++) {
      Float_t v = (Float_t)gain->GetBinContent(ix + 1, iy + 1);
      g_pixel_gain[ix][iy] = v;
      if (v > 0)
        n_good++;
    }
  file->Close();
  delete file;
  std::cout << "Loaded pixel gain table: " << n_good << " good / " << (nx * ny)
            << " pixels" << std::endl;
}

static inline Float_t PixelGain(Int_t ix, Int_t iy) {
  if (ix < 0 || iy < 0 || ix >= (Int_t)g_pixel_gain.size() ||
      iy >= (Int_t)g_pixel_gain[0].size())
    return 0.0f;
  return g_pixel_gain[ix][iy];
}

// Inverse-distance-weighted gain for an interaction at (x, y), using the
// GAIN_IDW_K_NEAREST nearest pixel centers (weight 1/d^2). Used by the
// total-energy tree, which keeps off-grid interactions instead of cutting
// them. Sets any_contributor_bad = kTRUE if any of the contributing centers
// is a bad/dead pixel (gain <= 0); the caller then drops the whole event.
static Float_t WeightedPixelGain(Float_t x, Float_t y,
                                 Bool_t &any_contributor_bad) {
  Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();
  Int_t k = Constants::GAIN_IDW_K_NEAREST;

  std::vector<Float_t> best_d2(k, 1e30f);
  std::vector<Int_t> best_ix(k, -1);
  std::vector<Int_t> best_iy(k, -1);

  for (Int_t ix = 0; ix < nx; ix++) {
    Float_t dx = x - Constants::PIXEL_CENTERS_X_MM[ix];
    for (Int_t iy = 0; iy < ny; iy++) {
      Float_t dy = y - Constants::PIXEL_CENTERS_Y_MM[iy];
      Float_t d2 = dx * dx + dy * dy;
      Int_t pos = k;
      while (pos > 0 && d2 < best_d2[pos - 1])
        pos--;
      if (pos < k) {
        for (Int_t m = k - 1; m > pos; m--) {
          best_d2[m] = best_d2[m - 1];
          best_ix[m] = best_ix[m - 1];
          best_iy[m] = best_iy[m - 1];
        }
        best_d2[pos] = d2;
        best_ix[pos] = ix;
        best_iy[pos] = iy;
      }
    }
  }

  any_contributor_bad = kFALSE;
  Float_t sum_w = 0;
  Float_t sum_wg = 0;
  for (Int_t m = 0; m < k; m++) {
    if (best_ix[m] < 0)
      continue;
    Float_t g = PixelGain(best_ix[m], best_iy[m]);
    if (g <= 0) {
      any_contributor_bad = kTRUE;
      return 0.0f;
    }
    Float_t d2 = best_d2[m];
    if (d2 < 1e-6f)
      d2 = 1e-6f;
    Float_t w = 1.0f / d2;
    sum_w += w;
    sum_wg += w * g;
  }
  if (sum_w <= 0)
    return 0.0f;
  return sum_wg / sum_w;
}

Bool_t IsOnPixelCenter(Float_t pos, const std::vector<Float_t> &centers) {
  for (Int_t i = 0; i < (Int_t)centers.size(); i++) {
    if (TMath::Abs(pos - centers[i]) <= Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      return kTRUE;
  }
  return kFALSE;
}

Int_t GetPixelIndex(Float_t pos, const std::vector<Float_t> &centers) {
  for (Int_t i = 0; i < (Int_t)centers.size(); i++) {
    if (TMath::Abs(pos - centers[i]) <= Constants::PIXEL_ACCEPT_HALFWIDTH_MM)
      return i;
  }
  return -1;
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
  return -1; // on boundary
}

std::vector<Float_t> GetTrimmedPixelCenters(const std::vector<Float_t> &centers,
                                            Int_t nEdgeSkip) {
  // centers is sorted. First half (negative values) is one crystal,
  // second half (positive values) is the other. Skip nEdgeSkip from
  // each end of each half.
  std::vector<Float_t> trimmed;
  Int_t halfSize = (Int_t)centers.size() / 2; // 11
  for (Int_t i = 0; i < halfSize; i++) {
    if (i >= nEdgeSkip && i < halfSize - nEdgeSkip)
      trimmed.push_back(centers[i]);
  }
  for (Int_t i = halfSize; i < (Int_t)centers.size(); i++) {
    Int_t j = i - halfSize;
    if (j >= nEdgeSkip && j < halfSize - nEdgeSkip)
      trimmed.push_back(centers[i]);
  }
  return trimmed;
}

void FilterDemo(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TFile *file =
        IO::OpenForWriting("filtered/" + filename + ".root", "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));

    Float_t energy = 0;
    Float_t x = 0, y = 0, z = 0;
    UInt_t eventTime = 0;
    Int_t liveTime = 0;
    Int_t nInteractions = 0;
    Int_t interaction = 0;

    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    tree->SetBranchAddress("zmm", &z);
    tree->SetBranchAddress("eventTime", &eventTime);
    tree->SetBranchAddress("liveTime", &liveTime);
    tree->SetBranchAddress("nInteractions", &nInteractions);
    tree->SetBranchAddress("interaction", &interaction);

    Int_t n_entries = tree->GetEntries();

    Int_t pos_n_bins = 1000;
    Int_t xy_range = 22;

    TH1F *X =
        new TH1F(PlottingUtils::GetRandomName(), "; X Position [mm]; Counts",
                 pos_n_bins, -xy_range, xy_range);

    Int_t z_n_bins = 500;
    Int_t z_min_mm = 0;

    TH2D *EZ = new TH2D(PlottingUtils::GetRandomName(),
                        "; Energy [keV]; Z Position [mm]",
                        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                        Constants::ZOOMED_XMAX, z_n_bins, z_min_mm, 10);

    TParameter<Float_t> *param =
        (TParameter<Float_t> *)file->Get("N42_RealTime_Total");

    if (!param) {
      std::cerr << "WARNING: Could not find N42_RealTime_Total parameter in "
                << filename << ".root" << std::endl;
    }

    TH1F *includedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    TH1F *excludedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    Bool_t in_excluded_region;

    Int_t event_time_prev = -1;
    Int_t delta_event_time = 0;

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      in_excluded_region = kFALSE;

      Int_t crystal = GetCrystalIndex(x, y);
      if (crystal < 0)
        in_excluded_region = kTRUE;

      if (nInteractions != 1)
        in_excluded_region = kTRUE;

      if (z < Constants::FILTER_DEPTH_MM)
        in_excluded_region = kTRUE;

      Int_t ix = GetPixelIndex(x, Constants::PIXEL_CENTERS_X_MM);
      Int_t iy = GetPixelIndex(y, Constants::PIXEL_CENTERS_Y_MM);
      if (ix < 0 || iy < 0)
        in_excluded_region = kTRUE;

      Float_t gain = PixelGain(ix, iy);
      if (gain <= 0)
        in_excluded_region = kTRUE;

      X->Fill(x);

      Float_t liveTime_us = liveTime * Constants::TENS_OF_NS_TO_S * 1e6;
      if (event_time_prev != -1)
        delta_event_time = eventTime - event_time_prev;
      if (liveTime_us < Constants::PILEUP_LIVETIME_THRESHOLD_US &&
          delta_event_time == 0)
        in_excluded_region = kTRUE;

      Float_t filledEnergy = (gain > 0) ? energy * gain : energy;

      if (filledEnergy > Constants::ZOOMED_XMIN &&
          filledEnergy < Constants::ZOOMED_XMAX) {
        if (z > z_min_mm)
          EZ->Fill(filledEnergy, z);
        if (in_excluded_region) {
          excludedSpectrum->Fill(filledEnergy);
        } else
          includedSpectrum->Fill(filledEnergy);
      }
      if (i != 0) {
        event_time_prev = eventTime;
      } else {
        event_time_prev = -1;
      };
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasX = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureAndDrawHistogram(X, kBlack);

    TSpectrum *spec = new TSpectrum(22);
    Int_t nFound = spec->Search(X, 1, "", 0.05);
    Double_t *centers = spec->GetPositionX();

    std::vector<Float_t> pixelCenters(centers, centers + nFound);
    std::sort(pixelCenters.begin(), pixelCenters.end());

    Float_t halfWidth = 0.02; // [mm]

    for (Int_t i = 0; i < (Int_t)pixelCenters.size(); i++) {
      TBox *box = new TBox(pixelCenters[i] - halfWidth, X->GetMinimum(),
                           pixelCenters[i] + halfWidth, X->GetMaximum());
      box->SetFillColorAlpha(kRed, 0.3);
      box->SetLineWidth(0);
      box->Draw();
    }

    PlottingUtils::SaveFigure(canvasX, filename + "_X", "filterDemo",
                              PlotSaveOptions::kLOG);

    TCanvas *canvasEZ = PlottingUtils::GetConfiguredCanvas();
    canvasEZ->SetRightMargin(0.2);
    EZ->Draw("COLZ");
    PlottingUtils::Configure2DHistogram(EZ, canvasEZ);
    canvasEZ->SetLogz(kFALSE);
    canvasEZ->Update();

    TBox *box = new TBox(Constants::ZOOMED_XMIN, z_min_mm,
                         Constants::ZOOMED_XMAX, Constants::FILTER_DEPTH_MM);
    box->SetLineColor(kBlack);
    box->SetLineWidth(3);
    box->SetLineStyle(2);
    box->SetFillColorAlpha(kBlack, 0.5);
    box->Draw();

    PlottingUtils::SaveFigure(canvasEZ, filename + "_ZvsE", "filterDemo",
                              PlotSaveOptions::kLINEAR);

    TCanvas *canvasRegions = PlottingUtils::GetConfiguredCanvas(kFALSE);

    Float_t maxY = TMath::Max(includedSpectrum->GetMaximum(),
                              excludedSpectrum->GetMaximum());

    PlottingUtils::ConfigureHistogram(includedSpectrum, kRed);
    includedSpectrum->GetYaxis()->SetRangeUser(0, maxY * 1.1);
    includedSpectrum->SetFillStyle(0);
    includedSpectrum->GetYaxis()->SetTitleOffset(1.5);
    includedSpectrum->Draw("HIST");

    PlottingUtils::ConfigureHistogram(excludedSpectrum, kBlack);
    excludedSpectrum->SetFillStyle(0);
    excludedSpectrum->Draw("HIST SAME");

    if (filename.Contains("Signal")) {
      TLine *line68 = new TLine(68.752, 0, 68.752, maxY * 1.1);
      line68->SetLineColor(kViolet);
      line68->SetLineWidth(2);
      line68->SetLineStyle(2);
      line68->Draw();

      TLatex *label68 =
          PlottingUtils::AddText("^{73m}Ge 68.75 keV #gamma ray", 0.47, 0.83);
      label68->SetTextColor(kViolet);
    }

    TLegend *leg = PlottingUtils::AddLegend(0.72, 0.9, 0.75, 0.88);
    leg->AddEntry(includedSpectrum, "Accepted", "l");
    leg->AddEntry(excludedSpectrum, "Rejected", "l");
    leg->Draw();

    canvasRegions->SetLeftMargin(0.2);
    PlottingUtils::SaveFigure(canvasRegions, filename + "_FilterSpectra",
                              "filterDemo", PlotSaveOptions::kLINEAR);

    std::cout << "Wrote histograms for " << filename << std::endl;

    EZ->Write("EZ", TObject::kOverwrite);
    X->Write("X", TObject::kOverwrite);
    canvasRegions->cd();
    canvasRegions->Write("AcceptedRejectedSpectra", TObject::kOverwrite);

    file->Close();
  }
}

// Build a single total-energy tree from bef_tree: one entry per event, summing
// the gain-corrected energy of every interaction (across quadrants). Drops the
// depth cut and the single-interaction requirement; keeps the pileup cut. An
// event is dropped if any interaction's IDW gain has a bad/dead contributor.
Bool_t BuildTotalEnergyTree(TFile *file) {
  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  if (!tree)
    return kFALSE;

  Float_t energy = 0;
  Float_t x = 0, y = 0, z = 0;
  UInt_t eventTime = 0;
  Int_t liveTime = 0;
  Int_t nInteractions = 0;
  Int_t interaction = 0;

  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("zmm", &z);
  tree->SetBranchAddress("eventTime", &eventTime);
  tree->SetBranchAddress("liveTime", &liveTime);
  tree->SetBranchAddress("nInteractions", &nInteractions);
  tree->SetBranchAddress("interaction", &interaction);

  Float_t outEnergy = 0;
  UInt_t outEventTime = 0;
  Int_t outLiveTime = 0;
  TTree *totalTree =
      new TTree("total_energy_tree", "Total deposited energy per event");
  totalTree->Branch("energykeV", &outEnergy, "energykeV/F");
  totalTree->Branch("eventTime", &outEventTime, "eventTime/i");
  totalTree->Branch("liveTime", &outLiveTime, "liveTime/I");

  Float_t totalLiveTime_s = 0;

  Float_t event_energy = 0;
  UInt_t event_eventTime = 0;
  Int_t event_liveTime = 0;
  Int_t event_n = 0;
  Int_t event_seen = 0;
  Bool_t event_dropped = kFALSE;

  Int_t event_time_prev = -1;
  Int_t delta_event_time = 0;

  Int_t n_entries = tree->GetEntries();
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    if (interaction == 0) {
      event_energy = 0;
      event_seen = 0;
      event_dropped = kFALSE;
      event_n = nInteractions;
      event_eventTime = eventTime;
      event_liveTime = liveTime;
    }

    Bool_t any_bad = kFALSE;
    Float_t gain = WeightedPixelGain(x, y, any_bad);
    if (any_bad)
      event_dropped = kTRUE;
    else
      event_energy += energy * gain;
    event_seen++;

    if (event_seen == event_n) {
      Float_t liveTime_us = event_liveTime * Constants::TENS_OF_NS_TO_S * 1e6;
      if (event_time_prev != -1)
        delta_event_time = (Int_t)event_eventTime - event_time_prev;
      Bool_t pileup = (liveTime_us < Constants::PILEUP_LIVETIME_THRESHOLD_US &&
                       delta_event_time == 0);

      if (!event_dropped && !pileup) {
        outEnergy = event_energy;
        outEventTime = event_eventTime;
        outLiveTime = event_liveTime;
        totalTree->Fill();
        totalLiveTime_s += event_liveTime * Constants::TENS_OF_NS_TO_S;
      }
      event_time_prev = (Int_t)event_eventTime;
    }
  }

  totalTree->Write("", TObject::kOverwrite);

  TParameter<Float_t> *ltParam =
      new TParameter<Float_t>("LiveTime_Total_s", totalLiveTime_s);
  ltParam->Write("", TObject::kOverwrite);

  {
    std::lock_guard<std::mutex> lock(print_mutex);
    std::cout << "  total_energy_tree: " << totalTree->GetEntries() << " events"
              << std::endl;
  }
  return kTRUE;
}

Bool_t FilterFile(TString filename) {
  TFile *file = IO::OpenForWriting("filtered/" + filename + ".root", "UPDATE");
  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));

  Float_t energy = 0;
  Float_t x = 0, y = 0, z = 0;
  UInt_t eventTime = 0;
  Int_t liveTime = 0;
  Int_t nInteractions = 0;
  Int_t interaction = 0;

  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("zmm", &z);
  tree->SetBranchAddress("eventTime", &eventTime);
  tree->SetBranchAddress("liveTime", &liveTime);
  tree->SetBranchAddress("nInteractions", &nInteractions);
  tree->SetBranchAddress("interaction", &interaction);

  Int_t n_entries = tree->GetEntries();

  Float_t outEnergy = 0;
  Float_t outX = 0, outY = 0, outZ = 0;
  UInt_t outEventTime = 0;
  Int_t outLiveTime = 0;

  Float_t filteredLiveTime_s[Constants::N_CRYSTALS];
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++)
    filteredLiveTime_s[c] = 0.0;

  TTree *filteredTrees[Constants::N_CRYSTALS];
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    filteredTrees[c] = new TTree(Form("crystal%d_filtered_tree", c),
                                 Form("Filtered events for crystal %d", c));
    filteredTrees[c]->Branch("energykeV", &outEnergy, "energykeV/F");
    filteredTrees[c]->Branch("xmm", &outX, "xmm/F");
    filteredTrees[c]->Branch("ymm", &outY, "ymm/F");
    filteredTrees[c]->Branch("zmm", &outZ, "zmm/F");
    filteredTrees[c]->Branch("eventTime", &outEventTime, "eventTime/i");
    filteredTrees[c]->Branch("liveTime", &outLiveTime, "liveTime/I");
  }

  Int_t event_time_prev = -1;
  Int_t delta_event_time = 0;

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    Int_t crystal = GetCrystalIndex(x, y);

    if (crystal < 0)
      continue;

    if (nInteractions != 1)
      continue;

    if (z < Constants::FILTER_DEPTH_MM)
      continue;

    Int_t ix = GetPixelIndex(x, Constants::PIXEL_CENTERS_X_MM);
    Int_t iy = GetPixelIndex(y, Constants::PIXEL_CENTERS_Y_MM);
    if (ix < 0 || iy < 0)
      continue;

    Float_t gain = PixelGain(ix, iy);
    if (gain <= 0)
      continue;

    Float_t gain_matched_energy = energy * gain;

    Float_t liveTime_us = liveTime * Constants::TENS_OF_NS_TO_S * 1e6;
    if (event_time_prev != -1)
      delta_event_time = eventTime - event_time_prev;
    if (liveTime_us < Constants::PILEUP_LIVETIME_THRESHOLD_US &&
        delta_event_time == 0)
      continue;

    outEnergy = gain_matched_energy;
    outX = x;
    outY = y;
    outZ = z;
    outEventTime = eventTime;
    outLiveTime = liveTime;

    filteredTrees[crystal]->Fill();
    filteredLiveTime_s[crystal] += liveTime * Constants::TENS_OF_NS_TO_S;

    if (i != 0) {
      event_time_prev = eventTime;
    } else {
      event_time_prev = -1;
    }
  }

  {
    std::lock_guard<std::mutex> lock(print_mutex);
    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      std::cout << filename << " crystal " << c << ": "
                << filteredTrees[c]->GetEntries() << " filtered" << std::endl;
    }
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    filteredTrees[c]->Write("", TObject::kOverwrite);
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TParameter<Float_t> *ltParam = new TParameter<Float_t>(
        Form("LiveTime_Filtered_Crystal%d_s", c), filteredLiveTime_s[c]);
    ltParam->Write("", TObject::kOverwrite);
  }

  BuildTotalEnergyTree(file);

  file->Close();
  return kTRUE;
}

void FilterEvents(std::vector<TString> filenames) {
  Int_t n_files = Int_t(filenames.size());
  Int_t n_workers = Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Filtering " << n_files << " files with " << n_workers
            << " workers." << std::endl;

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);

    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(
          std::async(std::launch::async, FilterFile, filenames[j]));
    }

    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      if (!result) {
        std::cerr << "FAILED: " << filenames[i + j] << std::endl;
      }
    }
  }
}

void Filter() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  ROOT::EnableThreadSafety();

  LoadPixelGainOrDie();

  std::vector<TString> filenames;
  filenames.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  filenames.push_back(Constants::POSTREACTOR_AM241_BA133_20260116);
  filenames.push_back(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116);
  FilterDemo(filenames);

  FilterEvents(Constants::ALL_DATASETS);
}
