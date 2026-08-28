#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>
#include <future>
#include <mutex>
#include <thread>

static std::mutex print_mutex;

Bool_t AddHistogram(TString filename) {
  TFile *file = IO::OpenForWriting("filtered/" + filename + ".root", "UPDATE");

  TString treeType = Constants::USE_FILTERED ? "filtered" : "unfiltered";

  TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                        Form("%s; Energy [keV]; Counts / %d eV",
                             filename.Data(), Constants::BIN_WIDTH_EV),
                        Constants::HIST_NBINS, Constants::HIST_XMIN,
                        Constants::HIST_XMAX);

  TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                              Form("%s; Energy [keV]; Counts / %d eV",
                                   filename.Data(), Constants::BIN_WIDTH_EV),
                              Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                              Constants::ZOOMED_XMAX);

  TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                            Form("%s; Energy [keV]; Counts / %d eV",
                                 filename.Data(), Constants::BIN_WIDTH_EV),
                            Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                            Constants::PEAK_XMAX);

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TString treeName = Form("crystal%d_%s_tree", c, treeType.Data());
    TTree *tree = static_cast<TTree *>(file->Get(treeName));
    if (!tree) {
      std::lock_guard<std::mutex> lock(print_mutex);
      std::cerr << filename << ": missing " << treeName << std::endl;
      continue;
    }

    Float_t energy = 0;
    tree->SetBranchAddress("energykeV", &energy);

    Int_t n_entries = tree->GetEntries();

    TH1F *cHist =
        new TH1F(PlottingUtils::GetRandomName(), "", Constants::HIST_NBINS,
                 Constants::HIST_XMIN, Constants::HIST_XMAX);
    TH1F *cZoomed =
        new TH1F(PlottingUtils::GetRandomName(), "", Constants::ZOOMED_NBINS,
                 Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);
    TH1F *cPeak =
        new TH1F(PlottingUtils::GetRandomName(), "", Constants::PEAK_NBINS,
                 Constants::PEAK_XMIN, Constants::PEAK_XMAX);

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      cHist->Fill(energy);

      if (energy > Constants::ZOOMED_XMIN && energy < Constants::ZOOMED_XMAX)
        cZoomed->Fill(energy);
      if (energy > Constants::PEAK_XMIN && energy < Constants::PEAK_XMAX)
        cPeak->Fill(energy);
    }

    cHist->Write(Form("hist_crystal%d", c), TObject::kOverwrite);
    cZoomed->Write(Form("zoomedHist_crystal%d", c), TObject::kOverwrite);
    cPeak->Write(Form("peakHist_crystal%d", c), TObject::kOverwrite);

    hist->Add(cHist);
    zoomedHist->Add(cZoomed);
    peakHist->Add(cPeak);

    delete cHist;
    delete cZoomed;
    delete cPeak;
  }

  {
    std::lock_guard<std::mutex> lock(print_mutex);
    std::cout << "Created histograms for " << filename << " (using " << treeType
              << " trees)" << std::endl;
  }

  PlottingUtils::ConfigureHistogram(hist, kP10Violet);
  PlottingUtils::ConfigureHistogram(zoomedHist, kP10Violet);
  PlottingUtils::ConfigureHistogram(peakHist, kP10Violet);

  hist->GetYaxis()->SetTitleOffset(1.2);
  zoomedHist->GetYaxis()->SetTitleOffset(1.2);
  peakHist->GetYaxis()->SetTitleOffset(1.2);

  hist->Write("hist", TObject::kOverwrite);
  zoomedHist->Write("zoomedHist", TObject::kOverwrite);
  peakHist->Write("peakHist", TObject::kOverwrite);

  {
    std::lock_guard<std::mutex> lock(print_mutex);
    std::cout << "Wrote histograms for " << filename << std::endl;
  }

  file->Close();
  return kTRUE;
}

void Histogram() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  ROOT::EnableThreadSafety();

  std::vector<TString> filenames = Constants::ALL_DATASETS;
  Int_t n_files = Int_t(filenames.size());
  Int_t n_workers = Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Processing " << n_files << " histograms with " << n_workers
            << " workers." << std::endl;

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);

    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(
          std::async(std::launch::async, AddHistogram, filenames[j]));
    }

    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      if (!result) {
        std::cerr << "FAILED: " << filenames[i + j] << std::endl;
      }
    }
  }

  std::cout << "All histograms complete." << std::endl;
}
