#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>

namespace Paths {
inline TString ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString macros_dir = gSystem->DirName(path);
  return gSystem->DirName(macros_dir);
}
} // namespace Paths

namespace Constants {

const PlotSaveFormat SAVE_FORMAT = PlotSaveFormat::kPDF;

// 14 bit digitizer
const Int_t ADC_MAX = 16384;

// Pulse height histogram (ADC units)
const Int_t PH_BIN_WIDTH = 20;
const Int_t PH_HIST_XMIN = 0;
const Int_t PH_HIST_XMAX = ADC_MAX;
const Int_t PH_HIST_NBINS = (PH_HIST_XMAX - PH_HIST_XMIN) / PH_BIN_WIDTH;

// Pulse integral histogram (a.u.)
const Int_t PI_BIN_WIDTH = 300;
const Int_t PI_HIST_XMIN = 0;
const Int_t PI_HIST_XMAX = 120000;
const Int_t PI_HIST_NBINS = (PI_HIST_XMAX - PI_HIST_XMIN) / PI_BIN_WIDTH;

// Light output histogram (keVee)
const Int_t LO_BIN_WIDTH = 10;
const Float_t LO_HIST_XMIN = 0;
const Float_t LO_HIST_XMAX = 2000;
const Int_t LO_HIST_NBINS = (LO_HIST_XMAX - LO_HIST_XMIN) / LO_BIN_WIDTH;

// Charge comparison histogram (a.u.)
const Int_t CC_HIST_NBINS = 250;
const Float_t CC_HIST_XMIN = 0.0;
const Float_t CC_HIST_XMAX = 1.0;

// Shape indicator histogram (a.u.)
const Int_t SI_HIST_NBINS = 1000;
const Float_t SI_HIST_XMIN = -0.5;
const Float_t SI_HIST_XMAX = 0.5;

const Float_t E_AM241_59KEV = 59.5409;
const Float_t E_CS137_662KEV = 661.7;
const Float_t E_NA22_511KEV = 511.0;
const Float_t E_NA22_1274KEV = 1274.5;
const Float_t E_CO60_1332KEV = 1332.5;

const TString AM241 = "Am241";
const TString CS137 = "Cs137";
const TString NA22 = "Na22";
const TString CO60 = "Co60";
const TString AM241_AND_CS137 = "Am241Cs137";
const TString AM241_AND_NA22 = "Am241Na22";
const TString AM241_AND_CO60 = "Am241Co60";

const std::vector<TString> ALL_OUTPUT_NAMES = {
    AM241, CS137, NA22, CO60, AM241_AND_CS137, AM241_AND_NA22, AM241_AND_CO60};

const std::vector<TString> SINGLE_OUTPUT_NAMES = {AM241, CS137, NA22, CO60};

const TString AM241_FILENAME = "DataR_Am241.root";
const TString CS137_FILENAME = "DataR_Cs137.root";
const TString NA22_FILENAME = "DataR_Na22.root";
const TString CO60_FILENAME = "DataR_Co60.root";
const TString AM241_AND_CS137_FILENAME = "DataR_Am241Cs137.root";
const TString AM241_AND_NA22_FILENAME = "DataR_Am241Na22.root";
const TString AM241_AND_CO60_FILENAME = "DataR_Am241Co60.root";

const std::vector<TString> ALL_FILENAMES = {
    AM241_FILENAME,         CS137_FILENAME,           NA22_FILENAME,
    CO60_FILENAME,          AM241_AND_CS137_FILENAME, AM241_AND_NA22_FILENAME,
    AM241_AND_CO60_FILENAME};

const std::vector<TString> SINGLE_FILENAMES = {AM241_FILENAME, CS137_FILENAME,
                                               NA22_FILENAME, CO60_FILENAME};

const TString AM241_LABEL = "Am-241";
const TString CS137_LABEL = "Cs-137";
const TString NA22_LABEL = "Na-22";
const TString CO60_LABEL = "Co-60";
const TString AM241_AND_CS137_LABEL = "Am-241 & Cs-137";
const TString AM241_AND_NA22_LABEL = "Am-241 & Na-22";
const TString AM241_AND_CO60_LABEL = "Am-241 & Co-60";

const std::vector<TString> ALL_FORMATTED_LABELS = {
    AM241_LABEL,         CS137_LABEL,           NA22_LABEL,
    CO60_LABEL,          AM241_AND_CS137_LABEL, AM241_AND_NA22_LABEL,
    AM241_AND_CO60_LABEL};

const std::vector<TString> SINGLE_FORMATTED_LABELS = {AM241_LABEL, CS137_LABEL,
                                                      NA22_LABEL, CO60_LABEL};

const FileProcessingConfig DEFAULT_PROCESSING_CONFIG = {
    .polarity = -1,
    .trigger_threshold = 0.2,
    .num_samples_baseline = 15,
    .pre_samples = 20,
    .post_samples = 200,
    .pre_gate = 5,
    .short_gate = 10,
    .long_gate = 250,
    .sample_waveforms_to_save = 5,
    .max_events = -1,
    .verbose = kTRUE,
    .store_waveforms = kTRUE,
};

const Int_t OPTIMAL_SHORT_GATE = 17;
const Int_t OPTIMAL_LONG_GATE = 150;

} // namespace Constants
#endif
