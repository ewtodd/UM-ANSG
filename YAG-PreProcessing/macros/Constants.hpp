#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <vector>

namespace Constants {

const PlotSaveFormat SAVE_FORMAT = PlotSaveFormat::kPNG;

// DT5742: 12-bit digitizer
const Int_t ADC_MAX = 4095;

const TString CS137_RAW =
    "/home/e-work/LabData/ANSG/CeYAG/21-10-2025/Cs137-Long/TR_0_0.dat";
const TString AM241_RAW =
    "/home/e-work/LabData/ANSG/CeYAG/21-10-2025/Am241-Weekend/TR_0_0.dat";

const std::vector<TString> ALL_RAW_FILEPATHS = {CS137_RAW, AM241_RAW};

const TString CS137 = "Cs137";
const TString CS137_LABEL = "Cs-137";

const TString AM241 = "Am241";
const TString AM241_LABEL = "Am-241";

const std::vector<TString> ALL_OUTPUT_NAMES = {CS137, AM241};
const std::vector<TString> ALL_FORMATTED_LABELS = {CS137_LABEL, AM241_LABEL};

const TString CS137_FILEPATH = "root_files/" + CS137 + "_raw.root";
const TString AM241_FILEPATH = "root_files/" + AM241 + "_raw.root";

const std::vector<TString> ALL_FILEPATHS = {CS137_FILEPATH, AM241_FILEPATH};

const FileProcessingConfig DEFAULT_PROCESSING_CONFIG = {
    .polarity = -1,
    .trigger_threshold = 0.3,
    .num_samples_baseline = 50,
    .pre_samples = 100,
    .post_samples = 800,
    .sample_waveforms_to_save = 5,
    .max_events = -1,
    .verbose = kTRUE,
    .store_waveforms = kTRUE,
    .input_format = InputFormat::kWAVEDUMP,
};

const PulseIntegralSettings DEFAULT_PULSE_INTEGRAL_SETTINGS = {
    .pre_gate = 20,
    .short_gate = 10,
    .long_gate = 1000,
};

} // namespace Constants

#endif
