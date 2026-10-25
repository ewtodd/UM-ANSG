#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <vector>

namespace Constants {

const Int_t ADC_MAX = 16384;

// Pulse Height histogram (ADC units)
const Int_t PH_BIN_WIDTH = 20;
const Int_t PH_HIST_XMIN = 0;
const Int_t PH_HIST_XMAX = ADC_MAX;
const Int_t PH_HIST_NBINS = (PH_HIST_XMAX - PH_HIST_XMIN) / PH_BIN_WIDTH;

// Pulse Integral histogram (a.u.)
const Int_t PI_BIN_WIDTH = 800;
const Int_t PI_HIST_XMIN = 0;
const Int_t PI_HIST_XMAX = 400000;
const Int_t PI_HIST_NBINS = (PI_HIST_XMAX - PI_HIST_XMIN) / PI_BIN_WIDTH;

// Light Output histogram (keVee)
const Int_t LO_BIN_WIDTH = 1;
const Float_t LO_HIST_XMIN = 0;
const Float_t LO_HIST_XMAX = 1250;
const Int_t LO_HIST_NBINS = (LO_HIST_XMAX - LO_HIST_XMIN) / LO_BIN_WIDTH;

const Float_t E_LA_33KEV = 33.4;
const Float_t E_AM241_59KEV = 59.5409;
const Float_t E_EU152_122KEV = 121.8;
const Float_t E_EU152_245KEV = 244.7;
const Float_t E_EU152_344KEV = 344.3;

const TString CALIBRATION_AM241 = "calibration_Am241";
const TString CALIBRATION_EU152 = "calibration_Eu152";
const TString BACKGROUND = "background";
const TString IRRADIATION_ONE = "irradiation_one";
const TString IRRADIATION_TWO = "irradiation_two";
const TString IRRADIATION_THREE = "irradiation_three";
const TString IRRADIATION_FOUR = "irradiation_four";
const std::vector<TString> ALL_OUTPUT_NAMES = {
    CALIBRATION_AM241, CALIBRATION_EU152, BACKGROUND,      IRRADIATION_ONE,
    IRRADIATION_TWO,   IRRADIATION_THREE, IRRADIATION_FOUR};
const std::vector<TString> IRRADIATION_DATASETS = {
    IRRADIATION_ONE, IRRADIATION_TWO, IRRADIATION_THREE, IRRADIATION_FOUR};

const TString CALIBRATION_AM241_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/59_5keV_calibration_300s/RAW/DataR_CH0@DT5730S_31017_59_5keV_calibration_300s.root";
const TString CALIBRATION_EU152_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/Europium_calibration_300s/RAW/DataR_CH0@DT5730S_31017_Europium_calibration_300s.root";
const TString BACKGROUND_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/day2/bkg_day2/RAW/DataR_CH0@DT5730B_969_bkg_day2.root";
const TString IRRADIATION_ONE_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_1/RAW/DataR_CH0@DT5730S_31017_irradiation_1.root";
const TString IRRADIATION_TWO_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_2/RAW/DataR_CH0@DT5730S_31017_irradiation_2.root";
const TString IRRADIATION_THREE_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_3/RAW/DataR_CH0@DT5730S_31017_irradiation_3.root";
const TString IRRADIATION_FOUR_FILEPATH =
    "/home/e-work/LabData/ANSG/78mBr/day2/irradiation_day2/RAW/DataR_CH0@DT5730B_969_irradiation_day2.root";
const std::vector<TString> ALL_FILEPATHS = {
    CALIBRATION_AM241_FILEPATH, CALIBRATION_EU152_FILEPATH,
    BACKGROUND_FILEPATH,        IRRADIATION_ONE_FILEPATH,
    IRRADIATION_TWO_FILEPATH,   IRRADIATION_THREE_FILEPATH,
    IRRADIATION_FOUR_FILEPATH};

const FileProcessingConfig DEFAULT_PROCESSING_CONFIG = {
    .polarity = -1,
    .trigger_threshold = 0.15,
    .num_samples_baseline = 10,
    .pre_samples = 18,
    .post_samples = 135,
    .sample_waveforms_to_save = 5,
    .max_events = -1,
    .verbose = kTRUE,
    .store_waveforms = kTRUE,
};

const PulseIntegralSettings DEFAULT_PULSE_INTEGRAL_SETTINGS = {
    .pre_gate = 5,
    .short_gate = 40,
    .long_gate = 220,
};
} // namespace Constants
#endif
