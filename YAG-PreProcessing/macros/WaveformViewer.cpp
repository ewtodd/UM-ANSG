#include "Constants.hpp"
#include "InitUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <iostream>
#include <vector>

void WaveformViewer() {
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT);

  PulseIntegralProcessor processor(Constants::DEFAULT_PROCESSING_CONFIG,
                                   Constants::DEFAULT_PULSE_INTEGRAL_SETTINGS);
  processor.ProcessFilesParallel(Constants::ALL_FILEPATHS,
                                 Constants::ALL_OUTPUT_NAMES);
}
