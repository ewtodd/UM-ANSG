// The PSI DT5780 CdTe waveforms (100 MS/s, positive pulses on a ~1000 ADC
// baseline, preamplifier decay 325 samples) to one tree of features per run,
// through Analysis-Utilities' PulseHeightProcessor: baseline from the leading
// samples, trigger at a fraction of the peak, pole-zero deconvolution to a
// step, a trapezoid of rise 200 and flat top 100 samples read late on its
// flat top (`trapezoid_height`), and the step's 10-90% rise (`rise_time`, in
// samples). The raw files are only read; the trees go to CdTe::FEATURES_DIR
// under the raw files' names, with a `.stats` summary each.
//
// Usage:
//   root -b -q 'ProcessWaveforms.cpp+("PSI_Beam_Ge_20251011_1715_1")'  one run
//   root -b -q 'ProcessWaveforms.cpp+'  every CdTe run, in parallel
#include "CdTeAnalysis.hpp"
#include "WaveformProcessingUtils.hpp"

#include <TList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

namespace CdTeProcessing {

// The trigger sits near sample 690 (1900 in the LongWf runs), about 1% of
// pulses earlier; a pulse that triggers at pre_samples starts a few samples
// before it, so the baseline ends 50 samples short of that.
FileProcessingConfig Config() {
  FileProcessingConfig config;
  config.polarity = 1;
  config.trigger_threshold = 0.15;
  config.num_samples_baseline = 350;
  config.pre_samples = 400;
  config.post_samples = 400;
  config.sample_waveforms_to_save = 0;
  config.verbose = kTRUE;
  config.store_waveforms = kFALSE;
  config.input_format = InputFormat::kCOMPASS;
  return config;
}

// The trigger is the 15% crossing, a few samples into the rise (up to ~10 for
// the slowest pulses used), so the height is read at 90% of the flat top:
// past the rise of the slowest step and before the fall.
PulseHeightSettings Filter() {
  PulseHeightSettings filter;
  filter.decay_time_samples = 325.0;
  filter.rise_samples = 200;
  filter.flat_top_samples = 100;
  filter.peaking_fraction = 0.9;
  filter.peak_average_samples = 1;
  filter.rise_plateau_start = 150;
  filter.rise_plateau_samples = 100;
  filter.rise_search_back = 50;
  return filter;
}

// The CdTe runs in the raw directory (DT5780 channel 1, PSI_ names).
std::vector<TString> RawRuns() {
  std::vector<TString> runs;
  const TString prefix = CdTe::FILE_PREFIX;
  TSystemDirectory dir(CdTe::RAW_DIR, CdTe::RAW_DIR);
  TList *files = dir.GetListOfFiles();
  if (!files)
    return runs;
  TIter next(files);
  while (TSystemFile *file = static_cast<TSystemFile *>(next())) {
    TString name = file->GetName();
    if (!name.BeginsWith(prefix + "PSI_") || !name.EndsWith(".root"))
      continue;
    name.Remove(0, prefix.Length());
    name.Remove(name.Length() - 5);
    runs.push_back(name);
  }
  delete files;
  return runs;
}

} // namespace CdTeProcessing

void ProcessWaveforms(const char *run = "") {
  using namespace CdTeProcessing;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/cdte/waveforms",
                                CdTe::FEATURES_DIR);
  PulseHeightProcessor processor(Config(), Filter());
  std::vector<TString> runs;
  if (TString(run).Length() > 0)
    runs.push_back(run);
  else
    runs = RawRuns();
  std::vector<TString> inputs, outputs;
  for (size_t i = 0; i < runs.size(); i++) {
    inputs.push_back(TString(CdTe::RAW_DIR) + CdTe::FILE_PREFIX + runs[i] +
                     ".root");
    outputs.push_back(TString(CdTe::FILE_PREFIX) + runs[i]);
  }
  processor.ProcessFilesParallel(inputs, outputs, 8);
}
