// The PSI DT5780 Si waveforms (100 MS/s, positive pulses on a ~590 ADC
// baseline, preamplifier decay 324 samples) to one tree of features per run,
// through Analysis-Utilities' PulseHeightProcessor with SiAnalysis.hpp's
// settings: baseline from the record's first samples up to 75 before the
// nominal trigger, trigger at 15% of the peak, pole-zero deconvolution to a
// step, a trapezoid read late on its flat top (`trapezoid_height`), and the
// step's 10-90% rise (`rise_time`, in samples). The raw files are only read;
// the trees go to Si::FEATURES_DIR under the raw files' names, with a
// `.stats` summary each. The baseline window depends on the record length,
// so the runs are processed in two groups.
//
// Usage:
//   root -b -q 'ProcessWaveforms.cpp+("PSI_Beam_Ge_20251011_1715_2")'  one run
//   root -b -q 'ProcessWaveforms.cpp+'  every Si run, in parallel
#include "SiAnalysis.hpp"

#include <TList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

namespace SiProcessing {

// The Si runs in the raw directory (DT5780 channel 0, PSI_ names; the
// DataR_ files are another acquisition setup).
std::vector<TString> RawRuns() {
  std::vector<TString> runs;
  const TString prefix = Si::FILE_PREFIX;
  TSystemDirectory dir(Si::RAW_DIR, Si::RAW_DIR);
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
  std::sort(runs.begin(), runs.end());
  return runs;
}

const Int_t WORKERS = 12;

} // namespace SiProcessing

void ProcessWaveforms(const char *run = "") {
  using namespace SiProcessing;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG, project_root + "/plots/studies/psi_si/waveforms",
      Si::FEATURES_DIR);
  std::vector<TString> runs;
  if (TString(run).Length() > 0)
    runs.push_back(run);
  else
    runs = RawRuns();
  const Int_t lengths[2] = {Si::SHORT_RECORD, Si::LONG_RECORD};
  for (Int_t l = 0; l < 2; l++) {
    std::vector<TString> inputs, outputs;
    for (size_t i = 0; i < runs.size(); i++) {
      Int_t length = Si::RecordLength(runs[i]);
      if (length != lengths[l]) {
        if (length != Si::SHORT_RECORD && length != Si::LONG_RECORD)
          std::cout << "SI-PROCESS skipping " << runs[i] << ": record length "
                    << length << std::endl;
        continue;
      }
      inputs.push_back(Si::RawFile(runs[i]));
      outputs.push_back(TString(Si::FILE_PREFIX) + runs[i]);
    }
    if (inputs.empty())
      continue;
    std::cout << "SI-PROCESS " << inputs.size() << " runs of " << lengths[l]
              << "-sample records" << std::endl;
    PulseHeightProcessor processor(Si::Config(lengths[l]), Si::Filter());
    processor.ProcessFilesParallel(inputs, outputs, WORKERS);
  }
}
