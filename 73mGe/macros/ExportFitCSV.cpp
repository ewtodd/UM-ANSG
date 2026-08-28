// Thin entry point: export the Ge-73m spectrum + fit + residuals to CSV for the
// collaborator's plotting. All logic lives in CalibrationLow.cpp's
// RunFitCSVExport() so the exported curve is byte-for-byte the published
// postcal fit (it loads the cached converged fit; no refit).
//
//   root -l ExportFitCSV.cpp
//
// Writes results/spectrum_CuShield10_20260113_{sig,bkg}.csv :
//   sig = Pb-Ka1, Pb-Ka2, Ge(68.75) spectrum (the figure of interest)
//   bkg = Pb-only background channel
// Each file: a smooth fit-curve block, then a binned data + residual-pull
// block.
#include "CalibrationLow.cpp"

void ExportFitCSV() { RunFitCSVExport(); }
