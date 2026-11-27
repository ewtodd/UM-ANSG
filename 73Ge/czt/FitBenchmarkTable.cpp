// The dissertation appendix's table of fit times on each RooFit backend, from
// FitBenchmark.cpp's results/fit_benchmark_<backend>_<cores>.result: the CPU
// build on one core, the CPU build on every core and the CUDA build, each
// fit's mean wall time and its speed-up over the CPU on one core. It prints
// what the caption quotes: the events fitted, the hardware, the cores the CPU
// fits kept busy and how far Ge moves between backends.
//
// Usage: root -b -q 'FitBenchmarkTable.cpp+'
#include "ProjectPaths.hpp"

#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const TString TABLE_PATH =
    TString(gSystem->HomeDirectory()) +
    "/Dissertation/appendices/appendix_au/fit_benchmark_table.tex";

// One result file: each line a key and the rest of the line its value.
struct BenchmarkRun {
  TString backend;
  Int_t cores = 0;
  std::map<std::string, std::string> values;
  Double_t Number(const std::string &key) const {
    std::map<std::string, std::string>::const_iterator it = values.find(key);
    return it == values.end() ? 0.0 : atof(it->second.c_str());
  }
  std::string Text(const std::string &key) const {
    std::map<std::string, std::string>::const_iterator it = values.find(key);
    return it == values.end() ? std::string("unknown") : it->second;
  }
};

BenchmarkRun ReadRun(const TString &path) {
  BenchmarkRun run;
  std::ifstream in(path.Data());
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    size_t space = line.find(' ');
    if (space == std::string::npos)
      continue;
    run.values[line.substr(0, space)] = line.substr(space + 1);
  }
  run.backend = run.Text("backend").c_str();
  run.cores = (Int_t)run.Number("cores");
  return run;
}

// CPU runs by their cores, fewest first, then the CUDA runs.
Bool_t TableOrder(const BenchmarkRun &a, const BenchmarkRun &b) {
  if (a.backend != b.backend)
    return a.backend == "cpu";
  return a.cores < b.cores;
}

TString RowLabel(const BenchmarkRun &run) {
  if (run.backend == "cuda")
    return "CUDA (GPU)";
  return TString::Format("CPU, %d thread%s", run.cores,
                         run.cores == 1 ? "" : "s");
}

} // namespace

void FitBenchmarkTable() {
  const TString results_dir = Paths::ProjectRootOf(__FILE__) + "/results";
  std::vector<BenchmarkRun> runs;
  void *dir = gSystem->OpenDirectory(results_dir);
  if (!dir) {
    std::cerr << "ERROR: cannot open " << results_dir << std::endl;
    return;
  }
  const char *entry = nullptr;
  while ((entry = gSystem->GetDirEntry(dir)) != nullptr) {
    TString name = entry;
    if (name.BeginsWith("fit_benchmark_") && name.EndsWith(".result"))
      runs.push_back(ReadRun(results_dir + "/" + name));
  }
  gSystem->FreeDirectory(dir);
  std::sort(runs.begin(), runs.end(), TableOrder);
  if (runs.empty() || runs[0].backend != "cpu") {
    std::cerr << "ERROR: no CPU benchmark to compare with; run FitBenchmark.cpp"
              << std::endl;
    return;
  }
  const BenchmarkRun &reference = runs[0];
  const char *fits[2] = {"direct", "insitu"};

  std::ofstream out(TABLE_PATH.Data());
  if (!out.is_open()) {
    std::cerr << "ERROR: cannot write " << TABLE_PATH << std::endl;
    return;
  }
  out << "% Written by UM-ANSG 73Ge/czt/FitBenchmarkTable.cpp from "
         "FitBenchmark.cpp's results."
      << std::endl;
  out << "\\begin{tabular}{@{}lcccc@{}}" << std::endl;
  out << "\t\\toprule" << std::endl;
  out << "\t& \\multicolumn{2}{c}{Direct fit, day~2} & "
         "\\multicolumn{2}{c}{Simultaneous fit, day~1} \\\\"
      << std::endl;
  out << "\t\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}" << std::endl;
  out << "\tBackend & Time [s] & Speed-up & Time [s] & Speed-up \\\\"
      << std::endl;
  out << "\t\\midrule" << std::endl;
  for (size_t r = 0; r < runs.size(); r++) {
    out << "\t" << RowLabel(runs[r]);
    for (Int_t f = 0; f < 2; f++) {
      const std::string key = std::string(fits[f]) + "_wall_mean";
      const Double_t time = runs[r].Number(key);
      const Double_t speedup = time > 0 ? reference.Number(key) / time : 0;
      out << " & " << std::fixed << std::setprecision(1) << time << " & "
          << speedup;
    }
    out << " \\\\" << std::endl;
  }
  out << "\t\\bottomrule" << std::endl;
  out << "\\end{tabular}" << std::endl;
  out.close();
  std::cout << "Wrote " << TABLE_PATH << std::endl;

  std::cout << "Events fitted: direct " << reference.Text("direct_events")
            << ", simultaneous " << reference.Text("insitu_sig_events")
            << " (signal) + " << reference.Text("insitu_bkg_events")
            << " (background)" << std::endl;
  for (size_t r = 0; r < runs.size(); r++) {
    std::cout << RowLabel(runs[r]) << ": CPU " << runs[r].Text("cpu_model")
              << ", GPU " << runs[r].Text("gpu_model");
    for (Int_t f = 0; f < 2; f++) {
      const std::string fit = fits[f];
      const Double_t shift_ev = 1000.0 * (runs[r].Number(fit + "_ge_kev") -
                                          reference.Number(fit + "_ge_kev"));
      std::cout << "; " << fit << " " << std::fixed << std::setprecision(1)
                << runs[r].Number(fit + "_wall_min") << "-"
                << runs[r].Number(fit + "_wall_max") << " s, "
                << runs[r].Number(fit + "_cores_busy") << " cores busy, Ge "
                << std::setprecision(5) << runs[r].Number(fit + "_ge_kev")
                << " keV (" << std::showpos << std::setprecision(3) << shift_ev
                << std::noshowpos << " eV from the reference)"
                << (runs[r].Number(fit + "_converged") > 0 ? "" : " [FAILED]");
    }
    std::cout << std::endl;
  }
}
