// The dissertation's table of what the pixel gain matching does to the CZT
// result: each run's Ge energy with no pixel gains (CalibrateAndFit.cpp with
// PIXEL_GAINS=off), with the 558 keV gains alone (PIXEL_MAP=off), and with
// the 558 keV gains and the day pixel map, as the analysis uses them (the
// dissertation calls the day map the local gain). Each run is quoted as the
// analysis quotes it (CombineGeResult.cpp): at the midpoint of the gain and
// offset transfer models, so each stage reads its results/ge_<tag>.result and
// its TRANSFER=offset partner's. Writes the table body into the dissertation,
// like the pole-zero chapter's spread_table.tex.
//
// Usage: root -b -q 'GainMatchingTable.cpp+'
#include "ProjectPaths.hpp"

#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace GainMatchingRows {

const Int_t N_STAGES = 3;
const char *const STAGE_TAGS[N_STAGES] = {
    "tau20_sigpb_esc_nogain", "tau20_sigpb_esc_pxoff", "tau20_sigpb_esc_pxday"};
const char *const OFFSET_STAGE_TAGS[N_STAGES] = {"tau20_sigpb_esc_oxfer_nogain",
                                                 "tau20_sigpb_esc_oxfer_pxoff",
                                                 "tau20_sigpb_esc_oxfer_pxday"};
const char *const STAGE_HEADERS[N_STAGES] = {"No pixel gains", "558~keV gains",
                                             "558~keV + local gains"};

// The rows in the order CalibrateAndFit.cpp writes them, each checked against
// its label so a reordering cannot shift a row.
const Int_t N_RUNS = 4;
const char *const RUN_LABELS[N_RUNS] = {
    "Cu_Shield_Signal_10%_(01/13)", "Cd_Shield_Signal_10%_(01/13)",
    "Cd_Shield_Signal_25%_(01/13)", "No_Shield_Signal_5%_(01/15)_[direct]"};
const char *const RUN_NAMES[N_RUNS] = {
    "Cu shield, 10\\%, day~1", "Cd shield, 10\\%, day~1",
    "Cd shield, 25\\%, day~1", "No shield, 5\\%, day~2"};
const Int_t N_INSITU = 3;

const TString TABLE_PATH = TString(gSystem->HomeDirectory()) +
                           "/Dissertation/chapters/chapter_ge_recoils/figures/"
                           "czt_gain_matching_table.tex";

struct RunResult {
  TString label;
  Double_t ge_kev = 0, chi2_ndf = 0;
};

// The rows of one result file: method ge_mu fit_err cal_err bkg_err chi2 ...
// label, '#' lines skipped.
Bool_t ReadResults(const TString &path, std::vector<RunResult> &rows) {
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return kFALSE;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    std::string method, label;
    Double_t ge = 0, fit_err = 0, cal_err = 0, bkg_err = 0, chi2 = 0;
    Double_t gain_err = 0, cal_off = 0, sys_err = 0, anchor_err = 0;
    fields >> method >> ge >> fit_err >> cal_err >> bkg_err >> chi2 >>
        gain_err >> cal_off >> sys_err >> anchor_err >> label;
    RunResult row;
    row.label = label.c_str();
    row.ge_kev = ge;
    row.chi2_ndf = chi2;
    rows.push_back(row);
  }
  return kTRUE;
}

// One result file's rows, each checked against its run's label.
Bool_t ReadStage(const TString &path, std::vector<RunResult> &rows) {
  if (!ReadResults(path, rows))
    return kFALSE;
  if ((Int_t)rows.size() != N_RUNS) {
    std::cerr << "ERROR: " << path << " has " << rows.size()
              << " rows, expected " << N_RUNS << std::endl;
    return kFALSE;
  }
  for (Int_t r = 0; r < N_RUNS; r++)
    if (rows[r].label != RUN_LABELS[r]) {
      std::cerr << "ERROR: row " << r << " of " << path << " is "
                << rows[r].label << ", expected " << RUN_LABELS[r] << std::endl;
      return kFALSE;
    }
  return kTRUE;
}

// The rms of the first n_rows runs' values about their mean, in eV: the
// three day 1 in-situ runs, or every run.
Double_t SpreadEv(const std::vector<RunResult> &rows, Int_t n_rows) {
  Double_t mean = 0;
  for (Int_t r = 0; r < n_rows; r++)
    mean += rows[r].ge_kev / n_rows;
  Double_t sum_sq = 0;
  for (Int_t r = 0; r < n_rows; r++)
    sum_sq += TMath::Power(rows[r].ge_kev - mean, 2);
  return 1000.0 * TMath::Sqrt(sum_sq / n_rows);
}

} // namespace GainMatchingRows

void GainMatchingTable() {
  using namespace GainMatchingRows;
  const TString results_dir = Paths::ProjectRootOf(__FILE__) + "/results";
  std::vector<std::vector<RunResult>> stages(N_STAGES);
  for (Int_t s = 0; s < N_STAGES; s++) {
    std::vector<RunResult> offset;
    if (!ReadStage(results_dir + "/ge_" + STAGE_TAGS[s] + ".result",
                   stages[s]) ||
        !ReadStage(results_dir + "/ge_" + OFFSET_STAGE_TAGS[s] + ".result",
                   offset))
      return;
    // Ge at the midpoint of the two transfer models; the fit's own error
    // and chi2/ndf are the gain model's, production's.
    for (Int_t r = 0; r < N_RUNS; r++)
      stages[s][r].ge_kev = 0.5 * (stages[s][r].ge_kev + offset[r].ge_kev);
  }

  std::ofstream out(TABLE_PATH.Data());
  if (!out) {
    std::cerr << "ERROR: cannot write " << TABLE_PATH << std::endl;
    return;
  }
  out << "% Written by UM-ANSG 73mGe/czt/GainMatchingTable.cpp; edit that, "
         "not this file."
      << std::endl;
  out << "% For the caption: Ge energies in keV, each at the midpoint of the "
         "gain and offset transfer models as the analysis quotes them; the "
         "spreads are the rms about their mean of the three day 1 in-situ "
         "values and of all four runs. The fits' chi2/ndf move without a "
         "consistent trend between the stages, so they are printed by the "
         "macro for the text, not tabulated."
      << std::endl;
  // Sized like the pole-zero chapter's spread_table.tex.
  out << "\\small" << std::endl;
  out << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  out << "\\begin{tabular}{@{}lccc@{}}" << std::endl;
  out << "\\toprule" << std::endl;
  out << "Run";
  for (Int_t s = 0; s < N_STAGES; s++)
    out << " & " << STAGE_HEADERS[s];
  out << " \\\\" << std::endl;
  out << "\\midrule" << std::endl;
  for (Int_t r = 0; r < N_RUNS; r++) {
    out << RUN_NAMES[r];
    for (Int_t s = 0; s < N_STAGES; s++)
      out << Form(" & %.4f", stages[s][r].ge_kev);
    out << " \\\\" << std::endl;
  }
  out << "\\midrule" << std::endl;
  out << "Spread of the day~1 runs";
  for (Int_t s = 0; s < N_STAGES; s++)
    out << Form(" & %.0f~eV", SpreadEv(stages[s], N_INSITU));
  out << " \\\\" << std::endl;
  out << "Spread of all runs";
  for (Int_t s = 0; s < N_STAGES; s++)
    out << Form(" & %.0f~eV", SpreadEv(stages[s], N_RUNS));
  out << " \\\\" << std::endl;
  out << "\\bottomrule" << std::endl;
  out << "\\end{tabular}" << std::endl;
  out.close();

  // For the text: the range of the fits' chi2/ndf over every run and stage.
  Double_t chi2_min = 1e30, chi2_max = -1e30;
  for (Int_t s = 0; s < N_STAGES; s++)
    for (Int_t r = 0; r < N_RUNS; r++) {
      chi2_min = TMath::Min(chi2_min, stages[s][r].chi2_ndf);
      chi2_max = TMath::Max(chi2_max, stages[s][r].chi2_ndf);
    }
  std::cout << Form("chi2/ndf over every run and stage: %.2f to %.2f", chi2_min,
                    chi2_max)
            << std::endl;
  std::cout << "Wrote " << TABLE_PATH << std::endl;
  std::ifstream written(TABLE_PATH.Data());
  std::string line;
  while (std::getline(written, line))
    std::cout << line << std::endl;
}
