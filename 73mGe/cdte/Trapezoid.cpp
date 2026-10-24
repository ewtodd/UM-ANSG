// Offline energy reconstruction for the PSI DT5780 CdTe waveforms (100 MS/s,
// positive pulses on a ~1000 ADC baseline, exponential preamp decay with
// tau = 325 samples). Per waveform:
//   - leading-edge trigger 100 ADC above the first 100 samples' mean;
//   - baseline: mean of up to 500 samples ending 20 before the trigger;
//   - pole-zero deconvolution to a step, y[n] = x[n] + (1/tau) sum_{j<n} x[j];
//   - energy = mean of y over k samples starting gap samples after the
//     trigger minus the mean over the k samples before it (a trapezoid of
//     rise k and flat top gap), for four (k, gap) settings; energy_x the same
//     on a baseline that includes an earlier pulse's decaying tail;
//   - the step's 10-90% rise time, its slope over the late window (zero when
//     tau is right), and the number of extra rising edges (pileup).
// Writes one tree "trap" per raw file into CdTe::TRAP_DIR, under the raw
// file's name; the raw files are only read.
//
// Usage:
//   root -b -q 'Trapezoid.cpp+("PSI_Beam_Ge_20251011_1715_1")'  one run
//   root -b -q 'Trapezoid.cpp+'          every CdTe run not yet processed
#include "CdTeAnalysis.hpp"

#include <TArrayS.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>

namespace CdTeTrapezoid {

const Double_t PREAMP_TAU_SAMPLES = 325.0;
const Int_t N_SHAPES = 4;
const Int_t SHAPE_K[N_SHAPES] = {100, 200, 300, 400};
const Int_t SHAPE_GAP[N_SHAPES] = {60, 100, 150, 200};

// Pole-zero deconvolution of the baseline-subtracted waveform to a step.
void Deconvolve(const TArrayS *samples, Double_t base, Double_t tau,
                std::vector<Double_t> &y) {
  Int_t n = samples->GetSize();
  y.assign(n, 0.0);
  Double_t running = 0;
  for (Int_t k = 0; k < n; k++) {
    Double_t x = samples->At(k) - base;
    y[k] = x + running / tau;
    running += x;
  }
}

// Step height between the k samples before the trigger and the k samples
// starting gap after it, for each shaping.
void TrapEnergies(const std::vector<Double_t> &y, Int_t trig, Int_t n,
                  Float_t *out) {
  for (Int_t s = 0; s < N_SHAPES; s++) {
    out[s] = -1;
    Int_t kk = SHAPE_K[s], gap = SHAPE_GAP[s];
    if (trig - kk < 0 || trig + gap + kk > n)
      continue;
    Double_t before = 0, after = 0;
    for (Int_t k = trig - kk; k < trig; k++)
      before += y[k];
    for (Int_t k = trig + gap; k < trig + gap + kk; k++)
      after += y[k];
    out[s] = (after - before) / kk;
  }
}

void ProcessFile(const TString &input, const TString &output, Double_t tau) {
  TFile *in = TFile::Open(input);
  if (!in || in->IsZombie()) {
    std::cerr << "ERROR: cannot open " << input << std::endl;
    delete in;
    return;
  }
  TTree *tree = static_cast<TTree *>(in->Get("Data_R"));
  if (!tree) {
    std::cerr << "ERROR: no Data_R tree in " << input << std::endl;
    in->Close();
    delete in;
    return;
  }
  TArrayS *samples = nullptr;
  UShort_t dpp = 0;
  ULong64_t timestamp = 0;
  UInt_t flags = 0;
  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("Energy", &dpp);
  tree->SetBranchAddress("Timestamp", &timestamp);
  tree->SetBranchAddress("Flags", &flags);

  TFile out(output, "RECREATE");
  TTree *t = new TTree("trap", "PSI offline trapezoid");
  Float_t energy[N_SHAPES], energy_x[N_SHAPES], baseline = 0, baseline_rms = 0,
                                                rise = 0, slope = 0;
  Float_t baseline_exp = 0, tail_amp = 0;
  Short_t trigger = 0, n_edges = 0;
  t->Branch("timestamp", &timestamp, "timestamp/l");
  t->Branch("dpp", &dpp, "dpp/s");
  t->Branch("flags", &flags, "flags/i");
  t->Branch("energy", energy, Form("energy[%d]/F", N_SHAPES));
  t->Branch("energy_x", energy_x, Form("energy_x[%d]/F", N_SHAPES));
  t->Branch("baseline", &baseline, "baseline/F");
  t->Branch("baseline_exp", &baseline_exp, "baseline_exp/F");
  t->Branch("tail_amp", &tail_amp, "tail_amp/F");
  t->Branch("baseline_rms", &baseline_rms, "baseline_rms/F");
  t->Branch("trigger", &trigger, "trigger/S");
  t->Branch("rise", &rise, "rise/F");
  t->Branch("slope", &slope, "slope/F");
  t->Branch("n_edges", &n_edges, "n_edges/S");

  std::vector<Double_t> y;
  Long64_t n_events = tree->GetEntries();
  for (Long64_t i = 0; i < n_events; i++) {
    tree->GetEntry(i);
    Int_t n = samples->GetSize();
    for (Int_t s = 0; s < N_SHAPES; s++) {
      energy[s] = -1;
      energy_x[s] = -1;
    }
    baseline_exp = 0;
    tail_amp = 0;
    rise = -1;
    slope = 0;
    n_edges = 0;
    trigger = -1;
    if (n < 400) {
      t->Fill();
      continue;
    }
    Double_t first = 0;
    for (Int_t k = 0; k < 100; k++)
      first += samples->At(k) / 100.0;
    Int_t trig = -1;
    for (Int_t k = 100; k < n; k++)
      if (samples->At(k) - first > 100.0) {
        trig = k;
        break;
      }
    if (trig < 0) {
      t->Fill();
      continue;
    }
    trigger = trig;
    Int_t b_end = TMath::Max(trig - 20, 10);
    Int_t b_start = TMath::Max(b_end - 500, 0);
    Double_t sum = 0, sum2 = 0;
    for (Int_t k = b_start; k < b_end; k++) {
      sum += samples->At(k);
      sum2 += (Double_t)samples->At(k) * samples->At(k);
    }
    Int_t nb = b_end - b_start;
    baseline = sum / nb;
    baseline_rms =
        TMath::Sqrt(TMath::Max(sum2 / nb - baseline * baseline, 0.0));

    // Baseline plus the decaying tail of an earlier pulse, B + C exp(-k/tau),
    // by least squares over the whole pre-trigger region.
    Double_t s0 = 0, s1 = 0, s2 = 0, sy = 0, sey = 0;
    for (Int_t k = 0; k < b_end; k++) {
      Double_t e = TMath::Exp(-k / tau);
      s0 += 1.0;
      s1 += e;
      s2 += e * e;
      sy += samples->At(k);
      sey += e * samples->At(k);
    }
    Double_t det = s0 * s2 - s1 * s1;
    baseline_exp = (det != 0) ? (sy * s2 - s1 * sey) / det : baseline;
    tail_amp = (det != 0) ? (s0 * sey - s1 * sy) / det : 0;

    Deconvolve(samples, baseline_exp, tau, y);
    TrapEnergies(y, trig, n, energy_x);
    Deconvolve(samples, baseline, tau, y);
    TrapEnergies(y, trig, n, energy);
    // Rise time on the step, from its plateau 150-250 samples after trigger.
    if (trig + 250 < n) {
      Double_t plateau = 0;
      for (Int_t k = trig + 150; k < trig + 250; k++)
        plateau += y[k] / 100.0;
      Double_t t10 = -1, t90 = -1;
      for (Int_t k = trig - 20; k < trig + 250 && k < n; k++) {
        if (k < 1)
          continue;
        if (t10 < 0 && y[k] >= 0.1 * plateau)
          t10 = k - 1 + (0.1 * plateau - y[k - 1]) / (y[k] - y[k - 1]);
        if (t90 < 0 && y[k] >= 0.9 * plateau) {
          t90 = k - 1 + (0.9 * plateau - y[k - 1]) / (y[k] - y[k - 1]);
          break;
        }
      }
      if (t10 >= 0 && t90 >= 0)
        rise = t90 - t10;
      // Late-window slope (per sample, in ADC), from a least-squares line.
      Int_t s_lo = trig + 300, s_hi = TMath::Min(n, trig + 900);
      if (s_hi - s_lo > 100) {
        Double_t sx = 0, sy_late = 0, sxx = 0, sxy = 0;
        Int_t m = s_hi - s_lo;
        for (Int_t k = s_lo; k < s_hi; k++) {
          sx += k;
          sy_late += y[k];
          sxx += (Double_t)k * k;
          sxy += k * y[k];
        }
        slope = (m * sxy - sx * sy_late) / (m * sxx - sx * sx);
      }
    }
    // Extra rising edges: a 10-sample difference of the step above 100 ADC,
    // outside the main edge.
    for (Int_t k = 10; k < n; k++) {
      if (k > trig - 20 && k < trig + 60)
        continue;
      if (y[k] - y[k - 10] > 100.0) {
        n_edges++;
        k += 60;
      }
    }
    t->Fill();
  }
  t->Write();
  out.Close();
  in->Close();
  delete in;
  std::cout << "TRAP " << gSystem->BaseName(input) << " processed " << n_events
            << std::endl;
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
    if (!name.BeginsWith(prefix) || !name.EndsWith(".root"))
      continue;
    name.Remove(0, prefix.Length());
    name.Remove(name.Length() - 5);
    runs.push_back(name);
  }
  delete files;
  return runs;
}

} // namespace CdTeTrapezoid

void Trapezoid(const char *run = "") {
  using namespace CdTeTrapezoid;
  gSystem->mkdir(CdTe::TRAP_DIR, kTRUE);
  std::vector<TString> runs;
  if (TString(run).Length() > 0)
    runs.push_back(run);
  else
    runs = RawRuns();
  for (size_t i = 0; i < runs.size(); i++) {
    TString input =
        TString(CdTe::RAW_DIR) + CdTe::FILE_PREFIX + runs[i] + ".root";
    TString output = CdTe::TrapFile(runs[i]);
    if (TString(run).Length() == 0 && !gSystem->AccessPathName(output)) {
      std::cout << "TRAP " << runs[i] << " already processed" << std::endl;
      continue;
    }
    ProcessFile(input, output, PREAMP_TAU_SAMPLES);
  }
}
