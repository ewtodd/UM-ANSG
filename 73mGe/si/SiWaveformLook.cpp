// A look at the PSI Si waveforms (Amptek Si, DT5780 channel 0), through
// Analysis-Utilities' PulseHeightProcessor: source runs against beam runs.
//   stage 1: process a sample of each run with the waveforms stored, the
//            preamplifier decay given (0: no deconvolution);
//   stage 2: the figures: example waveforms with their step and trapezoid,
//            records whose baseline sits at the bottom of the ADC range,
//            the average pulse and its decay constant, the baseline noise and
//            its spectrum, pile-up before and after the pulse, the rise time
//            against pulse height, and the trapezoid against the digitiser's
//            own energy.
//
// Usage: root -b -q 'SiWaveformLook.cpp+(1, 0)' then
//        root -b -q 'SiWaveformLook.cpp+(2)' for the decay constant, then
//        root -b -q 'SiWaveformLook.cpp+(1, <decay>)' and stage 2 again.
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"
#include "WaveformProcessingUtils.hpp"

#include <TArrayF.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TVirtualFFT.h>

#include <iomanip>
#include <iostream>
#include <vector>

namespace SiLook {

const char *const RAW_DIR = "/labdata/PSI/Si/";
const char *const FILE_PREFIX = "Data_CH0@DT5780P_14068_";
const char *const LOOK_DIR = "/labdata/ANSG/PSI/processed/Si/waveform_look";
const Double_t SAMPLE_NS = 10.0;
const Int_t PRE_SAMPLES = 400;
const Int_t POST_SAMPLES = 1000;
const Int_t BASELINE_SAMPLES = 350;
const Int_t EVENTS_PER_RUN = 20000;
// DPP channels per keV, from the earlier list-mode look (Am 59.54 keV near
// channel 1428).
const Double_t DPP_PER_KEV = 23.98;

struct LookRun {
  TString tag, label, run;
  Int_t color;
};

std::vector<LookRun> Runs() {
  std::vector<LookRun> runs;
  runs.push_back({"am", "Am-241 source, 10/10",
                  "PSI_Calibration_CdTe_Am_Si_Am_20251010_1422_1", kBlue + 1});
  runs.push_back({"ba", "Ba-133 source, 10/12",
                  "PSI_Calib_Long_CdTe_Cs_Si_Ba_20251012_0917_1", kCyan + 2});
  runs.push_back({"beam_ge", "Beam, Ge target, 10/11",
                  "PSI_Beam_Ge_20251011_1715_1", kRed + 1});
  runs.push_back({"beam_getarget", "Beam, Ge target, 10/10",
                  "PSI_Beam_GeTarget_20251010_1525_1", kOrange + 7});
  runs.push_back({"beam_nothing", "Beam, nothing in front, 10/13",
                  "PSI_Beam_NoPb_CdTe_Nothing_Si_Nothing_20251013_0652_1",
                  kGreen + 2});
  return runs;
}

void Process(Double_t decay_samples) {
  FileProcessingConfig config;
  config.polarity = 1;
  config.trigger_threshold = 0.15;
  config.num_samples_baseline = BASELINE_SAMPLES;
  config.pre_samples = PRE_SAMPLES;
  config.post_samples = POST_SAMPLES;
  config.sample_waveforms_to_save = 6;
  config.max_events = EVENTS_PER_RUN;
  config.store_waveforms = kTRUE;
  config.input_format = InputFormat::kCOMPASS;
  PulseHeightSettings filter;
  filter.decay_time_samples = decay_samples;
  filter.rise_samples = 200;
  filter.flat_top_samples = 100;
  filter.peaking_fraction = 0.9;
  filter.rise_plateau_start = 150;
  filter.rise_plateau_samples = 100;
  filter.rise_search_back = 50;
  PulseHeightProcessor processor(config, filter);
  std::vector<LookRun> runs = Runs();
  for (size_t i = 0; i < runs.size(); i++)
    processor.ProcessFile(TString(RAW_DIR) + FILE_PREFIX + runs[i].run +
                              ".root",
                          "si_" + runs[i].tag);
}

TString LookFile(const LookRun &run) {
  return TString(LOOK_DIR) + "/si_" + run.tag + ".root";
}

// A second pulse in the record: a 10-sample rise of ten baseline RMS after
// the pulse's own edge.
Bool_t HasEdgeAfterPulse(const TArrayF &samples, Float_t rms) {
  Double_t threshold = 10.0 * TMath::Max(rms, Float_t(1));
  for (Int_t k = PRE_SAMPLES + 60; k < samples.GetSize(); k++)
    if (samples.At(k) - samples.At(k - 10) > threshold)
      return kTRUE;
  return kFALSE;
}

// An earlier pulse's tail under the baseline: its first and last 100
// samples differing by three baseline RMS.
Bool_t HasSlopingBaseline(const TArrayF &samples, Float_t rms) {
  Double_t early = 0, late = 0;
  for (Int_t k = 0; k < 100; k++) {
    early += samples.At(k) / 100.0;
    late += samples.At(BASELINE_SAMPLES - 100 + k) / 100.0;
  }
  return TMath::Abs(early - late) > 3.0 * TMath::Max(rms, Float_t(1));
}

// One fill pass over a run's tree: baseline RMS, rise time against energy,
// and the trapezoid against the digitiser's energy (in keV on the digitiser's
// scale, the trapezoid scaled by gain).
void FillScalars(const LookRun &run, Double_t gain, TH1D *rms_hist,
                 TH2D *rise_hist, TH2D *trapezoid_hist) {
  TChain chain("features");
  chain.Add(LookFile(run));
  Float_t rms = 0, rise = 0, height = 0;
  Double_t dpp = 0;
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("baseline_rms", 1);
  chain.SetBranchStatus("rise_time", 1);
  chain.SetBranchStatus("trapezoid_height", 1);
  chain.SetBranchStatus("digitizer_energy", 1);
  chain.SetBranchAddress("baseline_rms", &rms);
  chain.SetBranchAddress("rise_time", &rise);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  for (Long64_t i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    Double_t e = dpp / DPP_PER_KEV;
    if (rms_hist)
      rms_hist->Fill(rms);
    if (rise_hist && rise > 0)
      rise_hist->Fill(e, SAMPLE_NS * rise);
    if (trapezoid_hist)
      trapezoid_hist->Fill(e, height / gain / DPP_PER_KEV);
  }
}

// The trapezoid's ADC per digitiser channel: the median ratio above
// 500 channels.
Double_t TrapezoidGain(const LookRun &run) {
  TChain chain("features");
  chain.Add(LookFile(run));
  Float_t height = 0;
  Double_t dpp = 0;
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  std::vector<Double_t> ratios;
  for (Long64_t i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    if (dpp > 500)
      ratios.push_back(height / dpp);
  }
  if (ratios.empty())
    return 1.0;
  return TMath::Median(ratios.size(), ratios.data());
}

// The average baseline-subtracted waveform of the pulses with a digitiser
// energy in [e_lo, e_hi] keV and no pile-up in the record, normalised to its
// maximum.
TGraph *AveragePulse(const LookRun &run, Double_t e_lo, Double_t e_hi,
                     Int_t &n_used) {
  TChain chain("features");
  chain.Add(LookFile(run));
  TArrayF *samples = nullptr;
  Double_t dpp = 0;
  Float_t rms = 0;
  chain.SetBranchAddress("Samples", &samples);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  chain.SetBranchAddress("baseline_rms", &rms);
  std::vector<Double_t> sum;
  n_used = 0;
  for (Long64_t i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    Double_t e = dpp / DPP_PER_KEV;
    if (e < e_lo || e > e_hi)
      continue;
    if (HasEdgeAfterPulse(*samples, rms) || HasSlopingBaseline(*samples, rms))
      continue;
    if (sum.empty())
      sum.assign(samples->GetSize(), 0.0);
    if (samples->GetSize() != (Int_t)sum.size())
      continue;
    for (Int_t k = 0; k < samples->GetSize(); k++)
      sum[k] += samples->At(k);
    n_used++;
  }
  TGraph *g = new TGraph();
  if (n_used == 0)
    return g;
  Double_t peak = 0;
  for (size_t k = 0; k < sum.size(); k++)
    peak = TMath::Max(peak, sum[k]);
  for (size_t k = 0; k < sum.size(); k++)
    g->SetPoint(k, SAMPLE_NS * (Double_t(k) - PRE_SAMPLES), sum[k] / peak);
  return g;
}

// The decay constant of an average pulse: an exponential fitted from 1 us
// after the trigger to the end of the record.
Double_t FitDecay(TGraph *pulse, TF1 *&fit) {
  Double_t t_lo = 1000.0, t_hi = SAMPLE_NS * (POST_SAMPLES - 5);
  fit = new TF1("decay", "[0]*exp(-x/[1])+[2]", t_lo, t_hi);
  fit->SetParameters(1.0, 3000.0, 0.0);
  pulse->Fit(fit, "QRN0");
  return fit->GetParameter(1);
}

void DrawAveragePulses(const std::vector<LookRun> &runs) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TLegend *leg = PlottingUtils::AddLegend(0.45, 0.9, 0.55, 0.9);
  Bool_t first = kTRUE;
  for (size_t i = 0; i < runs.size(); i++) {
    // Am 59.5 keV in the Am run, Pb K-alpha (72.8-75.0 keV) elsewhere.
    Bool_t am = (runs[i].tag == "am");
    Int_t n_used = 0;
    TGraph *g =
        AveragePulse(runs[i], am ? 58.0 : 72.0, am ? 61.0 : 76.0, n_used);
    if (n_used == 0)
      continue;
    TF1 *fit = nullptr;
    Double_t tau_ns = FitDecay(g, fit);
    std::cout << "SI-PULSE " << runs[i].tag << ": " << n_used
              << " pulses, decay " << std::fixed << std::setprecision(1)
              << tau_ns << " ns (" << tau_ns / SAMPLE_NS
              << " samples), offset at the end " << std::setprecision(4)
              << fit->GetParameter(2) << " of the peak" << std::endl;
    PlottingUtils::ConfigureGraph(
        g, runs[i].color, ";Time from trigger [ns];Average pulse (peak = 1)");
    g->SetMarkerSize(0);
    g->SetLineWidth(2);
    g->Draw(first ? "AL" : "L SAME");
    if (first)
      g->GetYaxis()->SetRangeUser(-0.1, 1.15);
    first = kFALSE;
    leg->AddEntry(g,
                  Form("%s (%s), #tau %.0f ns", runs[i].label.Data(),
                       am ? "Am 59.5 keV" : "Pb K#alpha", tau_ns),
                  "l");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "si_average_pulses", "",
                            PlotSaveOptions::kLINEAR);
}

void DrawBaselineNoise(const std::vector<LookRun> &runs) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TLegend *leg = PlottingUtils::AddLegend(0.5, 0.9, 0.62, 0.9);
  Double_t y_max = 0;
  std::vector<TH1D *> hists;
  for (size_t i = 0; i < runs.size(); i++) {
    TH1D *h = new TH1D("rms_" + runs[i].tag,
                       ";Baseline RMS, 350 samples [ADC];Pulses / 0.5 ADC", 160,
                       0, 80);
    h->SetDirectory(nullptr);
    FillScalars(runs[i], 1.0, h, nullptr, nullptr);
    if (h->Integral() > 0)
      h->Scale(1.0 / h->Integral());
    h->GetYaxis()->SetTitle("Fraction of pulses / 0.5 ADC");
    PlottingUtils::ConfigureHistogram(h, runs[i].color);
    y_max = TMath::Max(y_max, h->GetMaximum());
    hists.push_back(h);
    std::cout << "SI-BASELINE " << runs[i].tag << ": rms mean "
              << std::setprecision(2) << h->GetMean() << " ADC, median-ish "
              << h->GetBinCenter(h->GetMaximumBin()) << std::endl;
    leg->AddEntry(
        h, Form("%s: mean %.1f ADC", runs[i].label.Data(), h->GetMean()), "l");
  }
  for (size_t i = 0; i < hists.size(); i++) {
    hists[i]->SetMaximum(3.0 * y_max);
    hists[i]->SetMinimum(1e-5);
    hists[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "si_baseline_rms", "",
                            PlotSaveOptions::kLOG);
}

// The average power spectrum of the pre-trigger samples, from the stored
// waveforms with no pulse before the trigger.
void DrawNoiseSpectrum(const std::vector<LookRun> &runs) {
  Int_t n_fft = 256;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TLegend *leg = PlottingUtils::AddLegend(0.5, 0.9, 0.65, 0.9);
  Bool_t first = kTRUE;
  for (size_t i = 0; i < runs.size(); i++) {
    TChain chain("features");
    chain.Add(LookFile(runs[i]));
    TArrayF *samples = nullptr;
    chain.SetBranchAddress("Samples", &samples);
    std::vector<Double_t> power(n_fft / 2 + 1, 0.0);
    Int_t n_used = 0;
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &n_fft, "R2C ES K");
    for (Long64_t e = 0; e < chain.GetEntries() && n_used < 5000; e++) {
      chain.GetEntry(e);
      Double_t mean = 0;
      for (Int_t k = 0; k < n_fft; k++)
        mean += samples->At(k) / n_fft;
      for (Int_t k = 0; k < n_fft; k++)
        fft->SetPoint(k, samples->At(k) - mean);
      fft->Transform();
      for (Int_t k = 0; k <= n_fft / 2; k++) {
        Double_t re = 0, im = 0;
        fft->GetPointComplex(k, re, im);
        power[k] += (re * re + im * im) / n_fft;
      }
      n_used++;
    }
    delete fft;
    TGraph *g = new TGraph();
    for (Int_t k = 1; k <= n_fft / 2; k++)
      g->SetPoint(k - 1, k * 100.0 / n_fft, power[k] / n_used);
    PlottingUtils::ConfigureGraph(
        g, runs[i].color, ";Frequency [MHz];Baseline power [ADC^{2} per bin]");
    g->SetMarkerSize(0);
    g->SetLineWidth(2);
    g->Draw(first ? "AL" : "L SAME");
    first = kFALSE;
    leg->AddEntry(g, runs[i].label, "l");
    Int_t k_max = 1;
    for (Int_t k = 2; k <= n_fft / 2; k++)
      if (power[k] > power[k_max])
        k_max = k;
    std::cout << "SI-NOISE " << runs[i].tag << ": strongest component at "
              << std::setprecision(2) << k_max * 100.0 / n_fft << " MHz"
              << std::endl;
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "si_baseline_noise_spectrum", "",
                            PlotSaveOptions::kLOG);
}

// Pile-up seen in the stored record: a second rising edge after the pulse,
// or the tail of an earlier pulse sloping the pre-trigger region.
void CountPileup(const std::vector<LookRun> &runs) {
  for (size_t i = 0; i < runs.size(); i++) {
    TChain chain("features");
    chain.Add(LookFile(runs[i]));
    TArrayF *samples = nullptr;
    Float_t rms = 0;
    chain.SetBranchAddress("Samples", &samples);
    chain.SetBranchAddress("baseline_rms", &rms);
    Long64_t n = chain.GetEntries(), n_after = 0, n_before = 0;
    for (Long64_t e = 0; e < n; e++) {
      chain.GetEntry(e);
      if (HasEdgeAfterPulse(*samples, rms))
        n_after++;
      if (HasSlopingBaseline(*samples, rms))
        n_before++;
    }
    std::cout << "SI-PILEUP " << runs[i].tag << ": " << n << " pulses, "
              << std::setprecision(1) << 100.0 * n_after / n
              << "% with a second edge after the pulse, "
              << 100.0 * n_before / n << "% with a sloping pre-trigger baseline"
              << std::endl;
  }
}

// Records whose baseline sits near the bottom of the ADC range: their noise
// is clipped, so their baseline RMS is a few ADC instead of ~33.
const Float_t FLOOR_BASELINE_RMS = 8.0;

void DrawFloorBaselineExamples(const LookRun &run) {
  TChain chain("features");
  chain.Add(LookFile(run));
  TArrayF *samples = nullptr;
  Float_t rms = 0;
  Double_t dpp = 0;
  chain.SetBranchAddress("Samples", &samples);
  chain.SetBranchAddress("baseline_rms", &rms);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TLegend *leg = PlottingUtils::AddLegend(0.45, 0.9, 0.62, 0.9);
  const Int_t colors[4] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7};
  Int_t shown = 0;
  Long64_t n_floor = 0;
  for (Long64_t i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    if (rms >= FLOOR_BASELINE_RMS)
      continue;
    n_floor++;
    if (shown >= 4)
      continue;
    TGraph *g = new TGraph();
    for (Int_t k = 0; k < samples->GetSize(); k++)
      g->SetPoint(k, SAMPLE_NS * (k - PRE_SAMPLES), samples->At(k));
    PlottingUtils::ConfigureGraph(
        g, colors[shown], ";Time from trigger [ns];Above baseline [ADC]");
    g->SetMarkerSize(0);
    g->Draw(shown == 0 ? "AL" : "L SAME");
    leg->AddEntry(g, Form("Baseline RMS %.1f ADC, DPP channel %.0f", rms, dpp),
                  "l");
    shown++;
  }
  std::cout << "SI-FLOOR " << run.tag << ": " << n_floor << " of "
            << chain.GetEntries() << " records with baseline RMS below "
            << FLOOR_BASELINE_RMS << " ADC" << std::endl;
  if (shown == 0)
    return;
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "si_floor_baseline_examples_" + run.tag, "",
                            PlotSaveOptions::kLINEAR);
}

// The Am-241 59.54 keV line's width in the digitiser's energy and in the
// trapezoid (on the digitiser's scale), records with a normal baseline only.
void CompareResolution(const LookRun &run) {
  Double_t gain = TrapezoidGain(run);
  TChain chain("features");
  chain.Add(LookFile(run));
  Float_t rms = 0, height = 0;
  Double_t dpp = 0;
  chain.SetBranchAddress("baseline_rms", &rms);
  chain.SetBranchAddress("trapezoid_height", &height);
  chain.SetBranchAddress("digitizer_energy", &dpp);
  TH1D *h_dpp = new TH1D("res_dpp_" + run.tag, ";Energy [keV];Counts / 50 eV",
                         160, 55, 63);
  TH1D *h_trap = new TH1D("res_trap_" + run.tag, ";Energy [keV];Counts / 50 eV",
                          160, 55, 63);
  h_dpp->SetDirectory(nullptr);
  h_trap->SetDirectory(nullptr);
  for (Long64_t i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    if (rms < FLOOR_BASELINE_RMS)
      continue;
    h_dpp->Fill(dpp / DPP_PER_KEV);
    h_trap->Fill(height / gain / DPP_PER_KEV);
  }
  TH1D *hists[2] = {h_dpp, h_trap};
  const char *names[2] = {"digitiser DPP", "trapezoid"};
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  TLegend *leg = PlottingUtils::AddLegend(0.18, 0.62, 0.7, 0.9);
  for (Int_t k = 0; k < 2; k++) {
    TH1D *h = hists[k];
    Double_t top = h->GetBinCenter(h->GetMaximumBin());
    TF1 fit(Form("res_fit_%d", k), "gaus", top - 0.6, top + 0.6);
    fit.SetParameters(h->GetMaximum(), top, 0.3);
    h->Fit(&fit, "QRN0L");
    PlottingUtils::ConfigureHistogram(h, k == 0 ? kGray + 2 : kBlue + 1);
    h->Draw(k == 0 ? "HIST" : "HIST SAME");
    leg->AddEntry(h,
                  Form("%s: #sigma %.0f eV at %.2f keV", names[k],
                       1000 * fit.GetParameter(2), fit.GetParameter(1)),
                  "l");
    std::cout << "SI-RESOLUTION " << run.tag << " " << names[k] << ": sigma "
              << std::setprecision(1) << 1000 * fit.GetParameter(2) << " +/- "
              << 1000 * fit.GetParError(2) << " eV" << std::endl;
  }
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "si_am_resolution_" + run.tag, "",
                            PlotSaveOptions::kLINEAR);
}

void DrawRiseVsHeight(const LookRun &run) {
  TH2D *h = new TH2D("rise_" + run.tag, "", 240, 0, 120, 100, 0, 1000);
  h->SetDirectory(nullptr);
  FillScalars(run, 1.0, nullptr, h, nullptr);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      ";Digitiser energy [keV];Rise time 10-90% [ns];Pulses / 0.5 keV / 10 ns");
  PlottingUtils::SaveFigure(canvas, "si_rise_vs_energy_" + run.tag, "",
                            PlotSaveOptions::kLINEAR);
}

void DrawTrapezoidVsDpp(const LookRun &run) {
  Double_t gain = TrapezoidGain(run);
  TH2D *h = new TH2D("trap_" + run.tag, "", 240, 0, 120, 240, 0, 120);
  h->SetDirectory(nullptr);
  FillScalars(run, gain, nullptr, nullptr, h);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDraw2DHistogram(
      h, canvas,
      ";Digitiser energy [keV];Trapezoid on the same scale [keV];Pulses / "
      "0.5 keV / 0.5 keV");
  PlottingUtils::SaveFigure(canvas, "si_trapezoid_vs_dpp_" + run.tag, "",
                            PlotSaveOptions::kLINEAR);
  std::cout << "SI-GAIN " << run.tag << ": trapezoid ADC per DPP channel "
            << std::setprecision(4) << gain << std::endl;
}

} // namespace SiLook

void SiWaveformLook(Int_t stage, Double_t decay_samples = 0) {
  using namespace SiLook;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(
      PlotSaveFormat::kPNG, project_root + "/plots/studies/psi_si/waveforms",
      LOOK_DIR);
  if (stage == 1) {
    Process(decay_samples);
    return;
  }
  std::vector<LookRun> runs = Runs();
  DrawAveragePulses(runs);
  DrawBaselineNoise(runs);
  DrawNoiseSpectrum(runs);
  CountPileup(runs);
  for (size_t i = 0; i < runs.size(); i++) {
    DrawFloorBaselineExamples(runs[i]);
    DrawRiseVsHeight(runs[i]);
    DrawTrapezoidVsDpp(runs[i]);
  }
  CompareResolution(runs[0]);
}
