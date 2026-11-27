#ifndef INSITU_FIT_STATE_HPP
#define INSITU_FIT_STATE_HPP

// An in-situ fit's saved state (a .simroofits file) and the channel seeds
// built from it, for the macros that rerun production's in-situ fit from
// its state: InjectionTest.cpp (the testbed) and FitBenchmark.cpp.

#include "GainTransferFit.hpp"

#include <TMath.h>
#include <TString.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct ParamValue {
  Double_t value = 0, error = 0;
};

// An in-situ fit's window and parameters by name ("bkg:Mu1", "sig:Sigma3",
// ...), as a .simroofits file holds them.
struct FitState {
  Bool_t ok = kFALSE;
  Double_t range_lo = 0, range_hi = 0;
  std::map<std::string, ParamValue> params;
  Bool_t Has(const std::string &name) const { return params.count(name) > 0; }
  Double_t Get(const std::string &name, Double_t fallback = -1e30) const {
    std::map<std::string, ParamValue>::const_iterator it = params.find(name);
    return it == params.end() ? fallback : it->second.value;
  }
  Double_t Error(const std::string &name) const {
    std::map<std::string, ParamValue>::const_iterator it = params.find(name);
    return it == params.end() ? 0.0 : it->second.error;
  }
  void Set(const std::string &name, Double_t value, Double_t error) {
    params[name].value = value;
    params[name].error = error;
  }
};

inline TString Indexed(const TString &channel, const char *param, Int_t index) {
  return channel + ":" + param + TString::Format("%d", index);
}

inline Bool_t HasGeLine(const FitState &s) {
  return s.Has("sig:Mu3") && s.Has("sig:Sigma3") && s.Has("sig:GausAmplitude3");
}

inline FitState ReadFitState(const TString &path) {
  FitState s;
  std::ifstream in(path.Data());
  if (!in.is_open())
    return s;
  std::string token;
  if (!(in >> token) || token != "RANGE")
    return s;
  if (!(in >> s.range_lo >> s.range_hi))
    return s;
  std::string name;
  Double_t value, error;
  Int_t fixed;
  while (in >> name >> value >> error >> fixed)
    s.Set(name, value, error);
  s.ok = HasGeLine(s);
  return s;
}

// The in-situ model's own free parameters (GainTransferFit.hpp's
// FitBkgGeSimFromSeed): every Pb width and tail is bkg:Ka1's, so the signal
// channel adds only its centroids, yields, continuum and the Ge width.
inline std::vector<std::string> FreeStateParams() {
  const char *names[] = {"bkg:Mu1",
                         "bkg:Sigma1",
                         "bkg:GausAmplitude1",
                         "bkg:LowExpTailAmplitude1",
                         "bkg:LowExpTailRatio1",
                         "bkg:Mu2",
                         "bkg:GausAmplitude2",
                         "bkg:BkgConstant",
                         "sig:Mu1",
                         "sig:GausAmplitude1",
                         "sig:Mu2",
                         "sig:GausAmplitude2",
                         "sig:Mu3",
                         "sig:Sigma3",
                         "sig:GausAmplitude3",
                         "sig:BkgConstant"};
  std::vector<std::string> out(names, names + sizeof(names) / sizeof(names[0]));
  return out;
}

// Written in Analysis-Utilities' current convention, whose marker line says
// the exponential-tail amplitudes are full-range ratios.
inline Bool_t WriteFitState(const FitState &s, const TString &path) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    std::cerr << "ERROR: cannot write fit state " << path << std::endl;
    return kFALSE;
  }
  out << std::setprecision(17);
  out << "RANGE " << s.range_lo << " " << s.range_hi << std::endl;
  out << "ExpTailAmplitudeConvention 1 0 1" << std::endl;
  std::vector<std::string> names = FreeStateParams();
  for (size_t i = 0; i < names.size(); i++) {
    if (!s.Has(names[i]))
      continue;
    out << names[i] << " " << s.Get(names[i]) << " " << s.Error(names[i])
        << " 0" << std::endl;
  }
  return kTRUE;
}

inline FitState FitStateOf(const BkgGeSimResult &fit, Double_t range_lo,
                           Double_t range_hi) {
  FitState s;
  s.range_lo = range_lo;
  s.range_hi = range_hi;
  const FitResult *channels[2] = {&fit.bkg_channel, &fit.sig_channel};
  for (Int_t c = 0; c < 2; c++) {
    const std::vector<FitParameterDiagnostic> &diags =
        channels[c]->parameter_diagnostics;
    for (size_t i = 0; i < diags.size(); i++)
      s.Set(diags[i].name, diags[i].value, diags[i].error);
  }
  s.ok = fit.valid && HasGeLine(s);
  return s;
}

// Multiplies a yield by f, its error by sqrt(f), as for a run f times as long.
inline void ScaleYield(FitState &s, const std::string &name, Double_t f) {
  if (s.Has(name))
    s.Set(name, f * s.Get(name), TMath::Sqrt(f) * s.Error(name));
}

inline Double_t PbSigma(const FitState &s, Int_t index) {
  return s.Get(Indexed("bkg", "Sigma", index).Data(), s.Get("bkg:Sigma1"));
}

// A channel's seed for RunBkgGeSim, from the state its fit starts from. The
// loaded state overrides most of it; what it sets is the tail the 95.902 keV
// line is measured with for the escapes (the background channel's last peak),
// and anything the state lacks (an old state has no signal-run Pb centroids,
// which then start from the background run's).
inline FitResult ChannelSeed(const FitState &s, const TString &channel) {
  const Int_t n_peaks = (channel == "sig") ? 3 : 2;
  const Double_t tail_r = s.Get("bkg:LowExpTailAmplitude1", 0.0);
  const Double_t tau = s.Get("bkg:LowExpTailRatio1", 1.0);
  FitResult seed;
  seed.peaks.resize(n_peaks);
  for (Int_t i = 1; i <= n_peaks; i++) {
    PeakFitResult &peak = seed.peaks[i - 1];
    const Bool_t is_ge = (i == 3);
    peak.mu = s.Get(Indexed(channel, "Mu", i).Data(),
                    s.Get(Indexed("bkg", "Mu", i).Data()));
    peak.sigma = is_ge ? s.Get("sig:Sigma3") : PbSigma(s, i);
    peak.gaus_amplitude = s.Get(Indexed(channel, "GausAmplitude", i).Data());
    peak.low_exp_tail_amplitude = tail_r * peak.gaus_amplitude;
    peak.low_exp_tail_ratio = tau;
  }
  seed.bkg_constant = s.Get((channel + ":BkgConstant").Data(), 0.0);
  seed.lin_bkg_slope = s.Get((channel + ":BkgSlope").Data(), 0.0);
  seed.valid = kTRUE;
  return seed;
}

#endif // INSITU_FIT_STATE_HPP
