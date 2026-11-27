// Fits of the germanium ionization yield with the Lindhard model and with the
// excess-ionization ("modified Lindhard") model of Kavner et al.,
//   E_ee = Q_L(kappa, E_nr) E_nr + eta* W,  W = 2.96 eV,
// to the measurements in GeYieldData.hpp, with the 255 eV capture recoils
// measured against this work's 68.75 keV gamma ray.
//
// Each source's stated common uncertainties are nuisance parameters, one per
// common effect, each a standard normal pull shared by the points it moves:
// a shift of the reported ionization (a subtracted level or gamma energy, an
// energy-scale offset) or of the true recoil energy (a neutron-energy or
// recoil-energy scale). Where a source puts that part inside its printed
// errors, it is taken out of each point's error first. This is the fit form
// of a BLUE: with linear effects it equals a chi2 with the full covariance.
// The 68.75 keV gamma energy, which every capture-line ionization subtracts,
// takes that covariance form itself: it stays at its measured value and its
// uncertainty is a fully correlated term between those points, so a model's
// tension with them shows on the points and not as a moved gamma energy.
// A point-to-point recoil-energy error enters through the model's slope.
// Asymmetric errors use the side toward the model.
//
// Fits: all data and the HPGe data alone, each with the Lindhard model
// (kappa) and the excess-ionization model (kappa, eta*), correlated and, for
// comparison with Kavner et al., with every error as uncorrelated. All data
// includes Jones and Kraner's 1975 capture point, which Kavner et al. leave
// out of their fit without giving a reason.
//
// Reads results/recommended_value.result and results/recoil_energy.result;
// writes results/modified_lindhard.result.
//
// Usage: root -b -q 'ModifiedLindhard.cpp+'
#include "FigureStyle.hpp"
#include "GeIonizationModel.hpp"
#include "GeYieldData.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "ProjectPaths.hpp"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

using GeIonizationModel::IonizationKeV;

// The yield figures' plot width, wider than the standard 1200 px so the
// frame has room beside the long legend.
const Int_t YIELD_PLOT_WIDTH_PX = 1600;

enum Detector { kGeLi, kHPGe };

// How one nuisance moves one point, per unit pull: the true recoil energy
// and the reported ionization, both keV.
struct Coupling {
  Int_t nuisance;
  Double_t recoil_kev;
  Double_t ionization_kev;
};

// A measurement as a yield (reported ionization over reported recoil
// energy, a fraction) with its uncorrelated errors and its couplings.
struct YieldPoint {
  TString dataset;
  Detector detector;
  Double_t recoil_kev;
  Double_t yield;
  Double_t error_up;
  Double_t error_down;
  Double_t recoil_error_kev;
  std::vector<Coupling> couplings;
};

struct YieldData {
  std::vector<YieldPoint> points;
  std::vector<TString> nuisances;
  // The nuisance that enters as a covariance instead of a fitted pull.
  Int_t covariance_nuisance = -1;
};

// The 68.75 keV gamma ray the capture-line ionizations are measured against,
// the shift that puts Jones and Kraner's calibration on today's reference
// energies, and the mean capture recoil.
struct GammaReference {
  Double_t energy_kev;
  Double_t error_kev;
  Double_t jk_shift_kev;
  Double_t recoil_kev;
  Int_t nuisance;
};

Int_t AddNuisance(YieldData &data, const TString &name) {
  data.nuisances.push_back(name);
  return (Int_t)data.nuisances.size() - 1;
}

YieldPoint MakePoint(const TString &dataset, Detector detector,
                     Double_t recoil_kev, Double_t yield, Double_t error_up,
                     Double_t error_down) {
  YieldPoint point;
  point.dataset = dataset;
  point.detector = detector;
  point.recoil_kev = recoil_kev;
  point.yield = yield;
  point.error_up = error_up;
  point.error_down = error_down;
  point.recoil_error_kev = 0;
  return point;
}

Coupling MakeCoupling(Int_t nuisance, Double_t recoil_kev,
                      Double_t ionization_kev) {
  Coupling coupling;
  coupling.nuisance = nuisance;
  coupling.recoil_kev = recoil_kev;
  coupling.ionization_kev = ionization_kev;
  return coupling;
}

// An error with a common part removed, in quadrature.
Double_t WithoutCommonPart(Double_t error, Double_t common) {
  Double_t remainder = error * error - common * common;
  if (remainder <= 0) {
    std::cerr << "ERROR: a common part (" << common
              << ") exceeds the error it is inside (" << error << ")"
              << std::endl;
    return 0;
  }
  return TMath::Sqrt(remainder);
}

Double_t IonizationSlope(Double_t recoil_kev, Double_t kappa, Double_t eta) {
  Double_t step = 1e-4 * recoil_kev;
  return (IonizationKeV(recoil_kev + step, kappa, eta) -
          IonizationKeV(recoil_kev - step, kappa, eta)) /
         (2.0 * step);
}

void AddYieldTable(YieldData &data, const TString &dataset, Detector detector,
                   const GeYieldData::SourcePoint *table, Int_t n) {
  for (Int_t i = 0; i < n; i++)
    data.points.push_back(MakePoint(dataset, detector, table[i].recoil_kev,
                                    table[i].value, table[i].error_up,
                                    table[i].error_down));
}

void AddChasman1967(YieldData &data) {
  using namespace GeYieldData;
  Int_t level_nuisance[2] = {
      AddNuisance(data, "Chasman 1967 Ge72 690 keV level"),
      AddNuisance(data, "Chasman 1967 Ge74 596 keV level")};
  Int_t neutron_nuisance = AddNuisance(data, "Chasman 1967 neutron energy");
  Double_t neutron_fraction =
      CHASMAN1967_NEUTRON_RECOIL_ERR_KEV / CHASMAN1967_NEUTRON_RECOIL_KEV;
  for (Int_t i = 0; i < CHASMAN1967_N; i++) {
    const SourcePoint &source = CHASMAN1967[i];
    Double_t level_part = CHASMAN1967_LEVEL_ERR_KEV / source.recoil_kev;
    Double_t neutron_part = source.value * neutron_fraction;
    Double_t common =
        TMath::Sqrt(level_part * level_part + neutron_part * neutron_part);
    Double_t error = WithoutCommonPart(source.error_up, common);
    YieldPoint point = MakePoint("Chasman 1967", kGeLi, source.recoil_kev,
                                 source.value, error, error);
    point.couplings.push_back(MakeCoupling(level_nuisance[CHASMAN1967_LEVEL[i]],
                                           0, CHASMAN1967_LEVEL_ERR_KEV));
    point.couplings.push_back(MakeCoupling(
        neutron_nuisance, source.recoil_kev * neutron_fraction, 0));
    data.points.push_back(point);
  }
}

void AddChasman1968(YieldData &data) {
  using namespace GeYieldData;
  Int_t level_nuisance = AddNuisance(data, "Chasman 1968 Ge72 691.4 keV level");
  for (Int_t i = 0; i < CHASMAN1968_N; i++) {
    const SourcePoint &source = CHASMAN1968[i];
    Double_t error_kev =
        WithoutCommonPart(source.error_up, CHASMAN1968_LEVEL_ERR_KEV);
    YieldPoint point =
        MakePoint("Chasman 1968", kGeLi, source.recoil_kev,
                  source.value / source.recoil_kev,
                  error_kev / source.recoil_kev, error_kev / source.recoil_kev);
    point.couplings.push_back(
        MakeCoupling(level_nuisance, 0, CHASMAN1968_LEVEL_ERR_KEV));
    data.points.push_back(point);
  }
}

// Jones and Kraner's recoil-energy scale error at a recoil energy: the
// straight line through their two stated values.
Double_t Jones1971ScaleErrKeV(Double_t recoil_kev) {
  using namespace GeYieldData;
  Double_t slope =
      (JONES1971_SCALE_ERR_KEV[1] - JONES1971_SCALE_ERR_KEV[0]) /
      (JONES1971_SCALE_ANCHOR_KEV[1] - JONES1971_SCALE_ANCHOR_KEV[0]);
  return JONES1971_SCALE_ERR_KEV[0] +
         slope * (recoil_kev - JONES1971_SCALE_ANCHOR_KEV[0]);
}

// The ionizations put against this work's gamma energy: the measured line is
// their ionization plus the level energy they subtracted, moved onto today's
// calibration like their gamma-only values in RecommendedValue.cpp.
void AddJones1971(YieldData &data, const GammaReference &gamma) {
  using namespace GeYieldData;
  Int_t scale_nuisance = AddNuisance(data, "Jones 1971 recoil-energy scale");
  for (Int_t i = 0; i < JONES1971_N; i++) {
    const SourcePoint &source = JONES1971[i];
    Double_t ionization_kev = source.value + JONES1971_LEVEL_KEV +
                              gamma.jk_shift_kev - gamma.energy_kev;
    YieldPoint point = MakePoint("Jones 1971", kGeLi, source.recoil_kev,
                                 ionization_kev / source.recoil_kev,
                                 source.error_up / source.recoil_kev,
                                 source.error_down / source.recoil_kev);
    point.recoil_error_kev = JONES1971_RELATIVE_RECOIL_ERR_KEV;
    point.couplings.push_back(MakeCoupling(
        scale_nuisance, Jones1971ScaleErrKeV(source.recoil_kev), 0));
    point.couplings.push_back(MakeCoupling(gamma.nuisance, 0, gamma.error_kev));
    data.points.push_back(point);
  }
}

// A capture-line ionization at the mean recoil: the gamma + recoil line minus
// this work's gamma energy, the line first moved onto today's calibration if
// it was measured on Jones and Kraner's.
void AddCaptureLine(YieldData &data, const TString &dataset, Detector detector,
                    Double_t gamma_plus_recoil_kev, Double_t error_kev,
                    Bool_t on_jones_kraner_scale, const GammaReference &gamma) {
  Double_t line_kev = gamma_plus_recoil_kev;
  if (on_jones_kraner_scale)
    line_kev += gamma.jk_shift_kev;
  Double_t ionization_kev = line_kev - gamma.energy_kev;
  YieldPoint point = MakePoint(
      dataset, detector, gamma.recoil_kev, ionization_kev / gamma.recoil_kev,
      error_kev / gamma.recoil_kev, error_kev / gamma.recoil_kev);
  point.couplings.push_back(MakeCoupling(gamma.nuisance, 0, gamma.error_kev));
  data.points.push_back(point);
}

Double_t BarError(Double_t value, Double_t bar_end, Double_t fallback) {
  if (bar_end == GeYieldData::NO_BAR)
    return fallback;
  return TMath::Abs(bar_end - value);
}

void AddBaudis1998(YieldData &data) {
  using namespace GeYieldData;
  for (Int_t i = 0; i < BAUDIS1998_N; i++) {
    const FigurePoint &source = BAUDIS1998[i];
    data.points.push_back(MakePoint(
        "Baudis 1998", kHPGe, source.recoil_kev, source.value,
        BarError(source.value, source.bar_high, BAUDIS1998_MARKER_RADIUS),
        BarError(source.value, source.bar_low, BAUDIS1998_MARKER_RADIUS)));
  }
}

void AddCoGeNT(YieldData &data) {
  using namespace GeYieldData;
  for (Int_t i = 0; i < COGENT_N; i++) {
    const FigurePoint &source = COGENT[i];
    Double_t error_down = BarError(source.value, source.bar_low, 0);
    Double_t error_up = BarError(source.value, source.bar_high, error_down);
    data.points.push_back(MakePoint("CoGeNT", kHPGe, source.recoil_kev,
                                    source.value / 100.0, error_up / 100.0,
                                    error_down / 100.0));
  }
}

void AddCollar2021IronFilter(YieldData &data) {
  using namespace GeYieldData;
  Int_t intercept_nuisance =
      AddNuisance(data, "Collar 2021 calibration intercept");
  for (Int_t i = 0; i < COLLAR2021_FEF_N; i++) {
    const SourcePoint &source = COLLAR2021_FEF[i];
    Double_t intercept_percent =
        100.0 * COLLAR2021_INTERCEPT_ERR_KEV / source.recoil_kev;
    YieldPoint point = MakePoint(
        "Collar 2021 FeF", kHPGe, source.recoil_kev, source.value / 100.0,
        WithoutCommonPart(source.error_up, intercept_percent) / 100.0,
        WithoutCommonPart(source.error_down, intercept_percent) / 100.0);
    point.couplings.push_back(
        MakeCoupling(intercept_nuisance, 0, COLLAR2021_INTERCEPT_ERR_KEV));
    data.points.push_back(point);
  }
}

// The energy scale is an ionization offset of 12 eVee and each beam set's
// neutron energy a scale on its recoil energies (E_nr is proportional to
// E_n at a fixed angle). The band's half-width is the energy scale alone
// above about 1.5 keVnr; what it holds beyond that at low energy is taken as
// the trigger efficiency, which the paper calls dominant below 0.8 keVnr and
// negligible above 1.2 keVnr. The beam terms differ between sets at one
// energy, so a single band cannot hold them; they are taken from Table 1.
void AddBonhomme2022(YieldData &data) {
  using namespace GeYieldData;
  Int_t scale_nuisance = AddNuisance(data, "Bonhomme 2022 energy scale");
  Int_t trigger_nuisance =
      AddNuisance(data, "Bonhomme 2022 trigger efficiency");
  Int_t beam_nuisance[BONHOMME2022_N_BEAMS];
  for (Int_t b = 0; b < BONHOMME2022_N_BEAMS; b++)
    beam_nuisance[b] = AddNuisance(
        data, Form("Bonhomme 2022 beam %.1f keV", BONHOMME2022_BEAM_KEV[b]));
  for (Int_t i = 0; i < BONHOMME2022_N; i++) {
    const BonhommePoint &source = BONHOMME2022[i];
    YieldPoint point =
        MakePoint("Bonhomme 2022", kHPGe, source.recoil_kev, source.q,
                  source.bar_high - source.q, source.q - source.bar_low);
    point.couplings.push_back(
        MakeCoupling(scale_nuisance, 0, BONHOMME2022_ENERGY_SCALE_ERR_KEV));
    Double_t scale_part = BONHOMME2022_ENERGY_SCALE_ERR_KEV / source.recoil_kev;
    Double_t trigger_variance =
        source.band_half_width * source.band_half_width -
        scale_part * scale_part;
    if (trigger_variance > 0)
      point.couplings.push_back(
          MakeCoupling(trigger_nuisance, 0,
                       TMath::Sqrt(trigger_variance) * source.recoil_kev));
    Int_t beam = source.beam_set;
    point.couplings.push_back(
        MakeCoupling(beam_nuisance[beam],
                     source.recoil_kev * BONHOMME2022_BEAM_ERR_KEV[beam] /
                         BONHOMME2022_BEAM_KEV[beam],
                     0));
    data.points.push_back(point);
  }
}

// chi2 of the selected points; parameters: kappa, eta*, then one pull per
// nuisance (held at zero when not fitted).
// A point's residual (measured minus predicted yield), the variance of its
// own error, and its yield shift per unit of the covariance nuisance.
struct Residual {
  Double_t value;
  Double_t variance;
  Double_t shared;
};

// chi2 of the selected points; parameters: kappa, eta*, then one pull per
// nuisance (held at zero when not fitted). In the correlated form the
// covariance nuisance is never fitted: the points it acts on share it
// through their covariance V = D + s s^T (D their own variances, s their
// shifts), whose inverse is D^-1 - D^-1 s s^T D^-1 / (1 + s^T D^-1 s).
class YieldChi2 {
public:
  YieldChi2(const std::vector<const YieldPoint *> &points, Int_t n_nuisances,
            Int_t covariance_nuisance, Bool_t correlated)
      : points_(points), n_nuisances_(n_nuisances),
        covariance_nuisance_(covariance_nuisance), correlated_(correlated) {}

  Residual PointResidual(const YieldPoint &point,
                         const Double_t *parameters) const {
    Double_t kappa = parameters[0];
    Double_t eta = parameters[1];
    Double_t true_recoil_kev = point.recoil_kev;
    Double_t ionization_shift_kev = 0;
    if (correlated_) {
      for (size_t c = 0; c < point.couplings.size(); c++) {
        const Coupling &coupling = point.couplings[c];
        if (coupling.nuisance == covariance_nuisance_)
          continue;
        Double_t pull = parameters[2 + coupling.nuisance];
        true_recoil_kev += coupling.recoil_kev * pull;
        ionization_shift_kev += coupling.ionization_kev * pull;
      }
    }
    Double_t slope = IonizationSlope(true_recoil_kev, kappa, eta);
    Residual residual;
    residual.shared = 0;
    Double_t common_variance = 0;
    for (size_t c = 0; c < point.couplings.size(); c++) {
      const Coupling &coupling = point.couplings[c];
      Double_t effect =
          (slope * coupling.recoil_kev + coupling.ionization_kev) /
          point.recoil_kev;
      if (!correlated_)
        common_variance += effect * effect;
      else if (coupling.nuisance == covariance_nuisance_)
        residual.shared += effect;
    }
    Double_t predicted =
        (IonizationKeV(true_recoil_kev, kappa, eta) + ionization_shift_kev) /
        point.recoil_kev;
    Double_t yield_error =
        (predicted > point.yield) ? point.error_up : point.error_down;
    Double_t recoil_term = slope * point.recoil_error_kev / point.recoil_kev;
    residual.value = point.yield - predicted;
    residual.variance =
        yield_error * yield_error + recoil_term * recoil_term + common_variance;
    return residual;
  }

  // Measured minus predicted over the point's whole error: its own, and its
  // share of the covariance nuisance.
  Double_t PointPull(const YieldPoint &point,
                     const Double_t *parameters) const {
    Residual residual = PointResidual(point, parameters);
    return residual.value /
           TMath::Sqrt(residual.variance + residual.shared * residual.shared);
  }

  // chi2 of a set of points, those sharing the covariance nuisance through
  // their inverse covariance.
  Double_t SetChi2(const std::vector<const YieldPoint *> &points,
                   const Double_t *parameters) const {
    Double_t chi2 = 0, sum_rs = 0, sum_ss = 0;
    for (size_t i = 0; i < points.size(); i++) {
      Residual residual = PointResidual(*points[i], parameters);
      chi2 += residual.value * residual.value / residual.variance;
      sum_rs += residual.value * residual.shared / residual.variance;
      sum_ss += residual.shared * residual.shared / residual.variance;
    }
    return chi2 - sum_rs * sum_rs / (1.0 + sum_ss);
  }

  Double_t Evaluate(const Double_t *parameters) const {
    Double_t chi2 = SetChi2(points_, parameters);
    if (correlated_) {
      for (Int_t k = 0; k < n_nuisances_; k++)
        if (k != covariance_nuisance_)
          chi2 += parameters[2 + k] * parameters[2 + k];
    }
    return chi2;
  }

private:
  std::vector<const YieldPoint *> points_;
  Int_t n_nuisances_;
  Int_t covariance_nuisance_;
  Bool_t correlated_;
};

struct FitConfig {
  TString key;
  TString title;
  std::vector<TString> datasets;
  Bool_t hpge_only;
  Bool_t correlated;
};

struct FitOutcome {
  Bool_t valid;
  Int_t n_points;
  Int_t n_fitted_nuisances;
  Double_t kappa;
  Double_t kappa_err;
  Double_t eta;
  Double_t eta_err;
  Double_t chi2;
  Int_t ndf;
  std::vector<Double_t> best;
  std::vector<Double_t> pulls;
  std::vector<Bool_t> fitted;
  std::vector<TString> dataset_names;
  std::vector<Double_t> dataset_chi2;
  std::vector<Int_t> dataset_points;
};

Bool_t Selected(const YieldPoint &point, const FitConfig &config) {
  if (config.hpge_only && point.detector != kHPGe)
    return kFALSE;
  for (size_t d = 0; d < config.datasets.size(); d++)
    if (point.dataset == config.datasets[d])
      return kTRUE;
  return kFALSE;
}

FitOutcome RunFit(const YieldData &data, const FitConfig &config,
                  Bool_t excess_ionization) {
  FitOutcome outcome;
  outcome.valid = kFALSE;
  std::vector<const YieldPoint *> points;
  Int_t n_nuisances = (Int_t)data.nuisances.size();
  outcome.fitted.assign(n_nuisances, kFALSE);
  for (size_t i = 0; i < data.points.size(); i++) {
    if (!Selected(data.points[i], config))
      continue;
    points.push_back(&data.points[i]);
    for (size_t c = 0; c < data.points[i].couplings.size(); c++) {
      Int_t nuisance = data.points[i].couplings[c].nuisance;
      outcome.fitted[nuisance] =
          config.correlated && nuisance != data.covariance_nuisance;
    }
  }
  outcome.n_points = (Int_t)points.size();

  YieldChi2 chi2(points, n_nuisances, data.covariance_nuisance,
                 config.correlated);
  ROOT::Math::Functor functor(&chi2, &YieldChi2::Evaluate, 2 + n_nuisances);
  ROOT::Math::Minimizer *minimizer =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetFunction(functor);
  minimizer->SetStrategy(2);
  minimizer->SetTolerance(1e-4);
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetPrintLevel(0);
  minimizer->SetLimitedVariable(0, "kappa", 0.157, 0.001, 0.02, 1.0);
  if (excess_ionization)
    minimizer->SetVariable(1, "eta", 5.0, 0.1);
  else
    minimizer->SetFixedVariable(1, "eta", 0.0);
  outcome.n_fitted_nuisances = 0;
  for (Int_t k = 0; k < n_nuisances; k++) {
    std::string name = Form("pull_%d", k);
    if (outcome.fitted[k]) {
      minimizer->SetVariable(2 + k, name, 0.0, 0.1);
      outcome.n_fitted_nuisances++;
    } else {
      minimizer->SetFixedVariable(2 + k, name, 0.0);
    }
  }
  Bool_t converged = minimizer->Minimize();
  minimizer->Hesse();
  outcome.valid = converged && minimizer->Status() == 0;
  const Double_t *best = minimizer->X();
  const Double_t *errors = minimizer->Errors();
  outcome.kappa = best[0];
  outcome.kappa_err = errors[0];
  outcome.eta = best[1];
  outcome.eta_err = excess_ionization ? errors[1] : 0;
  outcome.chi2 = minimizer->MinValue();
  outcome.ndf = outcome.n_points - (excess_ionization ? 2 : 1);
  outcome.best.assign(best, best + 2 + n_nuisances);
  for (Int_t k = 0; k < n_nuisances; k++)
    outcome.pulls.push_back(best[2 + k]);
  // Each dataset's chi2 from its own points and their own covariance block.
  std::vector<std::vector<const YieldPoint *>> by_dataset;
  for (size_t i = 0; i < points.size(); i++) {
    Int_t index = -1;
    for (size_t d = 0; d < outcome.dataset_names.size(); d++)
      if (outcome.dataset_names[d] == points[i]->dataset)
        index = (Int_t)d;
    if (index < 0) {
      outcome.dataset_names.push_back(points[i]->dataset);
      by_dataset.push_back(std::vector<const YieldPoint *>());
      index = (Int_t)outcome.dataset_names.size() - 1;
    }
    by_dataset[index].push_back(points[i]);
  }
  for (size_t d = 0; d < by_dataset.size(); d++) {
    outcome.dataset_chi2.push_back(chi2.SetChi2(by_dataset[d], best));
    outcome.dataset_points.push_back((Int_t)by_dataset[d].size());
  }
  delete minimizer;
  return outcome;
}

void PrintFit(const YieldData &data, const FitConfig &config,
              const FitOutcome &lindhard, const FitOutcome &excess) {
  std::cout << std::endl;
  std::cout << "FIT " << config.key << ": " << config.title << " ("
            << excess.n_points << " points, " << excess.n_fitted_nuisances
            << " nuisance parameters)" << std::endl;
  std::cout << std::fixed;
  std::cout << "  Lindhard:           kappa " << std::setprecision(4)
            << lindhard.kappa << " +/- " << lindhard.kappa_err
            << "                    chi2 " << std::setprecision(1)
            << lindhard.chi2 << " / " << lindhard.ndf << "  p "
            << std::scientific << std::setprecision(1)
            << TMath::Prob(lindhard.chi2, lindhard.ndf) << std::fixed
            << (lindhard.valid ? "" : "  NOT CONVERGED") << std::endl;
  Double_t scale = TMath::Sqrt(excess.chi2 / excess.ndf);
  std::cout << "  Excess ionization:  kappa " << std::setprecision(4)
            << excess.kappa << " +/- " << excess.kappa_err << "  eta* "
            << std::setprecision(2) << excess.eta << " +/- " << excess.eta_err
            << "  chi2 " << std::setprecision(1) << excess.chi2 << " / "
            << excess.ndf << "  p " << std::scientific << std::setprecision(1)
            << TMath::Prob(excess.chi2, excess.ndf) << std::fixed
            << (excess.valid ? "" : "  NOT CONVERGED") << std::endl;
  Double_t delta = lindhard.chi2 - excess.chi2;
  std::cout << "  Delta chi2 " << std::setprecision(1) << delta << " ("
            << TMath::Sqrt(TMath::Max(delta, 0.0))
            << " sigma for one parameter); eta* error scaled by "
               "sqrt(chi2/ndf) = "
            << std::setprecision(2) << scale << ": +/- "
            << excess.eta_err * scale << std::endl;
  std::cout << "  chi2 by dataset (Lindhard / excess):";
  for (size_t d = 0; d < excess.dataset_names.size(); d++) {
    Double_t lindhard_chi2 = 0;
    for (size_t e = 0; e < lindhard.dataset_names.size(); e++)
      if (lindhard.dataset_names[e] == excess.dataset_names[d])
        lindhard_chi2 = lindhard.dataset_chi2[e];
    std::cout << (d == 0 ? " " : ", ") << excess.dataset_names[d] << " "
              << std::setprecision(1) << lindhard_chi2 << " / "
              << excess.dataset_chi2[d] << " (" << excess.dataset_points[d]
              << ")";
  }
  if (excess.n_fitted_nuisances > 0) {
    Double_t lindhard_constraints = 0, excess_constraints = 0;
    for (size_t k = 0; k < data.nuisances.size(); k++) {
      if (!excess.fitted[k])
        continue;
      lindhard_constraints += lindhard.pulls[k] * lindhard.pulls[k];
      excess_constraints += excess.pulls[k] * excess.pulls[k];
    }
    std::cout << ", pull constraints " << std::setprecision(1)
              << lindhard_constraints << " / " << excess_constraints;
  }
  std::cout << std::endl;
  if (excess.n_fitted_nuisances > 0) {
    std::cout << "  pulls (Lindhard / excess):";
    Bool_t first = kTRUE;
    for (size_t k = 0; k < data.nuisances.size(); k++) {
      if (!excess.fitted[k])
        continue;
      std::cout << (first ? " " : ", ") << data.nuisances[k] << " "
                << std::setprecision(2) << lindhard.pulls[k] << " / "
                << excess.pulls[k];
      first = kFALSE;
    }
    std::cout << std::endl;
  }
}

void PrintPoints(const YieldData &data) {
  std::cout << "Points (yield and uncorrelated errors in percent; common "
               "effects per unit pull)"
            << std::endl;
  std::cout << std::fixed;
  for (size_t i = 0; i < data.points.size(); i++) {
    const YieldPoint &point = data.points[i];
    std::cout << "POINT " << std::left << std::setw(23) << point.dataset
              << std::right << (point.detector == kHPGe ? " HPGe  " : " GeLi  ")
              << std::setprecision(4) << std::setw(9) << point.recoil_kev
              << " keV  " << std::setprecision(2) << std::setw(6)
              << 100.0 * point.yield << " +" << 100.0 * point.error_up << " -"
              << 100.0 * point.error_down;
    if (point.recoil_error_kev > 0)
      std::cout << "  recoil +/-" << std::setprecision(3)
                << point.recoil_error_kev << " keV";
    for (size_t c = 0; c < point.couplings.size(); c++) {
      const Coupling &coupling = point.couplings[c];
      std::cout << "  [" << data.nuisances[coupling.nuisance];
      if (coupling.recoil_kev != 0)
        std::cout << ": recoil " << std::setprecision(4) << coupling.recoil_kev
                  << " keV";
      if (coupling.ionization_kev != 0)
        std::cout << ": ionization " << std::setprecision(1)
                  << 1000.0 * coupling.ionization_kev << " eV";
      std::cout << "]";
    }
    std::cout << std::endl;
  }
}

struct DatasetStyle {
  const char *dataset;
  const char *label;
  const char *color_hex;
  Int_t marker;
  Double_t marker_size;
};

// Ge(Li) data in open markers, HPGe data in filled ones; a colour of its own
// for every dataset, none in the grey, azure and red of the curves.
const Int_t N_STYLES = 12;
const DatasetStyle STYLES[N_STYLES] = {
    {"Chasman 1965", "Chasman 1965, Ge(Li)", "#ff7f0e", 24, 1.0},
    {"Sattler 1966", "Sattler 1966, Ge(Li)", "#2ca02c", 26, 1.0},
    {"Chasman 1967", "Chasman 1967, Ge(Li)", "#9467bd", 25, 1.0},
    {"Chasman 1968", "Chasman 1968, Ge(Li)", "#e377c2", 32, 1.0},
    {"Jones 1971", "Jones 1971, Ge(Li)", "#17becf", 27, 1.3},
    {"Jones 1975", "Jones 1975 (n,#gamma), Ge(Li)", "#bcbd22", 30, 1.6},
    {"Baudis 1998", "Baudis 1998, HPGe", "#1f3a93", 20, 1.0},
    {"CoGeNT", "Barbeau 2007, HPGe", "#8c564b", 21, 1.0},
    {"Collar 2021 FeF", "Collar 2021 Fe filter, HPGe", "#d4a017", 22, 1.2},
    {"Collar 2021 capture", "Collar 2021 (n,#gamma), HPGe", "#6b0f1a", 33, 1.6},
    {"Bonhomme 2022", "Bonhomme 2022, HPGe", "#00695c", 23, 1.2},
    {"Kavner 2024 capture", "Kavner 2024 (n,#gamma), HPGe", "#000000", 29,
     1.9}};

// One graph per dataset among the selected points, in STYLES' order, with
// their legend labels, and the axis ranges that hold them.
std::vector<TGraphAsymmErrors *>
DatasetGraphs(const YieldData &data,
              const std::vector<const YieldPoint *> &selected,
              std::vector<TString> &labels, Double_t &x_lo, Double_t &x_hi,
              Double_t &y_lo, Double_t &y_hi) {
  std::vector<TGraphAsymmErrors *> graphs;
  x_lo = 1e9;
  x_hi = 0;
  y_lo = 1e9;
  y_hi = 0;
  for (Int_t s = 0; s < N_STYLES; s++) {
    TString dataset = STYLES[s].dataset;
    Bool_t in_selection = kFALSE;
    for (size_t i = 0; i < selected.size(); i++)
      if (selected[i]->dataset == dataset)
        in_selection = kTRUE;
    if (!in_selection)
      continue;
    TGraphAsymmErrors *graph = new TGraphAsymmErrors();
    for (size_t i = 0; i < data.points.size(); i++) {
      const YieldPoint &point = data.points[i];
      if (point.dataset != dataset)
        continue;
      Int_t n = graph->GetN();
      graph->SetPoint(n, point.recoil_kev, 100.0 * point.yield);
      graph->SetPointError(n, 0, 0, 100.0 * point.error_down,
                           100.0 * point.error_up);
      x_lo = TMath::Min(x_lo, point.recoil_kev);
      x_hi = TMath::Max(x_hi, point.recoil_kev);
      y_lo = TMath::Min(y_lo, 100.0 * (point.yield - point.error_down));
      y_hi = TMath::Max(y_hi, 100.0 * (point.yield + point.error_up));
    }
    Int_t color = TColor::GetColor(STYLES[s].color_hex);
    graph->SetMarkerStyle(STYLES[s].marker);
    graph->SetMarkerSize(STYLES[s].marker_size);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineWidth(PlottingUtils::GetLineWidth());
    graphs.push_back(graph);
    labels.push_back(STYLES[s].label);
  }
  x_lo = 0.7 * x_lo;
  x_hi = 1.4 * x_hi;
  y_lo = TMath::Max(0.0, 5.0 * TMath::Floor((y_lo - 2.0) / 5.0));
  y_hi = 5.0 * TMath::Ceil((y_hi + 0.08 * (y_hi - y_lo)) / 5.0);
  return graphs;
}

// A model curve on a logarithmic grid, yield in percent.
TGraph *YieldCurve(Double_t kappa, Double_t eta, Double_t x_lo, Double_t x_hi,
                   Int_t color) {
  const Int_t n = 400;
  TGraph *curve = new TGraph(n);
  for (Int_t i = 0; i < n; i++) {
    Double_t recoil_kev =
        x_lo * TMath::Power(x_hi / x_lo, (Double_t)i / (n - 1));
    curve->SetPoint(i, recoil_kev,
                    100.0 * IonizationKeV(recoil_kev, kappa, eta) / recoil_kev);
  }
  curve->SetLineColor(color);
  curve->SetLineWidth(PlottingUtils::GetLineWidth());
  return curve;
}

TGraph *PullGraph(const std::vector<const YieldPoint *> &points,
                  const YieldChi2 &chi2, const FitOutcome &outcome) {
  TGraph *pulls = new TGraph();
  for (size_t i = 0; i < points.size(); i++)
    pulls->SetPoint(pulls->GetN(), points[i]->recoil_kev,
                    chi2.PointPull(*points[i], &outcome.best[0]));
  return pulls;
}

// The data of one fit configuration, the theoretical Lindhard curve, both
// fits, and the pulls of both fits beneath.
void DrawYieldFigure(const YieldData &data, const FitConfig &config,
                     const FitOutcome &lindhard, const FitOutcome &excess,
                     const TString &output_name) {
  std::vector<const YieldPoint *> fitted;
  for (size_t i = 0; i < data.points.size(); i++)
    if (Selected(data.points[i], config))
      fitted.push_back(&data.points[i]);

  std::vector<TString> labels;
  Double_t x_lo = 0, x_hi = 0, y_lo = 0, y_hi = 0;
  std::vector<TGraphAsymmErrors *> graphs =
      DatasetGraphs(data, fitted, labels, x_lo, x_hi, y_lo, y_hi);

  TGraph *theory = YieldCurve(0.157, 0, x_lo, x_hi, kGray + 2);
  TGraph *lindhard_curve =
      YieldCurve(lindhard.kappa, 0, x_lo, x_hi, kAzure + 2);
  TGraph *excess_curve =
      YieldCurve(excess.kappa, excess.eta, x_lo, x_hi, kRed + 1);
  // Legend rows: the datasets, then each curve with its numbers on rows of
  // their own beneath it, so the legend is no wider than the dataset names.
  std::vector<TObject *> row_objects;
  std::vector<TString> row_options;
  for (size_t g = 0; g < graphs.size(); g++) {
    row_objects.push_back(graphs[g]);
    row_options.push_back("pe");
  }
  labels.push_back("Lindhard, #kappa = 0.157");
  row_objects.push_back(theory);
  row_options.push_back("l");
  labels.push_back("Lindhard fit");
  row_objects.push_back(lindhard_curve);
  row_options.push_back("l");
  labels.push_back(Form("#kappa = %.4f, #chi^{2}/ndf = %.3f", lindhard.kappa,
                        lindhard.chi2 / lindhard.ndf));
  row_objects.push_back(nullptr);
  row_options.push_back("");
  labels.push_back("Excess-ionization fit");
  row_objects.push_back(excess_curve);
  row_options.push_back("l");
  labels.push_back(Form("#kappa = %.4f, #eta* = %.1f #pm %.1f", excess.kappa,
                        excess.eta, excess.eta_err));
  row_objects.push_back(nullptr);
  row_options.push_back("");
  labels.push_back(Form("#chi^{2}/ndf = %.3f", excess.chi2 / excess.ndf));
  row_objects.push_back(nullptr);
  row_options.push_back("");

  const Int_t plot_h = 800, panel_h = 300;
  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels),
      plot_h + panel_h, kFALSE, YIELD_PLOT_WIDTH_PX);
  TPad *upper = nullptr, *lower = nullptr;
  PlottingUtils::SplitPadForPanel(plot, panel_h, upper, lower);
  upper->cd();
  upper->SetLogx();
  TH1D *frame = new TH1D("yield_frame_" + config.key, "", 1000, x_lo, x_hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(
      frame, kBlack, ";Nuclear recoil energy [keV];Ionization yield [%]");
  frame->SetMinimum(y_lo);
  frame->SetMaximum(y_hi);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw("AXIS");
  theory->Draw("L SAME");
  lindhard_curve->Draw("L SAME");
  excess_curve->Draw("L SAME");
  for (size_t g = 0; g < graphs.size(); g++)
    graphs[g]->Draw("P SAME");

  TLegend *legend =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  for (size_t r = 0; r < labels.size(); r++)
    legend->AddEntry(row_objects[r], labels[r], row_options[r]);

  lower->cd();
  lower->SetLogx();
  std::vector<const YieldPoint *> selected = fitted;
  YieldChi2 chi2(selected, (Int_t)data.nuisances.size(),
                 data.covariance_nuisance, config.correlated);
  TGraph *excess_pulls = PullGraph(selected, chi2, excess);
  TGraph *lindhard_pulls = PullGraph(selected, chi2, lindhard);
  TH1D *pull_frame = FigureStyle::DrawPullPanel(excess_pulls, x_lo, x_hi,
                                                "Nuclear recoil energy [keV]");
  pull_frame->GetXaxis()->SetNoExponent();
  excess_pulls->SetMarkerColor(kRed + 1);
  excess_pulls->SetLineColor(kRed + 1);
  lindhard_pulls->SetMarkerStyle(24);
  lindhard_pulls->SetMarkerSize(0.8);
  lindhard_pulls->SetMarkerColor(kAzure + 2);
  lindhard_pulls->Draw("P SAME");
  excess_pulls->Draw("P SAME");
  PlottingUtils::ScaleFigure(plot);
  // Room between the last x label and the legend, the same in both pads.
  upper->SetRightMargin(upper->GetRightMargin() + 0.03);
  lower->SetRightMargin(lower->GetRightMargin() + 0.03);
  FigureStyle::FitTitleToLabels(lower, pull_frame,
                                FigureStyle::PullPanelLabels());
  std::vector<TString> yield_labels;
  yield_labels.push_back(Form("%.0f", y_hi));
  FigureStyle::FitTitleToLabels(upper, frame, yield_labels);
  PlottingUtils::SaveFigure(canvas, output_name, "", PlotSaveOptions::kLINEAR);

  for (size_t i = 0; i < selected.size(); i++) {
    Double_t pull = excess_pulls->GetPointY((Int_t)i);
    if (TMath::Abs(pull) > 5.5)
      std::cout << "  off the pull panel (excess fit): " << selected[i]->dataset
                << " " << selected[i]->recoil_kev << " keV, pull "
                << std::setprecision(1) << pull << std::endl;
  }
}

// A capture-line point as its paper published it.
void SetPublished(YieldPoint &point, Double_t ionization_kev,
                  Double_t ionization_err_kev, Double_t recoil_kev) {
  point.recoil_kev = recoil_kev;
  point.yield = ionization_kev / recoil_kev;
  point.error_up = ionization_err_kev / recoil_kev;
  point.error_down = point.error_up;
}

// The measurements alone with Lindhard's prediction (kappa = 0.157), for
// introducing them before this work: the points of one fit configuration,
// the capture-line points as published, no fits, no pulls.
void DrawYieldDataFigure(const YieldData &data, const FitConfig &config,
                         const TString &output_name) {
  using namespace GeYieldData;
  YieldData published = data;
  for (size_t i = 0; i < published.points.size(); i++) {
    YieldPoint &point = published.points[i];
    if (point.dataset == "Jones 1975")
      SetPublished(point, JONES1975_PUBLISHED_IONIZATION_KEV,
                   JONES1975_PUBLISHED_IONIZATION_ERR_KEV,
                   JONES1975_PUBLISHED_RECOIL_KEV);
    if (point.dataset == "Collar 2021 capture")
      SetPublished(point, COLLAR2021_PUBLISHED_IONIZATION_KEV,
                   COLLAR2021_PUBLISHED_IONIZATION_ERR_KEV,
                   COLLAR2021_PUBLISHED_RECOIL_KEV);
    if (point.dataset == "Kavner 2024 capture")
      SetPublished(point, KAVNER2024_PUBLISHED_IONIZATION_KEV,
                   KAVNER2024_PUBLISHED_IONIZATION_ERR_KEV,
                   KAVNER2024_PUBLISHED_RECOIL_KEV);
  }
  std::vector<const YieldPoint *> selected;
  for (size_t i = 0; i < published.points.size(); i++)
    if (Selected(published.points[i], config))
      selected.push_back(&published.points[i]);
  std::vector<TString> labels;
  Double_t x_lo = 0, x_hi = 0, y_lo = 0, y_hi = 0;
  std::vector<TGraphAsymmErrors *> graphs =
      DatasetGraphs(published, selected, labels, x_lo, x_hi, y_lo, y_hi);
  TGraph *theory = YieldCurve(0.157, 0, x_lo, x_hi, kGray + 2);
  labels.push_back("Lindhard theory, #kappa = 0.157");

  TPad *plot = nullptr, *legend_pad = nullptr;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legend_pad, (Int_t)PlottingUtils::LegendWidthPx(labels), 800,
      kFALSE, YIELD_PLOT_WIDTH_PX);
  plot->cd();
  plot->SetLogx();
  TH1D *frame =
      new TH1D("yield_data_frame_" + config.key, "", 1000, x_lo, x_hi);
  frame->SetDirectory(nullptr);
  frame->SetStats(kFALSE);
  PlottingUtils::ConfigureHistogram(
      frame, kBlack, ";Nuclear recoil energy [keV];Ionization yield [%]");
  frame->SetMinimum(y_lo);
  frame->SetMaximum(y_hi);
  frame->GetXaxis()->SetNoExponent();
  frame->Draw("AXIS");
  theory->Draw("L SAME");
  for (size_t g = 0; g < graphs.size(); g++)
    graphs[g]->Draw("P SAME");
  TLegend *legend =
      PlottingUtils::AddSideLegend(plot, legend_pad, (Int_t)labels.size());
  for (size_t g = 0; g < graphs.size(); g++)
    legend->AddEntry(graphs[g], labels[g], "pe");
  legend->AddEntry(theory, labels.back(), "l");
  PlottingUtils::ScaleFigure(plot);
  // Room between the last x label and the legend.
  plot->SetRightMargin(plot->GetRightMargin() + 0.03);
  std::vector<TString> yield_labels;
  yield_labels.push_back(Form("%.0f", y_hi));
  FigureStyle::FitTitleToLabels(plot, frame, yield_labels);
  PlottingUtils::SaveFigure(canvas, output_name, "", PlotSaveOptions::kLINEAR);
}

// A dataset's tension with a best fit: the chi2 of its points plus the
// constraint terms of the nuisances that act on it alone (so that an offset
// a dataset's own nuisance absorbs still counts), and the p-value of that
// chi2 for its number of points.
void DatasetTension(const YieldData &data, const FitConfig &config,
                    const FitOutcome &fit, Int_t dataset_index, Double_t &chi2,
                    Double_t &p) {
  const TString &dataset = fit.dataset_names[dataset_index];
  chi2 = fit.dataset_chi2[dataset_index];
  for (size_t k = 0; k < data.nuisances.size(); k++) {
    if (!fit.fitted[k])
      continue;
    Bool_t acts_here = kFALSE, acts_elsewhere = kFALSE;
    for (size_t i = 0; i < data.points.size(); i++) {
      const YieldPoint &point = data.points[i];
      if (!Selected(point, config))
        continue;
      for (size_t c = 0; c < point.couplings.size(); c++) {
        if (point.couplings[c].nuisance != (Int_t)k)
          continue;
        if (point.dataset == dataset)
          acts_here = kTRUE;
        else
          acts_elsewhere = kTRUE;
      }
    }
    if (acts_here && !acts_elsewhere)
      chi2 += fit.pulls[k] * fit.pulls[k];
  }
  p = TMath::Prob(chi2, fit.dataset_points[dataset_index]);
}

struct AdaptiveOutcome {
  FitConfig config;
  FitOutcome fit;
  std::vector<TString> removed;
};

// Kavner et al.'s adaptive fit: while the fit's p-value is below 0.05, the
// dataset in the strongest tension with the best fit (the smallest
// DatasetTension p-value) is removed and the rest refitted.
AdaptiveOutcome RunAdaptiveFit(const YieldData &data, const FitConfig &start,
                               Bool_t excess_ionization) {
  const Double_t MIN_P = 0.05;
  AdaptiveOutcome outcome;
  outcome.config = start;
  while (kTRUE) {
    outcome.fit = RunFit(data, outcome.config, excess_ionization);
    Double_t fit_p = TMath::Prob(outcome.fit.chi2, outcome.fit.ndf);
    if (fit_p >= MIN_P || outcome.fit.dataset_names.size() <= 1)
      return outcome;
    Int_t worst = -1;
    Double_t worst_p = 2, worst_chi2 = 0;
    for (size_t d = 0; d < outcome.fit.dataset_names.size(); d++) {
      Double_t chi2 = 0, p = 0;
      DatasetTension(data, outcome.config, outcome.fit, (Int_t)d, chi2, p);
      if (p < worst_p) {
        worst = (Int_t)d;
        worst_p = p;
        worst_chi2 = chi2;
      }
    }
    const TString name = outcome.fit.dataset_names[worst];
    std::cout << "  fit chi2/ndf " << std::fixed << std::setprecision(3)
              << outcome.fit.chi2 / outcome.fit.ndf << " (p " << std::scientific
              << std::setprecision(1) << fit_p << "): remove " << name
              << ", chi2 " << std::fixed << std::setprecision(1) << worst_chi2
              << " for " << outcome.fit.dataset_points[worst] << " points (p "
              << std::scientific << std::setprecision(1) << worst_p << ")"
              << std::fixed << std::endl;
    outcome.removed.push_back(name);
    std::vector<TString> kept;
    for (size_t d = 0; d < outcome.config.datasets.size(); d++)
      if (outcome.config.datasets[d] != name)
        kept.push_back(outcome.config.datasets[d]);
    outcome.config.datasets = kept;
  }
}

TString JoinNames(const std::vector<TString> &names) {
  TString joined = names.empty() ? TString("none") : TString("");
  for (size_t n = 0; n < names.size(); n++)
    joined += (n == 0 ? "" : ", ") + names[n];
  return joined;
}

// One row of the dissertation's table: both models on one set of points.
struct TableRow {
  TString label;
  FitOutcome lindhard;
  FitOutcome excess;
};

TableRow MakeTableRow(const TString &label, const FitOutcome &lindhard,
                      const FitOutcome &excess) {
  TableRow row;
  row.label = label;
  row.lindhard = lindhard;
  row.excess = excess;
  return row;
}

// A dataset as the table's notes name it.
TString TableName(const TString &dataset) {
  if (dataset == "Collar 2021 FeF")
    return "Collar 2021 iron filter";
  return dataset;
}

// The adaptive subset of the excess-ionization model, both models fitted on
// it (returned as adaptive_row), and the Lindhard model's own adaptive subset.
void RunAdaptiveFits(const YieldData &data, const FitConfig &config,
                     std::ofstream &result, TableRow &adaptive_row) {
  std::cout << std::endl;
  std::cout << "ADAPTIVE " << config.key
            << ", excess-ionization model:" << std::endl;
  AdaptiveOutcome excess = RunAdaptiveFit(data, config, kTRUE);
  FitConfig subset = excess.config;
  subset.key = config.key + "_adaptive";
  subset.title = config.title +
                 ", excess-ionization adaptive subset "
                 "(removed: " +
                 JoinNames(excess.removed) + ")";
  FitOutcome lindhard_on_subset = RunFit(data, subset, kFALSE);
  PrintFit(data, subset, lindhard_on_subset, excess.fit);
  std::vector<TString> removed_names;
  for (size_t r = 0; r < excess.removed.size(); r++)
    removed_names.push_back(TableName(excess.removed[r]));
  adaptive_row =
      MakeTableRow("Adaptive subset", lindhard_on_subset, excess.fit);
  adaptive_row.label += " (without " + JoinNames(removed_names) + ")";

  std::cout << std::endl;
  std::cout << "ADAPTIVE " << config.key << ", Lindhard model:" << std::endl;
  AdaptiveOutcome lindhard = RunAdaptiveFit(data, config, kFALSE);
  std::cout << "  Lindhard adaptive subset (removed: "
            << JoinNames(lindhard.removed) << "): kappa "
            << std::setprecision(4) << lindhard.fit.kappa << " +/- "
            << lindhard.fit.kappa_err << ", chi2/ndf " << std::setprecision(3)
            << lindhard.fit.chi2 / lindhard.fit.ndf << " (p "
            << std::setprecision(2)
            << TMath::Prob(lindhard.fit.chi2, lindhard.fit.ndf) << ")"
            << std::endl;
  std::cout << "  its pulls beyond 2:";
  for (size_t k = 0; k < data.nuisances.size(); k++)
    if (lindhard.fit.fitted[k] && TMath::Abs(lindhard.fit.pulls[k]) > 2)
      std::cout << " " << data.nuisances[k] << " " << std::setprecision(2)
                << lindhard.fit.pulls[k];
  std::cout << std::endl;

  const TString &key = subset.key;
  result << "# " << key
         << " removed (excess-ionization model): " << JoinNames(excess.removed)
         << std::endl;
  result << "# " << key
         << " removed (Lindhard model): " << JoinNames(lindhard.removed)
         << std::endl;
  result << std::fixed << std::setprecision(5);
  result << key << "_h1_kappa " << excess.fit.kappa << std::endl;
  result << key << "_h1_kappa_err " << excess.fit.kappa_err << std::endl;
  result << std::setprecision(3);
  result << key << "_h1_eta " << excess.fit.eta << std::endl;
  result << key << "_h1_eta_err " << excess.fit.eta_err << std::endl;
  result << std::setprecision(2);
  result << key << "_h1_chi2 " << excess.fit.chi2 << std::endl;
  result << key << "_h1_ndf " << excess.fit.ndf << std::endl;
  result << key << "_h0_chi2 " << lindhard_on_subset.chi2 << std::endl;
  result << std::setprecision(5);
  result << key << "_lindhard_kappa " << lindhard.fit.kappa << std::endl;
  result << key << "_lindhard_kappa_err " << lindhard.fit.kappa_err
         << std::endl;
  result << std::setprecision(2);
  result << key << "_lindhard_chi2 " << lindhard.fit.chi2 << std::endl;
  result << key << "_lindhard_ndf " << lindhard.fit.ndf << std::endl;
}

TString TableCells(const TableRow &row) {
  const FitOutcome &lindhard = row.lindhard;
  const FitOutcome &excess = row.excess;
  Double_t delta = lindhard.chi2 - excess.chi2;
  return Form("%d & $%.4f \\pm %.4f$ & %.3f & $%.4f \\pm %.4f$ & $%.1f \\pm "
              "%.1f$ & %.3f & %.1f (%.1f$\\sigma$)",
              excess.n_points, lindhard.kappa, lindhard.kappa_err,
              lindhard.chi2 / lindhard.ndf, excess.kappa, excess.kappa_err,
              excess.eta, excess.eta_err, excess.chi2 / excess.ndf, delta,
              TMath::Sqrt(TMath::Max(delta, 0.0)));
}

// The dissertation's table of the fits: the correlated fits, then the same
// data with every error uncorrelated. A row's label goes in the first
// column, an adaptive subset's left-out datasets in a note for the caption.
void WriteTable(const TString &path, const std::vector<TableRow> &correlated,
                const std::vector<TableRow> &uncorrelated) {
  std::ofstream out(path.Data());
  if (!out) {
    std::cerr << "ERROR: cannot write " << path << std::endl;
    return;
  }
  out << "% Written by UM-ANSG 73mGe/evaluation/ModifiedLindhard.cpp; edit "
         "that, not this file."
      << std::endl;
  out << "% For the caption: N is the number of points; Lindhard fits kappa "
         "alone, the excess-ionization "
         "model kappa and eta* (electrons, W = 2.96 eV); Delta chi2 is the "
         "Lindhard chi2 minus the excess-ionization chi2 on the same points, "
         "its significance sqrt(Delta chi2). Correlated: each source's common "
         "uncertainties as shared nuisance parameters, the 68.75 keV gamma "
         "energy as a covariance. Uncorrelated: every error per point, as "
         "Kavner et al. fit."
      << std::endl;
  for (size_t r = 0; r < correlated.size(); r++)
    if (correlated[r].label.Contains("("))
      out << "% " << correlated[r].label << std::endl;
  out << "\\small" << std::endl;
  out << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  out << "\\begin{tabular}{@{}lccccccc@{}}" << std::endl;
  out << "\\toprule" << std::endl;
  out << " & & \\multicolumn{2}{c}{Lindhard} & "
         "\\multicolumn{3}{c}{Excess ionization} & \\\\"
      << std::endl;
  out << "\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}" << std::endl;
  out << "Data & $N$ & $\\kappa$ & $\\chi^2/\\mathrm{ndf}$ & $\\kappa$ & "
         "$\\eta^\\star$ & $\\chi^2/\\mathrm{ndf}$ & $\\Delta\\chi^2$ \\\\"
      << std::endl;
  out << "\\midrule" << std::endl;
  out << "\\multicolumn{8}{@{}l}{\\textit{Correlated}} \\\\" << std::endl;
  for (size_t r = 0; r < correlated.size(); r++) {
    TString label = correlated[r].label;
    if (label.Contains("("))
      label = label(0, label.Index(" ("));
    out << label << " & " << TableCells(correlated[r]) << " \\\\" << std::endl;
  }
  out << "\\midrule" << std::endl;
  out << "\\multicolumn{8}{@{}l}{\\textit{Uncorrelated}} \\\\" << std::endl;
  for (size_t r = 0; r < uncorrelated.size(); r++)
    out << uncorrelated[r].label << " & " << TableCells(uncorrelated[r])
        << " \\\\" << std::endl;
  out << "\\bottomrule" << std::endl;
  out << "\\end{tabular}" << std::endl;
  std::cout << "Wrote " << path << std::endl;
}

Bool_t ReadResult(const TString &path,
                  std::map<std::string, Double_t> &values) {
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "ERROR: cannot read " << path << std::endl;
    return kFALSE;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    std::string key;
    Double_t value = 0;
    if (fields >> key >> value)
      values[key] = value;
  }
  return kTRUE;
}

} // namespace

void ModifiedLindhard() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNGAndPDF,
                                project_root + "/plots/evaluation",
                                project_root + "/root_files");
  std::map<std::string, Double_t> recommended, recoil;
  if (!ReadResult(project_root + "/results/recommended_value.result",
                  recommended) ||
      !ReadResult(project_root + "/results/recoil_energy.result", recoil))
    return;

  YieldData data;
  GammaReference gamma;
  gamma.energy_kev = recommended["mean"];
  gamma.error_kev = recommended["mean_err"];
  gamma.jk_shift_kev = recommended["jk_shift"];
  gamma.recoil_kev = recoil["mean_total_recoil"] / 1000.0;
  gamma.nuisance = AddNuisance(data, "68.75 keV gamma energy (covariance)");
  data.covariance_nuisance = gamma.nuisance;
  std::cout << std::fixed << std::setprecision(6) << "Gamma "
            << gamma.energy_kev << " +/- " << gamma.error_kev
            << " keV, Jones and Kraner shift " << gamma.jk_shift_kev
            << " keV, capture recoil " << gamma.recoil_kev << " keV"
            << std::endl;

  using namespace GeYieldData;
  AddYieldTable(data, "Chasman 1965", kGeLi, CHASMAN1965, CHASMAN1965_N);
  AddYieldTable(data, "Sattler 1966", kGeLi, SATTLER1966, SATTLER1966_N);
  AddChasman1967(data);
  AddChasman1968(data);
  AddJones1971(data, gamma);
  AddCaptureLine(data, "Jones 1975", kGeLi, JONES1975_GAMMA_PLUS_RECOIL_KEV,
                 JONES1975_GAMMA_PLUS_RECOIL_ERR_KEV, kTRUE, gamma);
  AddBaudis1998(data);
  AddCoGeNT(data);
  AddCollar2021IronFilter(data);
  AddCaptureLine(data, "Collar 2021 capture", kHPGe,
                 COLLAR2021_GAMMA_PLUS_RECOIL_KEV,
                 COLLAR2021_GAMMA_PLUS_RECOIL_ERR_KEV, kFALSE, gamma);
  AddBonhomme2022(data);
  AddCaptureLine(data, "Kavner 2024 capture", kHPGe,
                 KAVNER2024_GAMMA_PLUS_RECOIL_KEV,
                 KAVNER2024_GAMMA_PLUS_RECOIL_ERR_KEV, kFALSE, gamma);
  PrintPoints(data);

  std::vector<TString> nominal;
  nominal.push_back("Chasman 1965");
  nominal.push_back("Sattler 1966");
  nominal.push_back("Chasman 1967");
  nominal.push_back("Chasman 1968");
  nominal.push_back("Jones 1971");
  nominal.push_back("Jones 1975");
  nominal.push_back("Baudis 1998");
  nominal.push_back("CoGeNT");
  nominal.push_back("Collar 2021 FeF");
  nominal.push_back("Collar 2021 capture");
  nominal.push_back("Bonhomme 2022");
  nominal.push_back("Kavner 2024 capture");

  std::vector<FitConfig> configs;
  FitConfig config;
  config.key = "all";
  config.title = "all data, correlated";
  config.datasets = nominal;
  config.hpge_only = kFALSE;
  config.correlated = kTRUE;
  configs.push_back(config);
  config.key = "hpge";
  config.title = "HPGe data, correlated";
  config.hpge_only = kTRUE;
  configs.push_back(config);
  config.key = "all_uncorrelated";
  config.title = "all data, every error uncorrelated";
  config.hpge_only = kFALSE;
  config.correlated = kFALSE;
  configs.push_back(config);
  config.key = "hpge_uncorrelated";
  config.title = "HPGe data, every error uncorrelated";
  config.hpge_only = kTRUE;
  configs.push_back(config);

  const TString result_path =
      project_root + "/results/modified_lindhard.result";
  std::ofstream result(result_path.Data());
  result << "# ModifiedLindhard.cpp: Lindhard (h0) and excess-ionization (h1) "
            "fits of the Ge ionization yield; eta* in electrons"
         << std::endl;
  std::vector<TableRow> correlated_rows, uncorrelated_rows;
  for (size_t f = 0; f < configs.size(); f++) {
    FitOutcome lindhard = RunFit(data, configs[f], kFALSE);
    FitOutcome excess = RunFit(data, configs[f], kTRUE);
    PrintFit(data, configs[f], lindhard, excess);
    const TString &key = configs[f].key;
    if (key == "all") {
      DrawYieldFigure(data, configs[f], lindhard, excess, "ge_yield_fit_all");
      DrawYieldDataFigure(data, configs[f], "ge_yield_data");
    }
    if (key == "hpge")
      DrawYieldFigure(data, configs[f], lindhard, excess, "ge_yield_fit_hpge");
    TString row_label = configs[f].hpge_only ? "HPGe data" : "All data";
    TableRow adaptive_row;
    if (key == "all" || key == "hpge" || key == "all_uncorrelated")
      RunAdaptiveFits(data, configs[f], result, adaptive_row);
    if (configs[f].correlated) {
      correlated_rows.push_back(MakeTableRow(row_label, lindhard, excess));
      if (key == "all")
        correlated_rows.push_back(adaptive_row);
    } else {
      uncorrelated_rows.push_back(MakeTableRow(row_label, lindhard, excess));
    }
    result << std::fixed << std::setprecision(5);
    result << key << "_points " << excess.n_points << std::endl;
    result << key << "_h0_kappa " << lindhard.kappa << std::endl;
    result << key << "_h0_kappa_err " << lindhard.kappa_err << std::endl;
    result << key << "_h1_kappa " << excess.kappa << std::endl;
    result << key << "_h1_kappa_err " << excess.kappa_err << std::endl;
    result << std::setprecision(3);
    result << key << "_h1_eta " << excess.eta << std::endl;
    result << key << "_h1_eta_err " << excess.eta_err << std::endl;
    result << std::setprecision(2);
    result << key << "_h0_chi2 " << lindhard.chi2 << std::endl;
    result << key << "_h1_chi2 " << excess.chi2 << std::endl;
    result << key << "_h1_ndf " << excess.ndf << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Wrote " << result_path << std::endl;
  WriteTable(TString(gSystem->HomeDirectory()) +
                 "/Dissertation/chapters/chapter_ge_recoils/figures/"
                 "ge_yield_fit_table.tex",
             correlated_rows, uncorrelated_rows);
}
