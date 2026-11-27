// The recoil of 73Ge after thermal neutron capture on 72Ge, and what it means
// for the ionization measured from the 68.75 keV summed peak: the recoil of
// the two cascades through the 68.75 keV level on current data, beside
// Jones and Kraner's 254.1 eV; how much of the 68.75 keV level's feeding those
// cascades explain; and the Lindhard ionization against the measured one,
// including at the largest recoil kinematics allows.
//
// Reads results/recommended_value.result (RecommendedValue.cpp) for the
// measured ionization; writes results/recoil_energy.result.
//
// Usage: root -b -q 'RecoilEnergy.cpp+'
#include "ProjectPaths.hpp"

#include <TMath.h>
#include <TString.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace {

// Atomic mass unit (CODATA 2018) and the 73Ge atomic mass (AME2020).
const Double_t U_KEV = 931494.10242;
const Double_t GE73_MASS_U = 72.923458956;

// ENSDF, Nuclear Data Sheets for A = 73 (Singh and Chen, Nucl. Data Sheets
// 158, 1 (2019), doi:10.1016/j.nds.2019.02.006): S(n) of 73Ge from the
// Adopted Levels; the primary gamma rays and the intensities per 100 neutron
// captures from 72Ge(n,g) E=thermal (primaries from 1991Is01); the secondary
// gamma rays, branchings and the 68.752 keV conversion coefficient from the
// Adopted Gammas.
const Double_t SN_KEV = 6782.94;
const Double_t GAMMA_68_KEV = 68.752;
const Double_t GAMMA_68_INTENSITY = 3.80; // per 100 captures
const Double_t GAMMA_68_ALPHA = 0.238;

struct Path {
  const char *name;
  Double_t primary_kev;
  Double_t primary_intensity; // per 100 captures
  Double_t secondary_kev[2];
  // The share of the intermediate levels' decays that go on to 68.752 keV.
  Double_t to_68_fraction;
};

// 931.73 keV level: 430.2 keV to 501.49 keV, which sends 432.7 keV (1.60 per
// 100 captures) to 68.752 keV and 488.0 keV (0.9) to 13.28 keV.
// 915.45 keV level: 561.6 keV to 353.70 keV (100) and 915.8 keV to the
// ground state (81), relative intensities from (n,n'g); 353.70 keV decays
// only by 284.8 keV to 68.752 keV.
const Int_t N_PATHS = 2;
const Path PATHS[N_PATHS] = {
    {"931.73 keV path", 5850.97, 2.21, {430.2, 432.7}, 1.60 / (1.60 + 0.9)},
    {"915.45 keV path", 5867.12, 0.75, {561.6, 284.8}, 100.0 / (100.0 + 81.0)}};

// Jones and Kraner, Phys. Rev. A 11, 1347 (1975), Sec. II: their primary
// recoils through the 931.5 and 915.2 keV levels, their total, and their
// Lindhard estimate of the stopping time, t = 9.59e-15 E^(1/6) s, E in eV.
const Double_t JK_PRIMARY_RECOIL_EV[N_PATHS] = {250.3, 251.7};
const Double_t JK_TOTAL_RECOIL_EV = 254.1;
const Double_t JK_STOPPING_COEFFICIENT_S = 9.59e-15;

// Lindhard ionization yield for Ge with the theoretical k = 0.157, in the
// parameterization of Lewin and Smith, Astropart. Phys. 6, 87 (1996).
const Double_t LINDHARD_K = 0.157;
const Double_t GE_Z = 32;

// An illustrative recoil mixture: a share of the 68.75 keV decays following
// other cascades with a lower mean recoil, and a typical HPGe width at 69 keV.
const Double_t MIXTURE_SHARE = 0.4;
const Double_t MIXTURE_RECOIL_EV = 235.0;
const Double_t HPGE_SIGMA_EV = 200.0;

Double_t RecoilEV(Double_t gamma_kev) {
  return 1000.0 * gamma_kev * gamma_kev / (2.0 * GE73_MASS_U * U_KEV);
}

Double_t LindhardQ(Double_t recoil_ev) {
  Double_t epsilon =
      11.5 * (recoil_ev / 1000.0) * TMath::Power(GE_Z, -7.0 / 3.0);
  Double_t g = 3.0 * TMath::Power(epsilon, 0.15) +
               0.7 * TMath::Power(epsilon, 0.6) + epsilon;
  return LINDHARD_K * g / (1.0 + LINDHARD_K * g);
}

// The recoil whose Lindhard ionization equals ionization_ev, by bisection.
Double_t LindhardRecoilForIonizationEV(Double_t ionization_ev) {
  Double_t low = 1.0, high = 10000.0;
  for (Int_t i = 0; i < 100; i++) {
    Double_t middle = 0.5 * (low + high);
    if (LindhardQ(middle) * middle < ionization_ev)
      low = middle;
    else
      high = middle;
  }
  return 0.5 * (low + high);
}

Bool_t ReadResult(const TString &path,
                  std::map<std::string, Double_t> &values) {
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "ERROR: cannot read " << path
              << "; run RecommendedValue.cpp first." << std::endl;
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

void RecoilEnergy() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  std::map<std::string, Double_t> recommended;
  if (!ReadResult(project_root + "/results/recommended_value.result",
                  recommended))
    return;
  Double_t ionization_ev = recommended["ionization_ev"];
  Double_t ionization_err_ev = recommended["ionization_err_ev"];

  std::cout << std::fixed;
  std::cout << "Recoil, E^2/2Mc^2 with M = M(73Ge):" << std::endl;
  Double_t primary_recoil[N_PATHS], total_recoil[N_PATHS];
  Double_t primary_weight[N_PATHS], feeding_weight[N_PATHS];
  Double_t primary_sum = 0, feeding_sum = 0;
  Double_t final_recoil = RecoilEV(GAMMA_68_KEV);
  for (Int_t p = 0; p < N_PATHS; p++) {
    primary_recoil[p] = RecoilEV(PATHS[p].primary_kev);
    Double_t secondary = final_recoil;
    for (Int_t s = 0; s < 2; s++)
      secondary += RecoilEV(PATHS[p].secondary_kev[s]);
    total_recoil[p] = primary_recoil[p] + secondary;
    primary_weight[p] = PATHS[p].primary_intensity;
    feeding_weight[p] = PATHS[p].primary_intensity * PATHS[p].to_68_fraction;
    primary_sum += primary_weight[p];
    feeding_sum += feeding_weight[p];
    std::cout << "  " << PATHS[p].name << ": primary " << std::setprecision(2)
              << PATHS[p].primary_kev << " keV -> " << primary_recoil[p]
              << " eV (Jones and Kraner " << std::setprecision(1)
              << JK_PRIMARY_RECOIL_EV[p] << "); later gamma rays "
              << std::setprecision(2) << secondary << " eV (largest "
              << RecoilEV(PATHS[p].secondary_kev[0]) << ", "
              << RecoilEV(PATHS[p].secondary_kev[1]) << "); total "
              << total_recoil[p] << " eV" << std::endl;
  }

  // Two weightings: by the primary intensities (the shares in the decay
  // scheme figure), and by how often each path reaches the 68.75 keV level.
  Double_t primary_share = primary_weight[0] / primary_sum;
  Double_t feeding_share = feeding_weight[0] / feeding_sum;
  Double_t mean_primary = 0, mean_total = 0;
  for (Int_t p = 0; p < N_PATHS; p++) {
    mean_primary += feeding_weight[p] / feeding_sum * primary_recoil[p];
    mean_total += feeding_weight[p] / feeding_sum * total_recoil[p];
  }
  std::cout << "Shares of the two paths: primary intensities "
            << std::setprecision(1) << 100.0 * primary_share << " / "
            << 100.0 * (1.0 - primary_share) << " %, reaching 68.75 keV "
            << 100.0 * feeding_share << " / " << 100.0 * (1.0 - feeding_share)
            << " %" << std::endl;
  std::cout << "Mean recoil of the 68.75 keV decays through the two paths: "
            << "primary " << std::setprecision(2) << mean_primary
            << " eV, with the later gamma rays " << mean_total
            << " eV (Jones and Kraner " << std::setprecision(1)
            << JK_TOTAL_RECOIL_EV << " eV)" << std::endl;

  // How much of the 68.75 keV level's population the two paths explain.
  Double_t decays_68 = GAMMA_68_INTENSITY * (1.0 + GAMMA_68_ALPHA);
  std::cout << "68.75 keV level: " << std::setprecision(2) << decays_68
            << " decays per 100 captures; the two paths feed " << feeding_sum
            << " (" << std::setprecision(0) << 100.0 * feeding_sum / decays_68
            << " %)" << std::endl;

  Double_t stopping_time_s =
      JK_STOPPING_COEFFICIENT_S * TMath::Power(mean_total, 1.0 / 6.0);
  std::cout << "Lindhard stopping time (Jones and Kraner): " << std::scientific
            << std::setprecision(1) << stopping_time_s << " s" << std::fixed
            << std::endl;

  // Lindhard against the measurement, up to the largest recoil possible:
  // the whole of S(n) in a single gamma ray.
  Double_t max_recoil = RecoilEV(SN_KEV);
  Double_t needed_recoil = LindhardRecoilForIonizationEV(ionization_ev);
  std::cout << std::endl
            << "Measured ionization (RecommendedValue.cpp): "
            << std::setprecision(1) << ionization_ev << " +/- "
            << ionization_err_ev << " eV" << std::endl;
  const Int_t N_RECOILS = 4;
  const char *const recoil_names[N_RECOILS] = {
      "Jones and Kraner", "two paths, current data", "illustrative mixture",
      "kinematic maximum"};
  Double_t mixture_recoil =
      (1.0 - MIXTURE_SHARE) * mean_total + MIXTURE_SHARE * MIXTURE_RECOIL_EV;
  Double_t recoils[N_RECOILS] = {JK_TOTAL_RECOIL_EV, mean_total, mixture_recoil,
                                 max_recoil};
  for (Int_t r = 0; r < N_RECOILS; r++) {
    Double_t q = LindhardQ(recoils[r]);
    std::cout << "  " << std::left << std::setw(26) << recoil_names[r]
              << std::right << std::setprecision(1) << recoils[r]
              << " eV: Lindhard Q " << std::setprecision(2) << 100.0 * q
              << " %, ionization " << std::setprecision(1) << q * recoils[r]
              << " eV; measured quenching factor " << std::setprecision(1)
              << 100.0 * ionization_ev / recoils[r] << " +/- "
              << 100.0 * ionization_err_ev / recoils[r] << " %" << std::endl;
  }
  Double_t max_lindhard = LindhardQ(max_recoil) * max_recoil;
  std::cout << "Measured over Lindhard at the kinematic maximum: "
            << std::setprecision(1)
            << (ionization_ev - max_lindhard) / ionization_err_ev
            << " sigma; Lindhard needs " << std::setprecision(0)
            << needed_recoil << " eV of recoil for the measured ionization"
            << std::endl;

  // The mixture's added width against a typical HPGe peak, both recoils
  // taken at the measured quenching factor.
  Double_t ionization_gap =
      ionization_ev / mean_total * (mean_total - MIXTURE_RECOIL_EV);
  Double_t mixture_rms =
      TMath::Sqrt(MIXTURE_SHARE * (1.0 - MIXTURE_SHARE)) * ionization_gap;
  Double_t broadening =
      TMath::Sqrt(HPGE_SIGMA_EV * HPGE_SIGMA_EV + mixture_rms * mixture_rms) -
      HPGE_SIGMA_EV;
  std::cout << "Illustrative mixture (" << std::setprecision(0)
            << 100.0 * MIXTURE_SHARE << " % at " << MIXTURE_RECOIL_EV
            << " eV): ionization spread " << std::setprecision(1) << mixture_rms
            << " eV rms, widening a " << std::setprecision(0) << HPGE_SIGMA_EV
            << " eV sigma by " << std::setprecision(3) << broadening << " eV"
            << std::endl;

  std::ofstream result((project_root + "/results/recoil_energy.result").Data());
  result << "# RecoilEnergy.cpp: 73Ge recoil after thermal capture on 72Ge "
            "(ENSDF 2019) and the Lindhard comparison; eV unless noted"
         << std::endl;
  result << std::fixed << std::setprecision(3);
  result << "primary_recoil_931 " << primary_recoil[0] << std::endl;
  result << "primary_recoil_915 " << primary_recoil[1] << std::endl;
  result << "primary_share_931 " << primary_share << std::endl;
  result << "feeding_share_931 " << feeding_share << std::endl;
  result << "mean_primary_recoil " << mean_primary << std::endl;
  result << "mean_total_recoil " << mean_total << std::endl;
  result << "feeding_fraction_of_68 " << feeding_sum / decays_68 << std::endl;
  result << "kinematic_max_recoil " << max_recoil << std::endl;
  result << "lindhard_ionization_at_max " << max_lindhard << std::endl;
  result << "lindhard_recoil_needed " << needed_recoil << std::endl;
  result << "quenching_factor_percent " << 100.0 * ionization_ev / mean_total
         << std::endl;
  result << "quenching_factor_err_percent "
         << 100.0 * ionization_err_ev / mean_total << std::endl;
}
