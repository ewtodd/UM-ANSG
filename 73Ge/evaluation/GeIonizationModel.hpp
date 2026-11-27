// The germanium ionization models that ModifiedLindhard.cpp fits: Lindhard's
// yield in the parameterization of Lewin and Smith, Astropart. Phys. 6, 87
// (1996), and the excess-ionization model of Kavner et al.,
// E_ee = Q_L(kappa, E_nr) E_nr + eta* W.
#ifndef GE_IONIZATION_MODEL_HPP
#define GE_IONIZATION_MODEL_HPP

#include <Rtypes.h>
#include <TMath.h>

namespace GeIonizationModel {

const Double_t GE_Z = 32;
// The pair-creation energy Kavner et al. use for Ge at 77 K.
const Double_t PAIR_ENERGY_KEV = 0.00296;

inline Double_t LindhardQ(Double_t recoil_kev, Double_t kappa) {
  Double_t epsilon = 11.5 * recoil_kev * TMath::Power(GE_Z, -7.0 / 3.0);
  Double_t g = 3.0 * TMath::Power(epsilon, 0.15) +
               0.7 * TMath::Power(epsilon, 0.6) + epsilon;
  return kappa * g / (1.0 + kappa * g);
}

inline Double_t IonizationKeV(Double_t recoil_kev, Double_t kappa,
                              Double_t eta) {
  return LindhardQ(recoil_kev, kappa) * recoil_kev + eta * PAIR_ENERGY_KEV;
}

} // namespace GeIonizationModel

#endif
