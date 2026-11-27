#ifndef REFERENCE_LINES_HPP
#define REFERENCE_LINES_HPP

// Line energies and K X-ray tables the CZT (czt/) and CdTe (cdte/) analyses
// share, in keV. One source per quantity:
//   decay gamma-ray energies: DDEP (LNHB), the nuclide's table;
//   capture gamma-ray energies: CapGam (NNDC);
//   X-ray energies: NIST SRD 128 (Deslattes et al. 2003), its direct
//   experimental value, else its theoretical one;
//   X-ray intensities: Table of Isotopes, 8th ed. (Firestone 1996), Table 7a
//   (per 100 K-shell vacancies),
//   https://grundpraktikum.physik.uni-saarland.de/gpalt/Anleitungen/Ergaenzungen/Versuche/R+oe+ntgenfluoreszenz/X_ray_Energies_Intensities.pdf
// SRD 128 lists no K-O lines; KO2,3 takes Table 7a's energy.

#include <RtypesCore.h>

namespace ReferenceLines {

// ------------------------------------------------------------ gamma lines

// http://www.lnhb.fr/nuclides/Am-241_tables.pdf: 59.5409(1)
const Double_t AM241 = 59.5409;
// http://www.lnhb.fr/nuclides/Ba-133_tables.pdf: 79.6142(19), 80.9979(11)
const Double_t BA133_79 = 79.6142;
const Double_t BA133_81 = 80.9979;
// 113Cd(n,g): the E2 transition in 114Cd whose K X-ray escapes fall among
// the lines near Ge. 95.902(3),
// https://www.nndc.bnl.gov/capgam/by-target.html?dsid=113CD%28N%2CG%29%20E%3DTHERMAL
const Double_t CD113_CAPTURE = 95.902;
// The measured line, ENSDF; for comparison only, never a fit input.
const Double_t GE73_ENSDF = 68.752;

// -------------------------------------------------------------- Pb K X-rays

// https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Pb&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV
// K-L2 72.80542(24), K-L3 74.97011(17).
const Double_t PB_KA2 = 72.80542;
const Double_t PB_KA1 = 74.97011;

// A line of an X-ray group: its energy and its intensity relative to the
// group's strongest line.
struct GroupLine {
  Double_t energy, weight;
};

// Pb K-beta groups: K-M3, K-M2, Kb5 (the mean of K-M4 and K-M5, which Table 7a
// gives one intensity), then K-N3, K-N2, Kb4 (K-N4,5), KO2,3. Intensities
// within each group from Table 7a, its Kb2 split between K-N2 and K-N3 as
// Kb3:Kb1.
const Int_t N_PB_KB1 = 3, N_PB_KB2 = 4;
const GroupLine PB_KB1[N_PB_KB1] = {
    {84.93908, 1.0}, {84.45045, 0.5215}, {85.4845, 0.0292}};
const GroupLine PB_KB2[N_PB_KB2] = {
    {87.3669, 1.0}, {87.238, 0.5215}, {87.589, 0.0350}, {87.911, 0.272}};

// ------------------------------------------- 95.902 keV line's K X-ray escapes

// The Cd and Te K X-rays whose escape leaves a line at the capture line's
// energy minus theirs: Cd Ka1, Ka2, Kb1', Kb2', Te Ka1, Ka2, Kb1', Kb2'. Kb1'
// is Kb1, Kb3 and Kb5 at their intensity-weighted energy, Kb2' Kb2, Kb4 and
// (Te) KO2,3; Kb5 and Kb4 at the mean of their two SRD 128 components.
// https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Cd&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV
// https://physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Te&trans=KL2&trans=KL3&trans=KM2&trans=KM3&trans=KM4&trans=KM5&trans=KN2&trans=KN3&trans=KN4&trans=KN5&lower=&upper=&units=eV
const Int_t N_ESCAPE_LINES = 8;
const Double_t ESCAPE_XRAY_KEV[N_ESCAPE_LINES] = {23.17398, 22.98405, 26.08491,
                                                  26.64430, 27.47257, 27.20199,
                                                  30.97997, 31.70765};
// Series of each line: 0 Cd K-alpha, 1 Cd K-beta, 2 Te K-alpha, 3 Te K-beta.
const Int_t ESCAPE_SERIES[N_ESCAPE_LINES] = {0, 0, 1, 1, 2, 2, 3, 3};
// Each line's share of its series is Geant4's, per detector
// (czt/GainTransferFit.hpp, cdte/CdTeAnalysis.hpp).

} // namespace ReferenceLines

#endif
