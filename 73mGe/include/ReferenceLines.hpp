#ifndef REFERENCE_LINES_HPP
#define REFERENCE_LINES_HPP

// Line energies and K X-ray tables the CZT (czt/) and CdTe (cdte/) analyses
// share, in keV.

#include <RtypesCore.h>

namespace ReferenceLines {

// ------------------------------------------------------------ gamma lines

const Double_t AM241 = 59.5409;
const Double_t BA133_79 = 79.6142;
const Double_t BA133_81 = 80.9979;
// 113Cd(n,g): the E2 transition in 114Cd whose K X-ray escapes fall among
// the lines near Ge.
const Double_t CD113_CAPTURE = 95.902;
// The measured line, ENSDF; for comparison only, never a fit input.
const Double_t GE73M_ENSDF = 68.752;

// -------------------------------------------------------------- Pb K X-rays

const Double_t PB_KA2 = 72.8042;
const Double_t PB_KA1 = 74.9694;

// A line of an X-ray group: its energy and its intensity relative to the
// group's strongest line.
struct GroupLine {
  Double_t energy, weight;
};

// Pb K-beta groups: NIST SRD 128 experimental energies (K-M3, K-M2, K-N3,
// K-N2, K-N4,5; Kb5 and KO2,3 from the Table of Isotopes), intensities within
// each group from the Table of Isotopes (its Kb2 split between K-N2 and K-N3
// as Kb3:Kb1).
const Int_t N_PB_KB1 = 3, N_PB_KB2 = 4;
const GroupLine PB_KB1[N_PB_KB1] = {
    {84.93908, 1.0}, {84.45045, 0.5215}, {85.470, 0.0292}};
const GroupLine PB_KB2[N_PB_KB2] = {
    {87.3669, 1.0}, {87.238, 0.5215}, {87.589, 0.0350}, {87.911, 0.272}};

// ------------------------------------------- 95.902 keV line's K X-ray escapes

// The Cd and Te K X-rays whose escape leaves a line at the capture line's
// energy minus theirs: Cd Ka1, Ka2, Kb1', Kb2', Te Ka1, Ka2, Kb1', Kb2'; Kb1'
// and Kb2' at their groups' intensity-weighted energies.
const Int_t N_ESCAPE_LINES = 8;
const Double_t ESCAPE_XRAY_KEV[N_ESCAPE_LINES] = {
    23.1739, 22.9843, 26.092, 26.650, 27.4726, 27.202, 30.988, 31.720};
// Series of each line: 0 Cd K-alpha, 1 Cd K-beta, 2 Te K-alpha, 3 Te K-beta.
const Int_t ESCAPE_SERIES[N_ESCAPE_LINES] = {0, 0, 1, 1, 2, 2, 3, 3};
// Each line's share of its series: DDEP emission intensities times the
// probability of escaping through the face of a slab.
const Double_t ESCAPE_SERIES_SHARE[N_ESCAPE_LINES] = {
    1.0 / 1.5242, 0.5242 / 1.5242, 0.3394 / 0.4035, 0.0641 / 0.4035,
    1.0 / 1.5275, 0.5275 / 1.5275, 0.3548 / 0.4349, 0.0801 / 0.4349};

} // namespace ReferenceLines

#endif
