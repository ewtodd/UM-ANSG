// The germanium ionization-yield measurements that ModifiedLindhard.cpp fits,
// as each source gives them. Values are as printed, in the source's own unit
// (yield as a fraction or a percentage, or ionization in keV); the macro turns
// them into yields and splits their uncertainties. Where a source has its
// values only in a figure, the numbers are readings of that figure and say
// so. The sources were checked point by point against Kavner et al.'s
// Fig. 6 (arXiv source of the super-Lindhard paper); the papers and that
// comparison are in ~/Papers/Reading/73Ge-PRC/QF-sources.
#ifndef GE_YIELD_DATA_HPP
#define GE_YIELD_DATA_HPP

#include <Rtypes.h>

namespace GeYieldData {

// A point with its upper and lower errors, in the source's unit.
struct SourcePoint {
  Double_t recoil_kev;
  Double_t value;
  Double_t error_up;
  Double_t error_down;
};

// A point read off a figure: the marker and the ends of its error bar.
struct FigurePoint {
  Double_t recoil_kev;
  Double_t value;
  Double_t bar_high;
  Double_t bar_low;
};
// A bar end not drawn, or hidden behind another marker.
const Double_t NO_BAR = -1.0;

// ---- Ge(Li) detectors ----

// C. Chasman, K. W. Jones, R. A. Ristinen, Phys. Rev. Lett. 15, 245 (1965),
// doi:10.1103/PhysRevLett.15.245, Table I p.247: E_R (keV) and eta(E)/E_R.
// Ge(Li) at 77 K (Fig. 2). The per-point error is two fixed terms in
// quadrature (+/-0.7 keV on eta(E), +/-1.2 keV on E_R; footnotes b and d);
// nothing is stated as common to the points.
const Int_t CHASMAN1965_N = 15;
const SourcePoint CHASMAN1965[CHASMAN1965_N] = {
    {20.2, 0.347, 0.042, 0.042}, {26.9, 0.299, 0.029, 0.029},
    {33.1, 0.295, 0.024, 0.024}, {39.3, 0.313, 0.020, 0.020},
    {45.0, 0.298, 0.017, 0.017}, {50.7, 0.314, 0.016, 0.016},
    {56.3, 0.320, 0.014, 0.014}, {61.9, 0.331, 0.013, 0.013},
    {67.5, 0.332, 0.012, 0.012}, {73.1, 0.335, 0.011, 0.011},
    {78.5, 0.336, 0.010, 0.010}, {84.1, 0.341, 0.010, 0.010},
    {89.4, 0.344, 0.009, 0.009}, {94.9, 0.343, 0.009, 0.009},
    {100.3, 0.352, 0.008, 0.008}};

// A. R. Sattler, F. L. Vook, J. M. Palms, Phys. Rev. 143, 588 (1966),
// doi:10.1103/PhysRev.143.588, Table I p.593: energy of the maximum Ge
// recoil (keV) and its pulse height relative to an electron of the same
// energy. Lithium-drifted Ge at liquid-nitrogen temperature (p.590). "The
// vertical error includes estimated systematic error in interpreting
// individual spectra" (p.592); nothing is stated as common to the points.
const Int_t SATTLER1966_N = 18;
const SourcePoint SATTLER1966[SATTLER1966_N] = {
    {21.4, 0.149, 0.057, 0.057},  {26.8, 0.175, 0.027, 0.027},
    {32.2, 0.215, 0.022, 0.022},  {32.2, 0.224, 0.022, 0.022},
    {37.5, 0.223, 0.030, 0.030},  {42.9, 0.234, 0.027, 0.027},
    {48.2, 0.249, 0.023, 0.023},  {54.3, 0.235, 0.021, 0.021},
    {80.4, 0.281, 0.015, 0.015},  {107.0, 0.274, 0.030, 0.030},
    {161.0, 0.312, 0.031, 0.031}, {268.0, 0.422, 0.045, 0.045},
    {321.0, 0.438, 0.039, 0.039}, {342.0, 0.460, 0.044, 0.044},
    {765.0, 0.555, 0.045, 0.045}, {916.0, 0.692, 0.045, 0.045},
    {921.0, 0.692, 0.045, 0.045}, {997.0, 0.680, 0.045, 0.045}};

// C. Chasman, K. W. Jones, R. A. Ristinen, J. T. Sample, Phys. Rev. 154, 239
// (1967), doi:10.1103/PhysRev.154.239, Table I p.242: E (keV) and
// eta(E)/E for the Ge72 690 keV and Ge74 596 keV levels, neutron detector at
// 90 and 150 deg. Ge(Li) at 77 K (Fig. 6). Table II (p.242) splits the
// printed error: the level energy, +/-1.0 keV of the +/-1.2 keV on eta(E),
// is common to the two points of a level; the mean neutron energy, +/-20 keV
// or +/-1.1 keV of the 103.1 keV recoil, is common to all four.
const Int_t CHASMAN1967_N = 4;
const SourcePoint CHASMAN1967[CHASMAN1967_N] = {{56.4, 0.346, 0.022, 0.022},
                                                {103.1, 0.365, 0.012, 0.012},
                                                {56.1, 0.317, 0.022, 0.022},
                                                {102.7, 0.346, 0.012, 0.012}};
// Which level each point is: 0 for Ge72 (690 keV), 1 for Ge74 (596 keV).
const Int_t CHASMAN1967_LEVEL[CHASMAN1967_N] = {0, 0, 1, 1};
const Double_t CHASMAN1967_LEVEL_ERR_KEV = 1.0;
const Double_t CHASMAN1967_NEUTRON_RECOIL_ERR_KEV = 1.1;
const Double_t CHASMAN1967_NEUTRON_RECOIL_KEV = 103.1;

// C. Chasman, K. W. Jones, H. W. Kraner, W. Brandt, Phys. Rev. Lett. 21,
// 1430 (1968), doi:10.1103/PhysRevLett.21.1430, Table I p.1432: mean recoil
// energy and energy lost to ionization, both keV (the printed ratio column
// is these rounded to two decimals). Ge(Li). The subtracted 691.4 +/- 0.1
// keV level energy (footnote 13) is common to all six; taken as inside the
// printed errors, as the same group's 1967 budget has its level energy.
const Int_t CHASMAN1968_N = 6;
const SourcePoint CHASMAN1968[CHASMAN1968_N] = {
    {10.0, 2.0, 0.2, 0.2}, {17.6, 4.0, 0.5, 0.5}, {19.4, 4.7, 0.5, 0.5},
    {22.5, 5.7, 0.5, 0.5}, {24.1, 6.5, 0.5, 0.5}, {29.8, 8.2, 0.5, 0.5}};
const Double_t CHASMAN1968_LEVEL_ERR_KEV = 0.1;

// K. W. Jones, H. W. Kraner, Phys. Rev. C 4, 125 (1971),
// doi:10.1103/PhysRevC.4.125, Table II p.129: recoil energy and energy lost
// to ionization, both keV; +/-0.020 keV, +/-0.04 keV at 0.96 keV (footnote
// b). 1 cc Ge(Li). The ionization is the measured line minus the 68.752 keV
// level energy of their Table I, both on their calibration (241Am and Pb
// K-alpha). Recoil energies: "+/-0.2 keV at E_R = 1.75 keV and +/-0.1 keV
// at E_R = 1.04. The relative uncertainty is estimated to be +/-0.05 keV"
// (footnote a): a common scale and a point-to-point part.
const Int_t JONES1971_N = 8;
const SourcePoint JONES1971[JONES1971_N] = {
    {0.96, 0.169, 0.04, 0.04},   {1.04, 0.224, 0.020, 0.020},
    {1.12, 0.276, 0.020, 0.020}, {1.22, 0.272, 0.020, 0.020},
    {1.35, 0.349, 0.020, 0.020}, {1.48, 0.334, 0.020, 0.020},
    {1.62, 0.435, 0.020, 0.020}, {1.75, 0.326, 0.020, 0.020}};
const Double_t JONES1971_LEVEL_KEV = 68.752;
const Double_t JONES1971_SCALE_ANCHOR_KEV[2] = {1.04, 1.75};
const Double_t JONES1971_SCALE_ERR_KEV[2] = {0.1, 0.2};
const Double_t JONES1971_RELATIVE_RECOIL_ERR_KEV = 0.05;

// K. W. Jones, H. W. Kraner, Phys. Rev. A 11, 1347 (1975),
// doi:10.1103/PhysRevA.11.1347, Table II p.1351: the gamma + recoil line,
// "average weighted by 1/sigma^2", keV, on their calibration. Ge(Li).
const Double_t JONES1975_GAMMA_PLUS_RECOIL_KEV = 68.7927;
const Double_t JONES1975_GAMMA_PLUS_RECOIL_ERR_KEV = 0.0034;

// ---- HPGe detectors ----

// L. Baudis et al., Nucl. Instrum. Meth. A 418, 348 (1998),
// doi:10.1016/S0168-9002(98)00782-7, Fig. 3 ("this work"): yield as a
// fraction. Figure only; read from the vector paths of the arXiv version
// (hep-ex/9901028, p.11), which is the same plot as the journal's. Five
// points have no drawn bar: it is shorter than the marker radius,
// BAUDIS1998_MARKER_RADIUS. n-type coaxial HPGe; the composition of the
// uncertainties is not described.
const Int_t BAUDIS1998_N = 26;
const FigurePoint BAUDIS1998[BAUDIS1998_N] = {
    {59.436, 0.3192, 0.3390, 0.2994},  {68.429, 0.3211, 0.3409, 0.3012},
    {75.485, 0.3173, 0.3372, 0.2969},  {77.424, 0.3353, 0.3452, 0.3254},
    {79.717, 0.3260, 0.3464, 0.3062},  {83.774, 0.3334, 0.3533, 0.3130},
    {86.068, 0.3483, 0.3582, 0.3384},  {88.712, 0.3384, 0.3582, 0.3180},
    {93.651, 0.3409, 0.3564, 0.3260},  {94.180, 0.3502, 0.3595, 0.3409},
    {103.525, 0.3471, 0.3551, 0.3390}, {105.643, 0.3421, 0.3570, 0.3273},
    {111.994, 0.3570, 0.3725, 0.3421}, {114.285, 0.3502, 0.3582, 0.3421},
    {123.280, 0.3570, 0.3650, 0.3489}, {127.336, 0.3452, 0.3601, 0.3303},
    {130.157, 0.3613, 0.3694, 0.3533}, {136.331, 0.3669, 0.3824, 0.3520},
    {137.390, 0.3644, 0.3725, 0.3564}, {146.209, 0.3830, 0.3991, 0.3669},
    {147.268, 0.3799, NO_BAR, NO_BAR}, {157.142, 0.3824, 0.3991, 0.3650},
    {158.377, 0.3780, NO_BAR, NO_BAR}, {173.720, 0.3799, NO_BAR, NO_BAR},
    {185.009, 0.3892, NO_BAR, NO_BAR}, {197.178, 0.4053, NO_BAR, NO_BAR}};
const Double_t BAUDIS1998_MARKER_RADIUS = 0.0068;

// P. S. Barbeau, J. I. Collar, O. Tench, JCAP 09 (2007) 009,
// doi:10.1088/1475-7516/2007/09/009, Fig. 4 inset (point-contact HPGe, 24 keV
// iron-filtered neutrons); the same points are in CoGeNT, Phys. Rev. D 88,
// 012002 (2013), Fig. 5. No paper prints them: read here from J. I. Collar,
// A. R. L. Kavner, C. M. Lewis, Phys. Rev. D 103, 122003 (2021), Fig. 9 (gray
// squares), yield in percent; they agree with both earlier figures. The
// lowest point, 0.65 keV, is left out: "an indirectly derived value assigned
// to scatters anywhere in a broad recoil energy range" (Collar 2021). The
// upper bar of the 0.986 keV point is hidden behind another marker.
const Int_t COGENT_N = 3;
const FigurePoint COGENT[COGENT_N] = {{0.7380, 18.48, 21.19, 15.75},
                                      {0.9864, 19.14, NO_BAR, 16.53},
                                      {1.2254, 19.77, 21.65, 17.96}};

// J. I. Collar, A. R. L. Kavner, C. M. Lewis, Phys. Rev. D 103, 122003 (2021),
// doi:10.1103/PhysRevD.103.122003, Fig. 8 panel labels p.7 (iron-filtered 24
// keV neutrons): quenching factor in percent. HPGe at 77 K (Fig. 9). Each
// bar holds the fit error and the energy-calibration error; the calibration's
// intercept, 49 eV "with a small uncertainty of 11 eV", is common to the four.
const Int_t COLLAR2021_FEF_N = 4;
const SourcePoint COLLAR2021_FEF[COLLAR2021_FEF_N] = {{0.45, 35.3, 3.5, 4.4},
                                                      {0.59, 31.7, 2.3, 2.2},
                                                      {0.75, 30.0, 2.2, 2.2},
                                                      {0.96, 25.1, 1.5, 1.5}};
const Double_t COLLAR2021_INTERCEPT_ERR_KEV = 0.011;

// The same paper, Sec. IV p.6: the gamma + recoil line in the OSURR data,
// "68.811 +/- 0.001 keV".
const Double_t COLLAR2021_GAMMA_PLUS_RECOIL_KEV = 68.811;
const Double_t COLLAR2021_GAMMA_PLUS_RECOIL_ERR_KEV = 0.001;

// A. Bonhomme et al., Eur. Phys. J. C 82, 815 (2022),
// doi:10.1140/epjc/s10052-022-10768-1, Fig. 14: the quenching factor q as a
// fraction, the repeated angles of each beam energy merged. The merged points
// are printed nowhere (the per-run Tables 5-6 are not merged); read here from
// the vector figure. HPGe at 88 K. The bars are uncorrelated (statistics and
// angle); the correlated uncertainties (beam energy, energy scale, trigger
// efficiency) are the figure's band, whose half-width at each point is read
// alongside. Table 4: the energy scale, +/-12 eVee, is fully correlated, the
// beam energy correlated within each beam set.
struct BonhommePoint {
  Int_t beam_set;
  Double_t recoil_kev;
  Double_t q;
  Double_t bar_high;
  Double_t bar_low;
  Double_t band_half_width;
};
const Int_t BONHOMME2022_N = 31;
const BonhommePoint BONHOMME2022[BONHOMME2022_N] = {
    {1, 0.60751, 0.16029, 0.17557, 0.14500, 0.0291},
    {2, 0.80022, 0.16535, 0.18076, 0.14995, 0.0177},
    {1, 0.80880, 0.17593, 0.19130, 0.16055, 0.0174},
    {0, 0.96614, 0.18053, 0.19388, 0.16719, 0.0135},
    {3, 0.97963, 0.17218, 0.19134, 0.15303, 0.0133},
    {2, 1.06325, 0.17504, 0.18889, 0.16118, 0.0120},
    {0, 1.16646, 0.18711, 0.19935, 0.17487, 0.0107},
    {1, 1.23341, 0.17053, 0.18077, 0.16030, 0.0101},
    {3, 1.30399, 0.18471, 0.19890, 0.17051, 0.0095},
    {0, 1.54146, 0.18613, 0.19461, 0.17765, 0.0079},
    {0, 1.58068, 0.19189, 0.19954, 0.18424, 0.0077},
    {2, 1.62167, 0.16789, 0.17836, 0.15743, 0.0075},
    {1, 1.72922, 0.18871, 0.19559, 0.18184, 0.0070},
    {1, 1.90071, 0.18559, 0.19415, 0.17704, 0.0064},
    {0, 1.92481, 0.18602, 0.19661, 0.17544, 0.0063},
    {3, 1.99017, 0.17207, 0.18244, 0.16171, 0.0061},
    {0, 1.99047, 0.19035, 0.19487, 0.18582, 0.0061},
    {1, 2.33124, 0.19984, 0.20663, 0.19306, 0.0052},
    {2, 2.49256, 0.18674, 0.19334, 0.18014, 0.0048},
    {3, 2.78920, 0.19465, 0.20175, 0.18754, 0.0043},
    {2, 3.01799, 0.19501, 0.20380, 0.18622, 0.0040},
    {1, 3.07239, 0.20155, 0.20584, 0.19726, 0.0039},
    {1, 3.78434, 0.21372, 0.21703, 0.21041, 0.0032},
    {3, 3.78921, 0.22079, 0.22980, 0.21179, 0.0032},
    {1, 3.91842, 0.21462, 0.21797, 0.21127, 0.0031},
    {2, 4.03596, 0.20510, 0.20969, 0.20050, 0.0030},
    {3, 4.95341, 0.21770, 0.22246, 0.21293, 0.0024},
    {2, 4.97229, 0.21853, 0.22217, 0.21489, 0.0024},
    {2, 5.14794, 0.21522, 0.21859, 0.21186, 0.0023},
    {3, 6.09577, 0.22554, 0.22944, 0.22164, 0.0020},
    {3, 6.31714, 0.23456, 0.23847, 0.23064, 0.0019}};
// Table 1: the adopted mean neutron energy of each beam set, keV (runs 3-4,
// 7-9, 5-6 and 10-11).
const Int_t BONHOMME2022_N_BEAMS = 4;
const Double_t BONHOMME2022_BEAM_KEV[BONHOMME2022_N_BEAMS] = {248.5, 490.3,
                                                              644.0, 790.0};
const Double_t BONHOMME2022_BEAM_ERR_KEV[BONHOMME2022_N_BEAMS] = {1.0, 3.0, 4.0,
                                                                  11.0};
const Double_t BONHOMME2022_ENERGY_SCALE_ERR_KEV = 0.012;

// A. R. L. Kavner, I. Jovanovic, Phys. Rev. D 110, 083043 (2024),
// doi:10.1103/PhysRevD.110.083043, Sec. IV D p.7: the gamma + recoil line
// with cuts, calibration choice and a nuisance term, "68.8158 +/- 0.0046
// keV". HPGe at 77 K.
const Double_t KAVNER2024_GAMMA_PLUS_RECOIL_KEV = 68.8158;
const Double_t KAVNER2024_GAMMA_PLUS_RECOIL_ERR_KEV = 0.0046;

// ---- The capture-line points as published, for showing the data before
// this work: each paper's ionization against its own gamma energy, and its
// recoil energy, keV ----

// Jones and Kraner 1975, Table II p.1351 ("Recoil ionization energy loss
// 39.2 +/- 5.5 eV") and p.1349 (254.1 eV recoil).
const Double_t JONES1975_PUBLISHED_IONIZATION_KEV = 0.0392;
const Double_t JONES1975_PUBLISHED_IONIZATION_ERR_KEV = 0.0055;
const Double_t JONES1975_PUBLISHED_RECOIL_KEV = 0.2541;
// Collar 2021, Sec. IV ("77 +/- 20 eV of ionization being produced by a
// 0.254 keVnr germanium NR").
const Double_t COLLAR2021_PUBLISHED_IONIZATION_KEV = 0.077;
const Double_t COLLAR2021_PUBLISHED_IONIZATION_ERR_KEV = 0.020;
const Double_t COLLAR2021_PUBLISHED_RECOIL_KEV = 0.254;
// Kavner 2024, Sec. IV D ("a nuclear recoil ionization of 63.8 +/- 8.4
// eVee") at 254.1 eVnr.
const Double_t KAVNER2024_PUBLISHED_IONIZATION_KEV = 0.0638;
const Double_t KAVNER2024_PUBLISHED_IONIZATION_ERR_KEV = 0.0084;
const Double_t KAVNER2024_PUBLISHED_RECOIL_KEV = 0.2541;

} // namespace GeYieldData

#endif
