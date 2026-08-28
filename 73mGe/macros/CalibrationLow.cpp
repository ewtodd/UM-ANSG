#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <RtypesCore.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

// calibration
const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_79 = 79.6142;
const Float_t E_BA133_81 = 80.9979;
// IUPAC convention: Ka1 is the STRONGER, HIGHER-energy line. These were
// previously named the other way round -- the energies were paired with the
// right measured centroids, so no result was wrong, but "Ka1" referred to the
// weak line, which is confusing in the code and wrong in paper text and figure
// captions. Note the intensity ratio (~1.7:1 in favour of Ka1) is why the
// Ka2 centroid carries the larger fit error (~27 eV vs ~4 eV), and Ka2 is the
// low-side anchor nearest the Ge line, so its precision matters most.
// Ordering elsewhere is by ascending energy, i.e. Ka2 first, then Ka1.
const Float_t E_PB_KA2 = 72.8042;
const Float_t E_PB_KA1 = 74.9694;
const Float_t E_PB_KB1 = 84.936;

// just used for seeding fits
const Float_t E_PB_KB2 = 87.32;
const Float_t E_GE_73M = 68.752;
const Float_t E_BKG_73KEV = 73.5;

// Anchor reference-energy uncertainties (keV), folded into the cal-fit x-errors
// so they propagate into the pol1 parameter covariance and thus the Ge energy.
// Pb Ka are known to <1 eV (Deslattes); Kb1 carries the ~2:1-vs-1.9:1 intensity
// -blend ambiguity (~6 eV).
//
// Am-241 is split out from the generic source-line bucket because under
// USE_AM_TRANSFER it is a calibration ANCHOR rather than a cross-check line, so
// its reference uncertainty sets how hard it pulls the fit. The evaluated line
// is 59.5409(1) keV -> 0.1 eV; the old lumped 10 eV would have down-weighted it
// 100x against Pb Ka and largely defeated the point of anchoring on it.
//
// dE_LINE remains the conservative ~10 eV lump, now used ONLY for the Ba-133
// lines of the 15th/16th line-cal days.
// TODO: replace with per-line evaluated Ba-133 uncertainties.
const Float_t dE_PB_KA = 0.001;
const Float_t dE_PB_KB1 = 0.006;
const Float_t dE_AM241 = 0.0001;
const Float_t dE_LINE = 0.010;

const Float_t CAL_RANGE_LOW = 61.5;
// High edge pulled in to 79.5 (was 80.5): the extra ~1 keV of continuum toward
// the Pb-Kb complex (~84.8) gave the Ge/Pb high-exp tails high-side structure
// to chase, railing the high-exp RATIO at its bound (~100 = near-infinite
// decay). Trimming it lets the high-exp tail settle at a finite ratio without
// removing the component (which is physical -- it absorbs the Pb-Ka2
// spillover). Still >4 keV above Pb-Ka1 (74.97), so the doublet + its tails
// stay in-window.
const Float_t CAL_RANGE_HIGH = 80.5;

// Lineshape model for the WHOLE chain (precal -> cal -> postcal): shared
// doublet tail on a flat continuum -- a free linear term is degenerate with the
// high-side pileup tail and unstable on these spectra. RESULT_TAG namespaces
// every output (calibrated trees, plots, fits, .result).
//
// CORRECTION: earlier comments here described the background/continuum
// systematic as a fit-range variation. It is not, and never was -- no fit range
// is varied anywhere in this macro (CAL_RANGE_LOW/HIGH are fixed). ge_bkg_err
// is |cal(precal Ge mu) - postcal Ge mu|: the disagreement between the
// calibration evaluated at the raw Ge centroid and what the post-cal fit
// actually returned. The rationale for treating that as a background term is
// argued at its computation site in ProcessPbAnchoredPair; note that it also
// absorbs seeding and shape-lock differences, so it is a conservative bound on
// the background alone.
const Bool_t TAIL_SHARED = kTRUE;

// Calibration scheme switch. This is the ONLY difference between the published
// analysis and the Am-anchored one, so any shift in the Ge energy is
// attributable to it alone.
//   kFALSE -> published scheme: per-run Pb-only self-cal (Ka1, Ka2, Kb1); Am
//             used solely for the VALIDATION and LINEARITY cross-checks.
//   kTRUE  -> the scheme the manuscript describes: the pair closest in time to
//             the Am run fits an Am+Pb reference line, and every other pair
//             transfers that line by its measured Pb-Ka gain (no refit).
const Bool_t USE_AM_TRANSFER = kTRUE;

// Polynomial degree of the Am-anchored reference calibration.
//
// A pol1 through Am + the three Pb lines gives chi2/ndf ~ 49 with a 165 eV miss
// at Kb1, and a pol2 fits all four to <4 eV. But that comparison is NOT clean
// evidence of detector nonlinearity: dE_AM241 pins Am at 0.1 eV against Pb's
// 1 eV, so the chi2 is dominated by Am declining to sit on the Pb line. The
// three Pb-only local slopes span just 0.7% (0.9975, 1.0046) while the Am->Ka1
// step is 1.4% away from them -- i.e. the apparent curvature is concentrated
// entirely in the one anchor that is weighted hardest. A biased Am centroid and
// a genuinely curved response are indistinguishable from this fit alone.
// Kept switchable so both can be run and compared.
// Include Pb-Kb1 (84.936) as a calibration anchor and as the third point of the
// gain transfer. It was adopted because it improved the calibration term from
// 17.8 to 12.6 eV -- but that was measured while Kb1 reached the weighted fit
// with a ZERO error bar (see GuardedMuError), so it was not fitting better, it
// was outvoting the other anchors. With honest errors that gain has to be
// re-earned. Kb1's saved states also predate TAIL_RATIO_MAX and carry railed
// tails, and Am + Ka1 + Ka2 without it are collinear at chi2/ndf 0.15 while
// any set containing Kb1 sits at 3.0-4.1.
const Bool_t USE_KB1_ANCHOR = kTRUE;
const Int_t REF_CAL_DEGREE = 2;

// The post-cal simultaneous fit can lock sigma / amplitude / tail state to the
// precal values (RunBkgGeSim's lock_postcal_state). With it kTRUE the lineshape
// is pinned, so the Ge centroid cannot respond to a change of calibration at
// all, and the reported fit error carries no shape uncertainty. That was a
// leftover test setting rather than a modelling choice: every other fit in this
// macro fixes only the structurally-disabled components (step amplitudes,
// background slope), which is 5 of 28 parameters, not 25 of 28.
// Lock the post-cal peak SHAPE (sigma + tails) to the precal values, leaving
// centroids and yields free. kTRUE is the conditioned default: the shape is
// determined on raw data where it is well constrained, and the post-cal fit
// then only has to place the peaks on the calibrated scale. Releasing the shape
// as well (kFALSE) puts ~23 free parameters on 4-14 M unbinned events, several
// of which sit railed at their bounds, and Minuit fails to build a Hessian.
// The free-shape variant is the natural lineshape systematic, run separately.
const Bool_t LOCK_POSTCAL_STATE = kFALSE;

// The post-cal simultaneous fit must NOT run in interactive mode.
//
// RooFitUtils' simultaneous path has three branches. Interactive with a saved
// .simroofits present REPLAYS the stored parameters and never fits at all
// (sim_valid is set straight from the load). Interactive with no saved file
// opens the GUI editor and waits for a human, then saves. Only the
// non-interactive branch actually minimises.
//
// Branch one is a trap here: after the calibration changes, a re-run looks
// entirely normal -- it prints chi2, regenerates every plot and writes a fresh
// .result -- while silently replaying the Ge centroid from the last hand-fit
// against the OLD calibration. The rising chi2 is the only symptom. The
// pre-cal fits can stay interactive (they act on raw data that does not change
// between calibration schemes); the post-cal fit is exactly the stage that must
// respond, so it is forced non-interactive and genuinely re-minimised.
const Bool_t POSTCAL_INTERACTIVE = kTRUE;

// Tie the Pb Ka doublet to a single Gaussian width in the simultaneous fits.
// See the LinkParameter call in RunBkgGeSim for why.
const Bool_t LINK_PB_KA_SIGMA = kTRUE;

// Lock the post-cal Pb centroids to their KNOWN energies (72.8042 / 74.9694)
// and let only the Ge float.
//
// Normally the Pb peaks are left free post-cal precisely so that where they
// land is a blind check. That check is currently failing in two runs: the
// doublet splays, with Ka2 pulled low and Ka1 pushed high, giving a separation
// 26-47 eV wider than the 2.1652 keV that atomic physics fixes it at.
//
// Locking them inverts the diagnostic. If the Ge centroid barely moves, the Pb
// splay is a local lineshape problem that does not propagate to the line of
// interest. If the Ge moves substantially, then whatever distorts the Pb
// doublet is also displacing the Ge, and the free-Pb numbers are contaminated
// by it. Either answer is worth having; the Ge is left free in both cases.
// Constrain the Pb Ka doublet SPACING to its tabulated value in the post-cal
// fit. This is deliberately not the same statement as LOCK_POSTCAL_PB, which
// pins both centroids to truth and so fixes the absolute scale as well: this
// leaves the pair free to slide together -- absorbing any residual gain or
// offset error, which is the calibration's job to own -- and forbids only the
// STRETCH.
//
// The stretch is the defect we can actually point at. Fitting the transfer on
// the Ka pair alone (two points, one gain parameter) gives chi2/ndf 4.3-4.6,
// i.e. Ka1 and Ka2 disagree about the gain by ~2 sigma, or ~7e-4, which is
// ~48 eV at the Ge energy -- the size of the Ge scatter. The post-cal doublet
// separation wanders over 13-62 eV run to run against a spacing tabulated to
// under 1 eV. The suspected mechanism is that the weak Ka2 sits on the strong
// Ka1's low-energy tail and the two share tail parameters, so an unconstrained
// tail is free to trade against the separation.
// Constrain the Pb Ka doublet SPACING in the PRE-cal fit to its tabulated
// value, scaled by the local gain. Default OFF, and NOT because it is
// unfinished -- because there is nothing to correct.
//
// The fitted spacing runs 14-23 eV wider than tabulated across the four pairs.
// That looked like a common-mode lineshape bias for most of one session, but
// the comparison being made was against the LITERATURE uncertainty on the
// spacing (<1 eV), which is the wrong denominator. The relevant uncertainty is
// this measurement's: the Pb centroids carry ~10.7 and ~6.0 eV, so the fitted
// separation is good to ~12 eV, and 14-23 eV is 1.2-1.9 sigma. Ordinary
// scatter. The shared sign across runs is not evidence either -- they share a
// lineshape model and a detector, so they are not independent draws.
//
// Kept because the diagnostic print is worth having and the machinery is
// correct, but do not treat the excess as a defect to chase.
//
// If it is ever switched on: the raw axis carries an unknown gain, so 2.1652
// cannot be imposed directly. The doublet's SUM fixes the local gain and the
// DIFFERENCE is then constrained to 2.1652 * g; sum and difference are
// near-independent, so the constraint does not feed on what it constrains.
// Needs a second pass, since g is unknown until the first fit has run.
const Bool_t CONSTRAIN_PRECAL_SEPARATION = kFALSE;
const Bool_t CONSTRAIN_PB_KA_SEPARATION = kFALSE;
// Uncertainty on the constraint, keV. The Deslattes Ka energies are each known
// to well under 1 eV; 1.5 eV on the difference is conservative. This stays a
// measurement rather than a hard equality -- a zero would have to fix a mu.
const Double_t dE_PB_KA_SEP = 0.0015;
const Bool_t LOCK_POSTCAL_PB = kFALSE;

// High-side exponential tail. Justified in an earlier comment as a pileup
// shoulder above Pb-Ka1, but it never takes an interior value: across six peak
// instances in the three hand-fitted runs it is either switched off (amplitude
// 1.7e-08, 6.6e-11, 2.5e-07, 2.2e-06) or railed at BOTH bounds at once
// (amplitude 0.49999999, ratio 99.9995). A 100 keV decay across a 19 keV window
// is not a tail, it is a flat pedestal -- and therefore degenerate with
// BkgConstant, which is why the fit slides it to a bound. The run where it
// rails is also the run with the worst Pb doublet distortion (+46.7 eV against
// +1.9 eV where it is off). The Ge peak also sits 4 keV BELOW the doublet, so a
// high-side pileup component on it has little physical motivation.
const Bool_t USE_HIGH_EXP_TAIL = kFALSE;

// Flat vs linear continuum. The flat-vs-linear energy difference is the
// preferred background systematic, but it was previously abandoned as unstable
// because "a free linear term is degenerate with the high-side tail" -- with
// USE_HIGH_EXP_TAIL off that degeneracy is removed, so the comparison should
// become meaningful again. Run both and take the difference.
const Bool_t USE_FLAT_BKG = kTRUE;

// Upper bound on the exponential tail decay lengths, in keV.
//
// The library default is 100, which exceeds the ~19 keV fit window: over that
// range the exponential is flat, so the "tail" is a constant pedestal and is
// degenerate with the background yield. The fits duly park it at the bound --
// HighExpTailRatio pinned at 100 with amplitude pinned at 0.5, repeatedly,
// hand-tuned fits included, and LowExpTailRatio wandering 1.4 to 88 across
// identical runs.
//
// A charge-collection or pileup shoulder decays over a few keV. Capping at 8
// forces the component to describe a real tail or fit to nothing, instead of
// absorbing continuum. This is the leading candidate for shrinking the 43 eV
// lineshape systematic, which exists only because the tail-on and tail-off
// variants cannot currently be told apart on fit quality.
// Fit the calibration peaks WITH a step (error-function shelf) as well as the
// tails. Incomplete charge collection produces a genuine step below a peak;
// with no step term the fit has to build one out of the tail components, which
// inflates both their amplitude and their decay length.
//
// Am-241 is the case that matters. Fit without a step it returns LowExpTail
// amplitude 0.4845 AND LowLinTail amplitude 0.3838 -- both large, which is the
// signature of two tails jointly standing in for a shelf they cannot represent
// -- and a decay of tau = 1.129 keV that we then took as the detector's true
// tail. If that tau is an artifact of the missing step, then imposing it on the
// Pb doublet over-corrects, which is what the cap-1.5 test looked like: the
// separations improved but did not close, and the two methods disagreed MORE.
const Bool_t USE_STEP_CAL_PEAKS = kFALSE;
const Double_t TAIL_RATIO_MAX = 8.0;

// Transfer model: affine (offset + gain) vs pure gain. Both now see the same
// three Pb lines (Ka2, Ka1, Kb1), so switching this isolates the MODEL from the
// choice of points -- adding a free parameter cannot worsen chi2 on fixed data,
// so any degradation is attributable to the points, not the model.
const Bool_t USE_AFFINE_TRANSFER = kFALSE;

// Rate-subtraction tail, FIXED (see FitRateSubtractedGe). Nominal is no tail;
// set FRAC to the in-situ lineshape's tail fraction to measure the systematic.
const Double_t RATESUB_TAIL_FRAC = 0.0;
const Double_t RATESUB_TAIL_TAU = 1.0;
// Diagnostic tail scan on the 01/14 rate-sub (largest dataset of the campaign,
// so the only one that can test the fixed-zero tail the others assume).
// Prints only; does not feed the result rows.
// Measure the 13th's pairs by rate subtraction as well, so every dataset in the
// campaign is measured BOTH ways and the six can also be combined under a
// single method. With one method across six datasets the in-situ/rate-sub
// difference stops being a systematic (12.2 eV, the second-largest term) and
// becomes a consistency check instead.
//
// Each pair is subtracted using ITS OWN transferred calibration, so the in-situ
// and rate-sub rows for a given run differ only in how the peak was fit --
// nothing else changes between them.
//
// Not obviously a win: rate subtraction discards the background run's counts
// and roughly doubles the variance, so its fit errors run 25-52 eV against
// 5-28 eV in-situ. The gain is cleanliness, and whether that beats the method
// systematic is what this measures.
const Bool_t RATESUB_ALL_DAYS = kFALSE;
const Bool_t RATESUB_TAIL_SCAN = kFALSE;

const TString RESULT_TAG = USE_AM_TRANSFER ? "amxfer" : "shared";

// Lineshape variant suffix for the .result FILENAME ONLY.
//
// Deliberately not folded into RESULT_TAG: that also drives OUT_SUBDIR and
// therefore where the .roofits/.simroofits live, so changing it would orphan
// every saved fit state (including hand-tuned ones) and silently send the fits
// back to cold starts. Only the results file is namespaced.
//
// The two variants are the lineshape systematic: the high-side exponential tail
// on versus off. They bracket the literature value (68.7335 vs 68.7768 against
// 68.752) with comparable fit quality, so neither can be rejected and the
// half-spread is the honest uncertainty. This restores the two-file lineshape
// mechanism CombineGeResult was built around and lost when CalibrationLowKa1
// was deleted.
const TString RESULT_VARIANT = USE_HIGH_EXP_TAIL ? "_hitail" : "_nohitail";
const TString OUT_SUBDIR = "calibrated_low_" + RESULT_TAG;

// When non-empty, RunBkgGeSim dumps the converged bkg+sig channel fits to
// <G_CSV_PREFIX>_{bkg,sig}.csv (data + smooth fit curve + residual pulls) for
// external plotting, then returns. Set by ExportFitCSV(); empty in normal runs.
TString G_CSV_PREFIX = "";

struct CalibrationData {
  std::vector<Float_t> mu, mu_errors, calibration_values_keV, reduced_chi2;
  // Reference-energy (y-axis) uncertainty per anchor; 0 -> use a tiny default.
  std::vector<Float_t> energy_errors;
  std::vector<TString> run_names;
};

struct BkgGeSimResult {
  FitResult bkg_channel;
  FitResult sig_channel;
  Bool_t valid = kFALSE;
};

struct DayResult {
  std::vector<TString> cal_labels;
  std::vector<TF1 *> cal_funcs;
  std::vector<TString> ge_labels;
  std::vector<Float_t> ge_mus, ge_errs, ge_chi2s;
  // Calibration-term energy uncertainty at the Ge peak, parallel to ge_* (the
  // pol1 covariance propagated to E(mu_Ge); the ge_errs are the Ge fit term).
  // ONE calibration term. The old orthogonal slope/offset pivot split was added
  // at a collaborator's request; it is not defined for the pol2 calibration the
  // Am anchor requires (three eigen-directions, not two), and a single total is
  // what we quote.
  std::vector<Float_t> ge_cal_errs;
  // The TRANSFER part of ge_cal_errs, parallel to ge_*. ge_cal_err is
  // hypot(reference_cal_err, gain_err), so the reference part -- COMMON to
  // every run on the day -- is sqrt(cal^2 - gain^2) and the transfer part is
  // independent per run. The combiner needs that split to build the covariance
  // correctly instead of assuming the whole calibration is correlated or none
  // of it is.
  std::vector<Float_t> ge_gain_errs;
  // |cal(precal Ge) - postcal Ge| per run, parallel to ge_* (NOT a fit-range
  // variation -- see the note at the top of this file).
  std::vector<Float_t> ge_bkg_errs;
  // Method 2 (rate-subtraction) Ge results, parallel to the ge_* vectors.
  std::vector<TString> rs_labels;
  std::vector<Float_t> rs_mus, rs_errs, rs_chi2s, rs_cal_errs, rs_gain_errs;
};

// Calibration line/curve plus its parameter covariance. npar is 2 for pol1 and
// 3 for pol2; the p2 terms stay zero for pol1 so every formula below reduces to
// the linear case exactly.
struct CalFit {
  TF1 *func = nullptr;
  Double_t p0 = 0, p1 = 1;
  Int_t npar = 2;
  Double_t var_p0 = 0, var_p1 = 0, cov_p0p1 = 0;
  Double_t var_p2 = 0, cov_p0p2 = 0, cov_p1p2 = 0;
};

struct PairCalResult {
  TF1 *cal_func = nullptr;
  Double_t ge_cal_err = 0; // calibration sigma propagated to the Ge energy

  Double_t ge_bkg_err = 0; // |cal(precal Ge) - postcal Ge|; see file-top note
  BkgGeSimResult postcal;
  TString pair_tag;
  std::vector<Float_t> pb_mus_precal;
  std::vector<Float_t> pb_mu_errs_precal;
  // Populated only on the Am reference pair: the master line that the day's
  // other pairs gain-transfer off. ge_gain_err is the transfer term already
  // folded into ge_cal_err, kept separately for diagnostics.
  CalFit ref_calfit;
  Bool_t has_ref_cal = kFALSE;
  Double_t ge_gain_err = 0;
};

struct LineCalConfig {
  TString date_label;
  TString am_run;
  Float_t am_lo, am_hi;
  TString ba_run;
  std::vector<TString> day_datasets;
  TString postcal_bkg;
  TString postcal_sig;
  TString ge_label;
  // Rate-subtraction fit window. Per-day because the residual quality differs:
  // the 01/15 subtraction is flat and clean across the full 64-73 keV span,
  // while the 01/16 residual is over-subtracted below ~66 keV (running to -5
  // counts/s) and unsettled above ~71.5 -- curvature a linear background term
  // cannot absorb, which biases the centroid.
  Float_t rs_lo = 62.0;
  Float_t rs_hi = 71.5;
  // Allow a sloped residual background per day. Default flat: see the note in
  // FitRateSubtractedGe about the slope competing with the tail.
  // Linear by default: the data prefers it. Forcing the 01/15 residual flat
  // costs real fit quality (chi2/ndf 1.36 -> 1.97) because that residual
  // genuinely slopes, and it moves the centroid by 35 eV (68.7002 -> 68.7348).
  // That flat-vs-linear difference is the rate-sub's dominant background
  // systematic, larger than anything else in this method's budget.
  //
  // The tail contributes nothing either way: tailFrac railed at 0 under BOTH
  // background models, with tau at its lower bound. After subtraction the
  // effective statistics are S/sqrt(S+B) with a large removed B, and a
  // low-amplitude broad tail is the first feature to become unresolvable. This
  // measurement therefore cannot constrain a lineshape and must not be used to
  // argue about one.
  Bool_t rs_linear_bkg = kTRUE;
};

std::vector<Double_t> LoadFiltered(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    return events;
  }
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    std::vector<Double_t> ch =
        RooFitUtils::LoadEventsFromTree(tree, "energykeV");
    events.insert(events.end(), ch.begin(), ch.end());
  }
  file->Close();
  delete file;
  return events;
}

std::vector<Double_t> LoadCalibrated(const TString &input_name) {
  std::vector<Double_t> events;
  TFile *file = IO::OpenForReading(OUT_SUBDIR + "/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << OUT_SUBDIR << "/" << input_name
              << ".root" << std::endl;
    return events;
  }
  TTree *tree = static_cast<TTree *>(file->Get("calibrated_tree"));
  if (tree)
    events = RooFitUtils::LoadEventsFromTree(tree, "depositedEnergykeV");
  file->Close();
  delete file;
  return events;
}

FitResult FitCalPeak(const std::vector<Double_t> &events,
                     const TString &input_name, const TString &peak_name,
                     Float_t fit_low, Float_t fit_high, Bool_t double_peak,
                     Double_t mu1, Double_t mu2, Bool_t interactive,
                     Bool_t link_sigma = kFALSE, Bool_t use_step = kFALSE) {
  if (events.empty())
    return {};
  RooFitUtils fitter(events, fit_low, fit_high, Constants::BIN_WIDTH_KEV, kTRUE,
                     use_step, kTRUE, kTRUE, kTRUE);
  fitter.SetTailRatioMax(TAIL_RATIO_MAX);
  if (interactive)
    fitter.SetInteractive();
  if (double_peak)
    return fitter.FitDoublePeak(input_name, peak_name, mu1, mu2, link_sigma);
  return fitter.FitSinglePeak(input_name, peak_name);
}

BkgGeSimResult RunBkgGeSim(
    const std::vector<Double_t> &bkg_events,
    const std::vector<Double_t> &sig_events,
    const std::vector<Double_t> &bkg_peak_mus, const TString &bkg_label,
    const TString &sig_label, const TString &fit_label, Bool_t interactive,
    Bool_t lock_bkg_peaks = kFALSE,
    const FitResult *precomputed_bkg_seed = nullptr,
    const FitResult *precomputed_sig_seed = nullptr,
    // lock_postcal_state now means ONLY "lock the peak SHAPE after
    // seeding". It used to be passed to all three of AddChannel's
    // bkg_yield_fixed / bkg_slope_fixed / lock_shape_after_seed slots,
    // which meant clearing it also freed the background SLOPE on a
    // model configured with a FLAT background -- a degenerate parameter
    // that helps wreck the Hessian. Yield is always free (normalisation
    // must adapt to the calibrated spectrum); slope is always fixed.
    Bool_t lock_postcal_state = kFALSE, Bool_t share_tail = kTRUE,
    Float_t fit_lo = CAL_RANGE_LOW, Float_t fit_hi = CAL_RANGE_HIGH,
    // Force the SIMULTANEOUS fit to actually minimise, independently of
    // the background seed fit. In interactive mode RooFitUtils replays a
    // saved .simroofits and never fits, so the Ge centroid is a stored
    // constant. The seed fit may legitimately keep replaying: it acts on
    // raw filtered data that does not change between calibrations.
    Bool_t force_sim_fit = kFALSE,
    // kTRUE when the x axis is already CALIBRATED keV (the post-cal fit), so
    // tabulated line spacings apply to it directly. Pre-cal the axis is raw and
    // carries the unknown gain, so no external spacing can be imposed there.
    // Separation constraint for the Pb Ka doublet, in the units of THIS fit's
    // x axis. Negative disables. Post-cal that axis is keV so the tabulated
    // 2.1652 applies directly; pre-cal it is raw, so the caller must scale by
    // the local gain first (see CONSTRAIN_PRECAL_SEPARATION).
    Double_t constrain_sep_delta = -1.0, Double_t constrain_sep_sigma = -1.0) {
  BkgGeSimResult out;
  out.valid = kFALSE;
  Int_t n_bkg = (Int_t)bkg_peak_mus.size();
  if (n_bkg < 1 || n_bkg > 2 || bkg_events.empty() || sig_events.empty())
    return out;

  const Bool_t kFlatBkg = USE_FLAT_BKG;
  const Bool_t kStep = kFALSE;
  const Bool_t kLowExp = kTRUE;
  // LOW-LINEAR TAIL DISABLED.
  //
  // LowExp and LowLin are two low-side components acting over the same few
  // hundred eV, and they trade against each other almost freely. Evidence:
  //  - AmLineshapeScan: with BOTH enabled the automated fit reaches chi2/ndf
  //    740; with either one alone, 1.67 (LowExp) or 1.99 (LowLin).
  //  - The fitted LowExpTailRatio takes 14.7 / 88.0 / 1.4 / 20.7 across four
  //    runs of the same line -- a decay length varying by a factor of 60 --
  //    and HighExpTailRatio rails at its 100 limit in two of them.
  //
  // The Ge tail is linked to bkg:Ka1's tail via share_tail, so it inherits that
  // instability directly. The Pb line has enormous statistics and still cannot
  // pin the shape, because the degeneracy is structural rather than
  // statistical. That instability is what makes the Ge centroid swing ~57 eV
  // between the shape-locked and shape-free post-cal fits.
  const Bool_t kLowLin = kTRUE;
  const Bool_t kHighExp = USE_HIGH_EXP_TAIL;

  FitResult seed;
  if (precomputed_bkg_seed) {
    seed = *precomputed_bkg_seed;
  } else {
    RooFitUtils bkg_fitter(bkg_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                           kFlatBkg, kStep, kLowExp, kLowLin, kHighExp);
    bkg_fitter.SetTailRatioMax(TAIL_RATIO_MAX);
    if (interactive)
      bkg_fitter.SetInteractive();
    TString seed_label = fit_label + "_BkgSeed";
    if (n_bkg == 1)
      seed = bkg_fitter.FitSinglePeak(bkg_label, seed_label);
    else
      seed = bkg_fitter.FitDoublePeak(bkg_label, seed_label, bkg_peak_mus[0],
                                      bkg_peak_mus[1]);
    if (!seed.valid) {
      std::cerr << "ERROR: bkg-only seed fit failed for " << bkg_label
                << std::endl;
      return out;
    }
  }

  RooFitUtils sim;
  sim.SetTailRatioMax(TAIL_RATIO_MAX);
  if (interactive) {
    sim.SetInteractive();
    // Saved parameters SEED a real minimisation instead of being adopted as the
    // answer. Keeps the hand-tuned starting point -- which the automated cold
    // start cannot reach (chi2/ndf ~1.4 seeded vs ~8 cold, edm 24000x
    // tolerance) -- while making the result reproducible and able to respond to
    // a changed calibration.
    if (force_sim_fit)
      sim.SetRefitAfterLoad();
  }

  std::vector<Double_t> sig_mus = bkg_peak_mus;
  sig_mus.push_back((Double_t)E_GE_73M);
  std::vector<Bool_t> bkg_fixed(n_bkg, lock_bkg_peaks);
  std::vector<Bool_t> sig_fixed(n_bkg + 1, kFALSE);
  for (Int_t i = 0; i < n_bkg; i++)
    sig_fixed[i] = lock_bkg_peaks;

  // Step shelf on the Pb/bkg peaks, off on the Ge peak (last sig peak).
  std::vector<Bool_t> sig_step(n_bkg + 1, kTRUE);
  sig_step[n_bkg] = kFALSE;

  sim.AddChannel("bkg", bkg_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg, bkg_peak_mus, kFlatBkg, kStep, kLowExp, kLowLin,
                 kHighExp, bkg_fixed, kFALSE, kTRUE, lock_postcal_state);
  sim.AddChannel("sig", sig_events, fit_lo, fit_hi, Constants::BIN_WIDTH_KEV,
                 n_bkg + 1, sig_mus, kFlatBkg, kStep, kLowExp, kLowLin,
                 kHighExp, sig_fixed, kFALSE, kTRUE, lock_postcal_state,
                 sig_step);
  // Share ONE low-side tail shape (LowExp + LowLin) across every Pb-Ka instance
  // AND the Ge peak, all tied to bkg:Ka1. Fit independently, Ka1 sits on Ka2's
  // tail and inflates its own (LowExp ratio railed ~5.3) to absorb the inter-
  // peak fill, so neither single Ka tail is trustworthy. Collapsing the doublet
  // to one tail recovers the value the combined lines support -- a clean ~73
  // keV tail. The Ge then borrows THAT instead of the inflated Ka1-only tail
  // (which over-corrected its centroid ~+0.05 keV high). Without the shared
  // tail the Ge tail otherwise rails to zero -> bare-Gaussian centroid biased
  // low, worse at low statistics. Ge mu/sigma/yield stay free; only the tail
  // SHAPE is borrowed.
  //
  // These are pushed BEFORE LinkPeakShape so they win first-match for the four
  // tail params; LinkPeakShape still ties mu/sigma/step/high-exp per line.
  // NOTE: routing every tail target straight to bkg:Ka1 (not chaining via
  // bkg:Ka2) avoids the ResolveOrCreate quirk where a link-resolved var is not
  // registered under its own key, which would orphan sig:Ka2's tail into a free
  // parameter. share_tail selects the lineshape model: kTRUE = one tail shared
  // across the Pb doublet + Ge (de-inflated); kFALSE = Ge tail tied to Ka1
  // only. Running both per pair and taking the Ge-mu spread gives the lineshape
  // systematic.
  // Only link components that actually exist: with kLowLin disabled the
  // LowLinTail* parameters are never built, and linking to a non-existent
  // source aborts the channel build.
  const char *tp[4] = {"LowExpTailAmplitude", "LowExpTailRatio",
                       "LowLinTailAmplitude", "LowLinTailSlope"};
  const Int_t n_tp = kLowLin ? 4 : 2;
  TString ge_idx = TString::Format("%d", n_bkg + 1);
  if (share_tail) {
    for (Int_t t = 0; t < n_tp; t++) {
      if (n_bkg == 2)
        sim.LinkParameter(TString("bkg:") + tp[t] + "2",
                          TString("bkg:") + tp[t] + "1");
      for (Int_t i = 1; i < n_bkg; i++)
        sim.LinkParameter(TString("sig:") + tp[t] +
                              TString::Format("%d", i + 1),
                          TString("bkg:") + tp[t] + "1");
      sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                        TString("bkg:") + tp[t] + "1");
    }
  } else {
    // Alternative lineshape: Ge tail tied to Ka1 only (no doublet share). The
    // Pb Ka tails fit independently; used to bound the lineshape systematic.
    for (Int_t t = 0; t < n_tp; t++)
      sim.LinkParameter(TString("sig:") + tp[t] + ge_idx,
                        TString("bkg:") + tp[t] + "1");
  }

  // Tie the Pb Ka doublet to ONE detector resolution. The two lines are 2.17
  // keV apart and share a resolution by construction, but their sigmas floated
  // independently, which is the same condition that wrecked the Ba-133 doublet
  // (separation compressed 39%, line cal chi2/ndf 213). Here it shows up as the
  // WEAK line -- peak 1, the 72.804 keV Ka2 -- inflating its width and trading
  // against the continuum: its centroid error is ~27 eV against Ka1's ~4 eV.
  // That matters because Ka2 is the low-side anchor nearest the Ge line and it
  // dominates the calibration-independent R_GEPB uncertainty (~80 eV).
  //
  // Direction is forced: LinkParameter requires the SOURCE to be built already,
  // and parameters are created in peak order, so only "2 follows 1" is legal
  // (asking for Sigma1 <- Sigma2 errors and then segfaults). That costs nothing
  // physically -- a linked sigma is a single shared free parameter constrained
  // by BOTH peaks' data, and the strong line dominates that constraint either
  // way; only which one is nominally free changes.
  if (constrain_sep_delta > 0 && n_bkg == 2)
    sim.ConstrainPeakSeparation("bkg", 1, 0, constrain_sep_delta,
                                constrain_sep_sigma > 0 ? constrain_sep_sigma
                                                        : dE_PB_KA_SEP);

  if (LINK_PB_KA_SIGMA && n_bkg == 2)
    sim.LinkParameter("bkg:Sigma2", "bkg:Sigma1");

  for (Int_t i = 0; i < n_bkg; i++)
    sim.LinkPeakShape("sig", i, "bkg", i);

  sim.SeedChannel("bkg", seed);
  if (precomputed_sig_seed)
    sim.SeedChannel("sig", *precomputed_sig_seed);

  std::vector<FitResult> results = sim.FitSimultaneous(sig_label, fit_label);
  if (results.size() < 2)
    return out;
  out.bkg_channel = results[0];
  out.sig_channel = results[1];
  out.valid = results[1].valid;

  // CSV export hook: dump the converged channels for external plotting. Each
  // call writes <prefix>_<chan>_spectrum.csv (data+fit+residual per bin) and
  // <prefix>_<chan>_fitcurve.csv (smooth overlay line).
  if (!G_CSV_PREFIX.IsNull()) {
    sim.DumpChannelCSV("bkg", G_CSV_PREFIX + "_bkg");
    sim.DumpChannelCSV("sig", G_CSV_PREFIX + "_sig");
  }
  return out;
}

void AddCalPoint(CalibrationData &cal_data, const TString &name, Float_t mu,
                 Float_t mu_err, Float_t true_e, Float_t chi2,
                 Float_t energy_err = 0) {
  cal_data.run_names.push_back(name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_err);
  cal_data.calibration_values_keV.push_back(true_e);
  cal_data.reduced_chi2.push_back(chi2);
  cal_data.energy_errors.push_back(energy_err);
}

// Cramer-Rao floor on a fitted centroid. A Gaussian of width sigma holding N
// counts cannot locate its centre better than sigma/sqrt(N); a fit reporting
// less than that has a broken covariance, not a better measurement.
//
// Pb-Kb1 was reaching the calibration with a centroid error under 0.05 eV --
// 260x smaller than the STRONGER Ka2 line at 13.1 eV, on a weaker peak. Its
// saved lineshape state predates TAIL_RATIO_MAX (LowExpTailRatio 30-78 against
// a cap of 8, several tail amplitudes railed at the 0.5 bound), and parameters
// sitting on their bounds collapse Hesse. Both AddCalPoint and the gain
// transfer weight by 1/mu_error^2, so that one railed anchor was effectively
// defining the calibration line while every other point contributed pure chi2.
// That is the transfer chi2/ndf of 3.5-4.1 -- against chi2/ndf 0.15 for Am plus
// the Ka pair on their own, which are collinear to well within their errors.
//
// This is a floor, not a correction. It cannot rescue a centroid that is in the
// wrong place; it only stops a broken error bar from outvoting good anchors. A
// peak that trips this warning still wants its saved state re-tuned under the
// cap -- the floor makes the fit honest, not right.
Float_t GuardedMuError(const TString &name, const PeakFitResult &pk,
                       const FitResult &fr) {
  if (!(pk.sigma > 0) || !(pk.gaus_amplitude > 0))
    return pk.mu_error;
  const Float_t floor_err = pk.sigma / TMath::Sqrt(pk.gaus_amplitude);
  if (pk.mu_error > floor_err)
    return pk.mu_error;
  std::cerr << "WARNING: " << name << " centroid error " << pk.mu_error * 1000
            << " eV is below the Poisson floor " << floor_err * 1000
            << " eV (sigma " << pk.sigma << ", N " << pk.gaus_amplitude
            << ", covQual " << fr.cov_qual << ") -- raising to the floor. "
            << "Re-tune this fit's saved state under TAIL_RATIO_MAX."
            << std::endl;
  return floor_err;
}

void PrintCalSummary(const CalibrationData &cal_data,
                     const TString &date_label) {
  std::cout << "Calibration points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); i++) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0)
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
    std::cout << std::endl;
  }
}

// Calibration fit result carrying the pol1 covariance so the Ge energy
// uncertainty from calibration can be propagated analytically:
//   E(mu) = p0 + p1*mu
//   var(E) = var(p0) + mu^2 var(p1) + 2 mu cov(p0,p1)   (cal term)
//          + (p1 * sigma_mu_Ge)^2                       (Ge fit term, added
//          later)

// Calibration energy uncertainty at a given precalibrated mu (cal term only).
//   var(E) = var_p0 + mu^2 var_p1 + 2 mu cov(p0,p1)
// Full quadratic-form propagation: var(E) = sum_ij mu^(i+j) cov(pi,pj).
// For pol1 the p2 terms are zero and this is exactly
// var_p0 + mu^2 var_p1 + 2 mu cov_p0p1, as before.
Double_t CalEnergyError(const CalFit &c, Double_t mu) {
  Double_t var = c.var_p0 + mu * mu * c.var_p1 + 2.0 * mu * c.cov_p0p1;
  if (c.npar >= 3) {
    Double_t mu2 = mu * mu;
    var += mu2 * mu2 * c.var_p2 + 2.0 * mu2 * c.cov_p0p2 +
           2.0 * mu2 * mu * c.cov_p1p2;
  }
  return (var > 0) ? std::sqrt(var) : 0.0;
}

// Pivot decomposition of the calibration error into ORTHOGONAL slope and offset
// terms at the energy of interest, so they may be quoted separately AND summed
// in quadrature without dropping the (p0,p1) covariance. Re-expand the line
// about the pivot x0 = -cov/var_p1, where slope and intercept are uncorrelated:
//   offset = sqrt(var_p0 - cov^2/var_p1)        (intercept error at the pivot)
//   slope  = |mu - x0| * sqrt(var_p1)           (slope error x lever arm to mu)
// By construction offset^2 + slope^2 == CalEnergyError(c, mu)^2 exactly.

// degree 1 -> pol1 (the per-run Pb self-cal). degree 2 -> pol2, required once
// Am-241 joins the Pb lines: the CZT response is measurably nonlinear over
// 59.5-85 keV (local slope runs 0.9841 -> 0.9975 -> 1.0046 between consecutive
// anchors), so a straight line cannot pass through Am and the three Pb lines at
// once -- it misses Kb1 by ~165 eV. A quadratic fits all four to <4 eV.
CalFit CreateAndSavePol1Cal(const CalibrationData &cal_data,
                            const TString &date_label, Bool_t fix_p0 = kFALSE,
                            Float_t p0_value = 0, Int_t degree = 1) {
  Int_t n = (Int_t)cal_data.mu.size();
  // Errors in BOTH axes: x = precal mu fit error, y = reference-energy error.
  std::vector<Float_t> ex(cal_data.mu_errors);
  std::vector<Float_t> ey(n);
  for (Int_t i = 0; i < n; i++) {
    Float_t e = (i < (Int_t)cal_data.energy_errors.size())
                    ? cal_data.energy_errors[i]
                    : 0.0f;
    ey[i] = (e > 0)
                ? e
                : 1e-4f; // tiny floor so the point isn't infinitely weighted
    if (!(ex[i] > 0))
      ex[i] = 1e-4f;
  }
  TGraphErrors *graph = new TGraphErrors(n, cal_data.mu.data(),
                                         cal_data.calibration_values_keV.data(),
                                         ex.data(), ey.data());
  // Draw range spans all anchor points (the Pb K-beta at ~85 keV sits above
  // CAL_RANGE_HIGH; Am at ~59.5 below CAL_RANGE_LOW) so the fit line covers
  // them.
  Double_t draw_lo = CAL_RANGE_LOW, draw_hi = CAL_RANGE_HIGH;
  for (Int_t i = 0; i < n; i++) {
    draw_lo = std::min({draw_lo, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
    draw_hi = std::max({draw_hi, (Double_t)cal_data.mu[i],
                        (Double_t)cal_data.calibration_values_keV[i]});
  }
  draw_lo -= 3.0;
  draw_hi += 3.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      graph, kBlue, "; Precalibrated Energy [keV]; Deposited Energy [keV]");
  graph->GetXaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetYaxis()->SetRangeUser(draw_lo, draw_hi);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetMarkerStyle(5);
  graph->SetMarkerSize(2);
  graph->Draw("AP");

  TF1 *cal = new TF1("cal_low_" + date_label, (degree >= 2) ? "pol2" : "pol1",
                     draw_lo, draw_hi);
  cal->SetParameter(0, 0);
  cal->SetParameter(1, 1);
  if (degree >= 2)
    cal->SetParameter(2, 0);
  // Drift datasets hold the reference offset and only let the slope (gain)
  // float so the measured Pb doublet lands on truth.
  if (fix_p0)
    cal->FixParameter(0, p0_value);
  cal->SetNpx(1000);
  // "S" returns the fit result (covariance); "EX0"? no -- we WANT x-errors
  // used, so plain "E S". ROOT's linear-fitter handles x-errors via effective
  // variance.
  TFitResultPtr fr = graph->Fit(cal, "E S Q");
  cal->Draw("SAME");

  CalFit out;
  out.func = cal;
  out.npar = cal->GetNpar();
  out.p0 = cal->GetParameter(0);
  out.p1 = cal->GetParameter(1);
  if (fr.Get()) {
    out.var_p0 = fr->CovMatrix(0, 0);
    out.var_p1 = fr->CovMatrix(1, 1);
    out.cov_p0p1 = fr->CovMatrix(0, 1);
    if (out.npar >= 3) {
      out.var_p2 = fr->CovMatrix(2, 2);
      out.cov_p0p2 = fr->CovMatrix(0, 2);
      out.cov_p1p2 = fr->CovMatrix(1, 2);
    }
  } else {
    out.var_p0 = cal->GetParError(0) * cal->GetParError(0);
    out.var_p1 = cal->GetParError(1) * cal->GetParError(1);
    if (out.npar >= 3)
      out.var_p2 = cal->GetParError(2) * cal->GetParError(2);
  }
  // Fit quality of the calibration itself. With Am + 3 Pb anchors a pol1 gives
  // chi2/ndf ~ 49 and a pol2 ~ 0: the nonlinearity is not subtle, so print it.
  if (fr.Get() && fr->Ndf() > 0)
    std::cout << "CAL-FIT " << date_label << " (pol" << (out.npar - 1) << ", "
              << n << " anchors): chi2/ndf = " << std::fixed
              << std::setprecision(2) << fr->Chi2() / fr->Ndf() << std::endl;

  PlottingUtils::SaveFigure(canvas, "calibration_low_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return out;
}

// AFFINE transfer: mu_run = a + b * mu_ref.
//
// The pure-gain model (mu_run = g * mu_ref) assumes the drift between runs is
// purely multiplicative. It is not: fitted against the two Pb Ka lines it gives
// chi2/ndf of 3.2-4.6 on ONE degree of freedom, i.e. the two lines disagree
// about the gain at 3-4.5 sigma. That disagreement is what the sqrt(chi2/ndf)
// inflation was converting into a 12-19 eV transfer uncertainty -- an honest
// penalty for a model that cannot describe the data.
//
// Allowing an offset as well as a gain gives the drift a shape it can actually
// take. Two Ka lines alone would determine a and b exactly (zero dof, no error
// estimate), so Pb Kb1 is included as a third point: 3 points, 2 parameters,
// 1 dof -- enough to fit the model AND report a residual uncertainty.
struct AffineResult {
  Double_t a = 0.0; // offset [keV]
  Double_t b = 1.0; // gain
  Double_t var_a = 0, var_b = 0, cov_ab = 0;
  Double_t chi2_ndf = 0;
  Int_t n = 0;
  Bool_t valid = kFALSE;
};

// Weighted straight-line fit of this run's Pb centroids against the reference
// run's, with errors on both axes.
AffineResult PbAffineVsReference(const std::vector<Float_t> &mu_run,
                                 const std::vector<Float_t> &err_run,
                                 const std::vector<Float_t> &mu_ref,
                                 const std::vector<Float_t> &err_ref) {
  AffineResult out;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  std::vector<Double_t> x, y, ex, ey;
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    x.push_back(mu_ref[i]);
    y.push_back(mu_run[i]);
    ex.push_back(i < err_ref.size() && err_ref[i] > 0 ? err_ref[i] : 1e-4);
    ey.push_back(i < err_run.size() && err_run[i] > 0 ? err_run[i] : 1e-4);
  }
  out.n = (Int_t)x.size();
  if (out.n < 3)
    return out; // need a spare point for a meaningful residual

  TGraphErrors g((Int_t)x.size(), x.data(), y.data(), ex.data(), ey.data());
  TF1 lin("pb_affine", "pol1", x.front() - 5.0, x.back() + 5.0);
  lin.SetParameters(0.0, 1.0);
  TFitResultPtr fr = g.Fit(&lin, "E S Q N");
  if (!fr.Get())
    return out;

  out.a = lin.GetParameter(0);
  out.b = lin.GetParameter(1);
  out.var_a = fr->CovMatrix(0, 0);
  out.var_b = fr->CovMatrix(1, 1);
  out.cov_ab = fr->CovMatrix(0, 1);
  out.chi2_ndf = (fr->Ndf() > 0) ? fr->Chi2() / fr->Ndf() : 0.0;
  // Same honesty as before: if the affine model still does not describe the
  // points, inflate rather than pretend.
  if (out.chi2_ndf > 1.0) {
    out.var_a *= out.chi2_ndf;
    out.var_b *= out.chi2_ndf;
    out.cov_ab *= out.chi2_ndf;
  }
  out.valid = kTRUE;
  return out;
}

// Transfer the reference calibration through the affine map. With
// mu_run = a + b*mu_ref, the reference scale is read at mu_ref = (mu - a)/b, so
// E_run(mu) = f_ref((mu - a)/b). For f_ref = p0 + p1 x + p2 x^2 that is again a
// polynomial in mu, with coefficients:
//   q0 = p0 - p1 a/b + p2 a^2/b^2
//   q1 = p1/b - 2 p2 a/b^2
//   q2 = p2/b^2
CalFit AffineScaleCal(const CalFit &ref, const AffineResult &t,
                      const TString &label) {
  TF1 *src = ref.func;
  Int_t npar = src->GetNpar();
  Double_t lo = 0, hi = 0;
  src->GetRange(lo, hi);
  TF1 *scaled = new TF1("cal_low_" + label, src->GetExpFormula(), lo, hi);

  Double_t p0 = src->GetParameter(0);
  Double_t p1 = (npar > 1) ? src->GetParameter(1) : 0.0;
  Double_t p2 = (npar > 2) ? src->GetParameter(2) : 0.0;
  Double_t a = t.a, b = (t.b != 0 ? t.b : 1.0);

  scaled->SetParameter(0, p0 - p1 * a / b + p2 * a * a / (b * b));
  if (npar > 1)
    scaled->SetParameter(1, p1 / b - 2.0 * p2 * a / (b * b));
  if (npar > 2)
    scaled->SetParameter(2, p2 / (b * b));
  scaled->SetNpx(1000);

  CalFit out;
  out.func = scaled;
  out.npar = npar;
  out.p0 = scaled->GetParameter(0);
  out.p1 = scaled->GetParameter(1);
  // Reference covariance carried through, as for the pure-gain case.
  out.var_p0 = ref.var_p0;
  out.var_p1 = ref.var_p1 / (b * b);
  out.cov_p0p1 = ref.cov_p0p1 / b;
  if (npar >= 3) {
    Double_t b2 = b * b;
    out.var_p2 = ref.var_p2 / (b2 * b2);
    out.cov_p0p2 = ref.cov_p0p2 / b2;
    out.cov_p1p2 = ref.cov_p1p2 / (b2 * b);
  }
  return out;
}

// Uncertainty from the transfer itself, propagated to the Ge energy.
// E = f_ref(z) with z = (mu - a)/b, so dE/da = -f'(z)/b and dE/db = -f'(z) z/b.
Double_t AffineTransferEnergyError(const CalFit &ref_cal, Double_t mu_raw,
                                   const AffineResult &t) {
  if (!ref_cal.func || !t.valid || t.b == 0)
    return 0.0;
  Double_t z = (mu_raw - t.a) / t.b;
  Double_t fp = ref_cal.func->Derivative(z);
  Double_t dda = -fp / t.b;
  Double_t ddb = -fp * z / t.b;
  Double_t var =
      dda * dda * t.var_a + ddb * ddb * t.var_b + 2.0 * dda * ddb * t.cov_ab;
  return (var > 0) ? std::sqrt(var) : 0.0;
}

struct GainResult {
  Double_t g = 1.0; // error-weighted mean relative gain
  Double_t sigma_g = 0.0;
  Double_t chi2_ndf = 0.0;
  Int_t n = 0;
};

// Relative gain of a run's Pb-Ka doublet against the reference run's, under the
// multiplicative-drift model: g = <mu_run / mu_ref>_w. The Pb energies are
// physically identical in both runs, so the ratio is pure gain. Weights combine
// both runs' centroid errors in quadrature on the ratio. When Ka1 and Ka2
// disagree by more than their errors allow (chi2/ndf > 1), sigma_g is inflated
// by sqrt(chi2/ndf): the drift was not purely multiplicative and the transfer
// deserves to be penalised for it.
GainResult PbGainVsReferenceWithError(const std::vector<Float_t> &mu_run,
                                      const std::vector<Float_t> &err_run,
                                      const std::vector<Float_t> &mu_ref,
                                      const std::vector<Float_t> &err_ref) {
  GainResult out;
  Double_t sw = 0, swg = 0;
  std::vector<Double_t> gs, ws;
  size_t n = std::min(mu_run.size(), mu_ref.size());
  for (size_t i = 0; i < n; i++) {
    if (mu_run[i] <= 0 || mu_ref[i] <= 0)
      continue;
    Double_t g = (Double_t)mu_run[i] / (Double_t)mu_ref[i];
    Double_t fr = (err_run[i] > 0) ? err_run[i] / mu_run[i] : 0.0;
    Double_t ff = (err_ref[i] > 0) ? err_ref[i] / mu_ref[i] : 0.0;
    Double_t rel = std::sqrt(fr * fr + ff * ff);
    Double_t sg = (rel > 0) ? g * rel : 1e-9;
    Double_t w = 1.0 / (sg * sg);
    sw += w;
    swg += w * g;
    gs.push_back(g);
    ws.push_back(w);
  }
  if (sw <= 0)
    return out;
  out.g = swg / sw;
  out.sigma_g = std::sqrt(1.0 / sw);
  out.n = (Int_t)gs.size();
  if (out.n > 1) {
    Double_t chi2 = 0;
    for (size_t i = 0; i < gs.size(); i++)
      chi2 += ws[i] * (gs[i] - out.g) * (gs[i] - out.g);
    out.chi2_ndf = chi2 / (out.n - 1);
    if (out.chi2_ndf > 1.0)
      out.sigma_g *= std::sqrt(out.chi2_ndf);
  }
  return out;
}

// Transfer a reference calibration to a drift run by pure gain rescaling:
// E_run(mu) = f_ref(mu / g). For f_ref = sum_k p_k mu^k this is exactly
// sum_k (p_k / g^k) mu^k, so coefficient k scales by g^-k. No refit -- the
// reference line is carried over intact and only its scale moves. Deliberately
// kept pol1: ge_bkg_err is defined as |cal(precal Ge) - postcal Ge|, which is
// only interpretable as a background term while the calibration stays linear.
CalFit GainScaleCal(const CalFit &ref, Double_t g, const TString &label) {
  TF1 *src = ref.func;
  Int_t npar = src->GetNpar();
  Double_t lo = 0, hi = 0;
  src->GetRange(lo, hi);
  TF1 *scaled = new TF1("cal_low_" + label, src->GetExpFormula(), lo, hi);
  Double_t gk = 1.0;
  for (Int_t k = 0; k < npar; k++) {
    scaled->SetParameter(k, src->GetParameter(k) / gk);
    gk *= g;
  }
  scaled->SetNpx(1000);
  CalFit out;
  out.func = scaled;
  out.p0 = scaled->GetParameter(0);
  out.p1 = scaled->GetParameter(1);
  // Reference covariance carried through the rescaling. Coefficient k scales by
  // g^-k, so cov(pi,pj) -> cov(pi,pj) / g^(i+j).
  out.npar = npar;
  out.var_p0 = ref.var_p0;
  out.var_p1 = ref.var_p1 / (g * g);
  out.cov_p0p1 = ref.cov_p0p1 / g;
  if (npar >= 3) {
    Double_t g2 = g * g;
    out.var_p2 = ref.var_p2 / (g2 * g2);
    out.cov_p0p2 = ref.cov_p0p2 / g2;
    out.cov_p1p2 = ref.cov_p1p2 / (g2 * g);
  }
  if (npar > 0)
    scaled->SetParError(0, std::sqrt(std::max(0.0, out.var_p0)));
  if (npar > 1)
    scaled->SetParError(1, std::sqrt(std::max(0.0, out.var_p1)));
  return out;
}

// Gain-transfer uncertainty propagated to the Ge energy. For E = f_ref(mu/g),
// dE/dg = -(mu/g^2) f_ref'(mu/g).
Double_t GainTransferEnergyError(const CalFit &ref_cal, Double_t mu_raw,
                                 Double_t g, Double_t sigma_g) {
  if (!(sigma_g > 0) || !(g > 0) || !ref_cal.func)
    return 0.0;
  Double_t fprime = ref_cal.func->Derivative(mu_raw / g);
  return std::fabs(-mu_raw / (g * g) * fprime) * sigma_g;
}

// Run start time (BEF header timeStart) for a dataset, used to pick the Pb
// dataset closest in time to the Am-241 reference run.
Float_t LoadRunStartTime(const TString &input_name) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name
              << ".root for run start time" << std::endl;
    if (file) {
      file->Close();
      delete file;
    }
    return 0;
  }
  Float_t t = 0;
  TTree *hdr = static_cast<TTree *>(file->Get("bef_header_tree"));
  if (hdr && hdr->GetEntries() > 0) {
    UInt_t timeStart = 0;
    hdr->SetBranchAddress("timeStart", &timeStart);
    hdr->GetEntry(0);
    t = (Float_t)timeStart;
  }
  file->Close();
  delete file;
  return t;
}

void ApplyPol1Cal(const std::vector<TString> &input_names, TF1 *cal,
                  const TString &date_label) {
  TFile *cal_file = IO::OpenForWriting(
      OUT_SUBDIR + "/calibration_function_low_" + date_label + ".root");
  cal->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  for (size_t i = 0; i < input_names.size(); i++) {
    TString input_name = input_names[i];
    TFile *in_file = IO::OpenForReading("filtered/" + input_name + ".root");
    if (!in_file || in_file->IsZombie()) {
      std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
                << std::endl;
      delete in_file;
      continue;
    }
    TFile *out_file =
        IO::OpenForWriting(OUT_SUBDIR + "/" + input_name + ".root");

    Float_t deposited_energy = 0;
    TTree *cal_tree = new TTree("calibrated_tree", "Calibrated events");
    cal_tree->SetDirectory(out_file);
    cal_tree->Branch("depositedEnergykeV", &deposited_energy,
                     "depositedEnergykeV/F");

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TTree *tree = static_cast<TTree *>(
          in_file->Get(Form("crystal%d_filtered_tree", c)));
      if (!tree)
        continue;
      Float_t energy = 0;
      tree->SetBranchAddress("energykeV", &energy);
      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        deposited_energy = cal->Eval(energy);
        cal_tree->Fill();
      }
    }

    out_file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    std::cout << "Calibrated " << input_name << " -> " << OUT_SUBDIR << "/"
              << input_name << ".root" << std::endl;
    out_file->Close();
    delete out_file;
    in_file->Close();
    delete in_file;
  }
}

void AddGeResult(DayResult &result, const BkgGeSimResult &r,
                 const TString &label, Float_t cal_err = 0,
                 Float_t gain_err = 0, Float_t bkg_err = 0) {
  if (!r.valid || r.sig_channel.peaks.empty())
    return;
  const PeakFitResult &ge = r.sig_channel.peaks.back();
  result.ge_labels.push_back(label);
  result.ge_mus.push_back(ge.mu);
  result.ge_errs.push_back(ge.mu_error);
  result.ge_chi2s.push_back(r.sig_channel.reduced_chi2);
  result.ge_cal_errs.push_back(cal_err);
  result.ge_gain_errs.push_back(gain_err);
  result.ge_bkg_errs.push_back(bkg_err);
}

// ---------------------------------------------------------------------------
// Method 2: live-time rate subtraction (15th/16th line-cal days).
//
// An independent cross-check of the in-situ simultaneous fit. The signal and
// background runs are calibrated PER CRYSTAL straight from the filtered file,
// each crystal's histogram normalized by THAT crystal's filtered live time
// (counts/s), then signal - background per crystal, summed over crystals. The
// Pb/continuum that the in-situ method fits through cancels in the difference,
// leaving the Ge peak on a near-flat residual -- a different background
// treatment, hence a genuinely independent systematic.
// ---------------------------------------------------------------------------

Float_t LoadCrystalLiveTime(TFile *file, Int_t crystal) {
  TParameter<Float_t> *p = static_cast<TParameter<Float_t> *>(
      file->Get(Form("LiveTime_Filtered_Crystal%d_s", crystal)));
  if (!p) {
    std::cerr << "WARNING: missing LiveTime_Filtered_Crystal" << crystal << "_s"
              << std::endl;
    return 0;
  }
  return p->GetVal();
}

// Per-crystal calibrated rate histogram (counts / live-time-seconds), summed
// over crystals, for one run. cal maps filtered energykeV -> deposited keV.
TH1F *BuildRateHist(const TString &input_name, TF1 *cal, const TString &hname,
                    Float_t lo, Float_t hi) {
  TFile *file = IO::OpenForReading("filtered/" + input_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << input_name << ".root"
              << std::endl;
    delete file;
    return nullptr;
  }
  Int_t nbins = (Int_t)std::lround((hi - lo) / Constants::BIN_WIDTH_KEV);
  TH1F *sum = new TH1F(hname, ";Deposited Energy [keV]; Rate [counts/s]", nbins,
                       lo, hi);
  sum->SetDirectory(nullptr);
  sum->Sumw2();

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TTree *tree =
        static_cast<TTree *>(file->Get(Form("crystal%d_filtered_tree", c)));
    if (!tree)
      continue;
    Float_t lt = LoadCrystalLiveTime(file, c);
    if (lt <= 0)
      continue;
    // Per-crystal counts histogram, then scale by 1/live-time so Sumw2 carries
    // the Poisson error of the raw counts through the rate normalization.
    TH1F *ch = new TH1F(hname + Form("_c%d", c), "", nbins, lo, hi);
    ch->SetDirectory(nullptr);
    ch->Sumw2();
    Float_t energy = 0;
    tree->SetBranchAddress("energykeV", &energy);
    Int_t n = (Int_t)tree->GetEntries();
    for (Int_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      ch->Fill(cal->Eval(energy));
    }
    ch->Scale(1.0 / lt);
    sum->Add(ch);
    delete ch;
  }
  file->Close();
  delete file;
  return sum;
}

// Fit the rate-subtracted (signal - background) residual Ge peak with the
// binned FittingUtils backend (a subtracted histogram has no event list for
// the unbinned RooFit path). Flat residual background absorbs imperfect
// cancellation. Returns the Ge peak mu/error via the result struct.
struct RateSubResult {
  Float_t mu = -1, mu_err = -1, chi2 = -1;
  Bool_t valid = kFALSE;
};

RateSubResult
FitRateSubtractedGe(const TString &sig_run, const TString &bkg_run, TF1 *cal,
                    const TString &tag, Bool_t interactive, Float_t lo = 62.0,
                    Float_t hi = 71.5, Bool_t linear_bkg = kFALSE,
                    Bool_t free_tail = kFALSE, Double_t tail_frac = -1.0,
                    Double_t tail_tau = -1.0, Double_t bkg_scale = 1.0) {
  RateSubResult out;

  TH1F *sig = BuildRateHist(sig_run, cal, "rate_sig_" + tag, lo, hi);
  TH1F *bkg = BuildRateHist(bkg_run, cal, "rate_bkg_" + tag, lo, hi);
  if (!sig || !bkg) {
    delete sig;
    delete bkg;
    return out;
  }
  TH1F *res = static_cast<TH1F *>(sig->Clone("rate_residual_" + tag));
  res->SetDirectory(nullptr);
  // Background scale. Live-time normalisation makes this 1.0 by construction,
  // but the two runs are not related by a pure scale: the Ge source gammas
  // fluoresce the Pb shielding, so the SIGNAL run carries induced Pb K x-rays
  // with no counterpart in the background run. That is an ADDITIVE line
  // component in one run only -- k cannot remove it (raising k to kill the Pb
  // residue over-subtracts the continuum everywhere else), which is why the fit
  // window excludes the Pb region instead. k is kept as a cross-check on the
  // live-time normalisation itself, fit on the clean window.
  res->Add(bkg, -bkg_scale);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawHistogram(res, kBlack);
  PlottingUtils::SaveFigure(canvas, "rate_residual_" + tag, "fits",
                            PlotSaveOptions::kLINEAR);

  // CHI-SQUARE fit, NOT FittingUtils.
  //
  // FittingUtils fits with ROOT option "L" -- a Poisson log-likelihood, which
  // assumes bin contents are integer COUNTS. This histogram holds counts/SECOND
  // (values of order 1) and, after subtraction, goes NEGATIVE. Poisson on that
  // is meaningless: a bin holding 2.5 gets an implied error of sqrt(2.5) ~ 1.6
  // instead of its true ~0.03, so parameter errors come out ~20x too large and
  // the quoted chi2 is not a chi2 at all (0.024 on the 01/15 residual, against
  // visible pulls of order 1). BuildRateHist already calls Sumw2(), so the
  // correctly propagated errors are available -- a chi2 fit uses them, copes
  // with negative bins, and gives an interpretable goodness of fit.
  //
  // LINEAR background, not flat: the 01/16 residual does not cancel to a
  // constant (it runs from about -2 counts/s at 64 keV to +2 at 72), and a flat
  // model both fits badly and drags the centroid high. The slope fits to ~0
  // where the subtraction is clean, so this costs the 01/15 case nothing.
  Double_t edge_lo = res->Integral(res->FindBin(lo + 0.2), res->FindBin(66.0));
  Double_t nedge = res->FindBin(66.0) - res->FindBin(lo + 0.2) + 1;
  Double_t bkg0 = (nedge > 0) ? edge_lo / nedge : 0.0;
  Int_t pk_bin = res->GetMaximumBin();
  Double_t pk_val = res->GetBinContent(pk_bin) - bkg0;
  Double_t sigma0 = 0.85;
  // Gaussian core PLUS a low-side exponential tail, as an exponentially
  // modified Gaussian (the exponential convolved with the resolution, which is
  // what the tail physically is).
  //
  // The tail is NOT optional. A bare Gaussian biases the centroid LOW when the
  // peak is tailed -- the shared-tail comment in RunBkgGeSim says so directly
  // -- and the in-situ model this result is compared against carries tails. An
  // earlier version of this function fitted a bare Gaussian, which is very
  // likely why the rate-sub sat ~46 eV below the in-situ once the in-situ shape
  // was freed. Both methods must describe the same lineshape or the
  // "method systematic" is just a comparison of two different biases.
  //
  //   [0] area   [1] mu   [2] sigma   [3] bkg0   [4] bkg slope
  //   [5] tail fraction   [6] tail decay length tau [keV]
  const char *kEmgLowTail =
      "(1-[5])*[0]*TMath::Gaus(x,[1],[2],1)"
      " + [5]*[0]/(2*[6])*TMath::Exp((x-[1])/[6] + [2]*[2]/(2*[6]*[6]))"
      "   *TMath::Erfc((x-[1])/(TMath::Sqrt2()*[2]) + [2]/(TMath::Sqrt2()*[6]))"
      " + [3] + [4]*(x-68.75)";
  TF1 *model = new TF1("ratesub_model_" + tag, kEmgLowTail, lo, hi);
  model->SetParNames("Area", "Mu", "Sigma", "Bkg0", "BkgSlope", "TailFrac",
                     "TailTau");
  model->SetParameters(pk_val * sigma0 * 2.5066, E_GE_73M, sigma0, bkg0, 0.0,
                       0.20, 1.0);
  model->SetParLimits(1, lo + 1.0, hi - 1.0);
  model->SetParLimits(2, 0.2, 3.0);
  // TAIL FIXED, NOT FLOATED.
  //
  // tailFrac railed at a bound in every configuration tried -- 0 with a linear
  // background, 0 with a flat one, 0.6 with the 01/16 slope free -- because
  // after subtraction the effective statistics (S/sqrt(S+B) with a large
  // removed B) cannot resolve a low-amplitude broad tail. A parameter pinned at
  // a bound also reports a meaningless error: floating it gave the 01/15
  // centroid an uncertainty of 1.3 eV, down from 33.6 eV, and that single bogus
  // number then dominated the five-run combination by weight.
  //
  // So the tail is FIXED. The nominal is no tail, which is what this data
  // prefers wherever it can express a preference. The systematic is then the
  // shift when it is instead fixed at the in-situ lineshape's value -- an
  // explicit, quotable number rather than a fit artefact.
  // free_tail overrides the above for a single call (diagnostic scans only);
  // tail_frac/tail_tau < 0 mean "use the global fixed nominal".
  if (free_tail) {
    model->SetParLimits(5, 0.0, 0.6);
    model->SetParLimits(6, 0.2, 5.0);
  } else {
    model->FixParameter(5, tail_frac >= 0 ? tail_frac : RATESUB_TAIL_FRAC);
    model->FixParameter(6, tail_tau > 0 ? tail_tau : RATESUB_TAIL_TAU);
  }
  // FLAT background by default. A linear term and a low-side tail trade against
  // each other almost freely across a 9 keV window, which is how the 01/16 fit
  // ended up with tailFrac railed at its 0.6 ceiling AND a slope of 1.09 -- two
  // parameters describing the same structure. Fixing the slope to zero leaves
  // the tail as the only low-side degree of freedom, so the fit has to commit.
  if (!linear_bkg)
    model->FixParameter(4, 0.0);
  TFitResultPtr fr = res->Fit(model, "S Q R");
  Bool_t ok = fr.Get() && fr->IsValid() && model->GetNDF() > 0;
  if (ok) {
    out.mu = model->GetParameter(1);
    out.mu_err = model->GetParError(1);
    out.chi2 = model->GetChisquare() / model->GetNDF();
    out.valid = kTRUE;
    std::cout << "RATESUB-FIT " << tag << ": mu = " << std::fixed
              << std::setprecision(4) << out.mu << " +/- " << out.mu_err
              << " keV   sigma = " << model->GetParameter(2)
              << "   tailFrac = " << model->GetParameter(5)
              << "   tailTau = " << model->GetParameter(6)
              << "   bkg0 = " << model->GetParameter(3)
              << "   slope = " << model->GetParameter(4)
              << "   k = " << bkg_scale
              << "   chi2/ndf = " << std::setprecision(2) << out.chi2
              << std::endl;
  } else {
    std::cerr << "WARNING: rate-sub chi2 fit failed for " << tag << std::endl;
  }
  // Same styling as every other fit in this analysis: model over the data with
  // the residual pull panel underneath. logy = kFALSE because a subtracted
  // residual legitimately goes negative (the 01/16 spectrum does).
  const Int_t kNPts = 400;
  Double_t xstep = (hi - lo) / (kNPts - 1);
  TGraph *total_graph = new TGraph(kNPts);
  TGraph *bkg_graph = new TGraph(kNPts);
  for (Int_t i = 0; i < kNPts; i++) {
    Double_t xv = lo + i * xstep;
    total_graph->SetPoint(i, xv, model->Eval(xv));
    bkg_graph->SetPoint(
        i, xv, model->GetParameter(3) + model->GetParameter(4) * (xv - 68.75));
  }
  bkg_graph->SetLineStyle(2);
  std::vector<TGraph *> components = {bkg_graph};

  // Clamp the drawn range to the fit window and set the y-range from the
  // in-range bins only. The UNDERFLOW bin holds every event below 64 keV --
  // tens of thousands of counts/s against a signal of order 1 -- and if it is
  // drawn it dictates the y-scale and flattens the spectrum into the axis.
  res->GetXaxis()->SetRangeUser(lo, hi);
  Double_t ymin = 1e30, ymax = -1e30;
  for (Int_t b = res->FindBin(lo); b <= res->FindBin(hi - 1e-6); b++) {
    Double_t v = res->GetBinContent(b);
    Double_t e = res->GetBinError(b);
    if (v - e < ymin)
      ymin = v - e;
    if (v + e > ymax)
      ymax = v + e;
  }
  Double_t pad = 0.12 * (ymax - ymin);
  res->SetMinimum(ymin - pad);
  res->SetMaximum(ymax + pad);
  PlottingUtils::PlotFitWithResiduals(res, total_graph, components, lo, hi,
                                      "ratesub_fit_" + tag, "fits",
                                      "Rate-subtracted " + tag, kFALSE);
  delete sig;
  delete bkg;
  delete res;
  return out;
}

// Process a background/signal pair anchored on the Pb-Ka doublet.
//
// is_reference = kTRUE: this pair defines the day's reference line. Under
//   USE_AM_TRANSFER it is fit from Am-241 (am_run) + the measured Pb points, so
//   68.75 sits BETWEEN anchors (Am 59.5 below, Pb 72.8/75.0/84.9 above) instead
//   of below all of them. With the toggle off it self-cals on Pb alone and Am
//   is only cross-checked.
// is_reference = kFALSE: drift dataset. Only the Pb doublet is measured. Under
//   USE_AM_TRANSFER the reference line is rescaled by this run's Pb-Ka gain
//   (see ref_cal below); with the toggle off the run self-cals on Pb alone.
//
// ref_p0 is VESTIGIAL. It is threaded in from the day functions but never read
//   here, and CreateAndSavePol1Cal is only ever called with fix_p0 = kFALSE.
//   The "held reference offset" that older comments in this file described was
//   never actually wired up: every pair got an independent free-offset Pb fit,
//   including the four that produced the published 68.7606. Gain transfer
//   supersedes that idea rather than reviving it; the parameter is kept only so
//   the day-function signatures stay stable.
// ref_cal / ref_pb_mus / ref_pb_errs are supplied only under USE_AM_TRANSFER,
// and only for drift pairs: the run then rescales the reference line by its own
// Pb-Ka gain instead of fitting a calibration of its own.
PairCalResult ProcessPbAnchoredPair(
    const TString &bkg_run, const TString &sig_run, const TString &date_label,
    const TString &pair_tag, Bool_t interactive, const TString &am_run,
    Bool_t is_reference, Float_t ref_p0, const CalFit *ref_cal = nullptr,
    const std::vector<Float_t> *ref_pb_mus = nullptr,
    const std::vector<Float_t> *ref_pb_errs = nullptr) {
  PairCalResult out;
  out.pair_tag = pair_tag;

  std::vector<Double_t> bkg_events = LoadFiltered(bkg_run);
  std::vector<Double_t> sig_events = LoadFiltered(sig_run);
  std::vector<Double_t> pb_mus = {(Double_t)E_PB_KA2, (Double_t)E_PB_KA1};
  TString cal_label = date_label + "_" + pair_tag;

  FitResult am;
  Bool_t have_am = kFALSE;
  if (is_reference && !am_run.IsNull()) {
    std::vector<Double_t> am_events = LoadFiltered(am_run);
    am = FitCalPeak(am_events, am_run, "Am_59.5keV", 51, 71, kFALSE, 0, 0,
                    interactive, kFALSE, USE_STEP_CAL_PEAKS);
    if (am.valid && !am.peaks.empty()) {
      const PeakFitResult &p = am.peaks.at(0);
      std::cout << "AM-SHAPE step=" << (USE_STEP_CAL_PEAKS ? "on " : "off")
                << "  chi2/ndf " << am.reduced_chi2 << "  sigma " << p.sigma
                << "  step " << p.step_amplitude << "  lowExpAmp "
                << p.low_exp_tail_amplitude << "  lowExpTau "
                << p.low_exp_tail_ratio << "  lowLinAmp "
                << p.low_lin_tail_amplitude << "  lowLinSlope "
                << p.low_lin_tail_slope << std::endl;
    }
    have_am = am.valid;
    if (!have_am)
      std::cerr << "WARNING: Am-241 fit failed for " << am_run
                << "; reference line will use Pb peaks only" << std::endl;
  }

  // Precal uses the SAME tail model as everything downstream: the Pb centroids
  // it returns set the calibration, so the chain must be internally consistent.
  // force_sim_fit = kTRUE: the raw Ge centroid this whole measurement rests on
  // must be MEASURED, not replayed from a past GUI session. The background seed
  // fit above still replays its saved state, which is fine -- it only touches
  // raw data.
  BkgGeSimResult pre =
      RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_run, sig_run,
                  "PreCal_" + pair_tag, interactive, kFALSE, nullptr, nullptr,
                  kFALSE, TAIL_SHARED, CAL_RANGE_LOW, CAL_RANGE_HIGH, kTRUE);
  if (!pre.valid)
    return out;

  // Second pass with the spacing constrained; see CONSTRAIN_PRECAL_SEPARATION.
  // Same fit_label deliberately, so the hand-tuned saved state still seeds this
  // fit -- a cold start here does not converge. That does mean this pass
  // re-saves the state, so snapshot it before running a comparison.
  if (pre.bkg_channel.peaks.size() >= 2) {
    Double_t mu_lo = pre.bkg_channel.peaks.at(0).mu;
    Double_t mu_hi = pre.bkg_channel.peaks.at(1).mu;
    Double_t g = (mu_lo + mu_hi) / (Double_t)(E_PB_KA2 + E_PB_KA1);
    Double_t delta_raw = (Double_t)(E_PB_KA1 - E_PB_KA2) * g;
    // Tabulated spacing error, plus the gain uncertainty propagated onto the
    // spacing: a 5e-4 gain error moves a 2.17 keV spacing by about 1.1 eV.
    Double_t sig_raw = TMath::Sqrt(dE_PB_KA_SEP * dE_PB_KA_SEP +
                                   (delta_raw * 5e-4) * (delta_raw * 5e-4));
    std::cout << "PRECAL-SEP " << pair_tag << ": fitted " << std::fixed
              << std::setprecision(5) << (mu_hi - mu_lo) << "  expected "
              << delta_raw << " +/- " << sig_raw << "  (excess "
              << std::setprecision(1) << ((mu_hi - mu_lo) - delta_raw) * 1000
              << " eV)" << std::setprecision(4) << std::endl;
    if (CONSTRAIN_PRECAL_SEPARATION) {
      BkgGeSimResult pre2 =
          RunBkgGeSim(bkg_events, sig_events, pb_mus, bkg_run, sig_run,
                      "PreCal_" + pair_tag, interactive, kFALSE, nullptr,
                      nullptr, kFALSE, TAIL_SHARED, CAL_RANGE_LOW,
                      CAL_RANGE_HIGH, kTRUE, delta_raw, sig_raw);
      if (pre2.valid)
        pre = pre2;
      else
        std::cerr << "WARNING: constrained pre-cal fit failed for " << pair_tag
                  << "; keeping the unconstrained result." << std::endl;
    }
  }

  for (size_t i = 0; i < pre.bkg_channel.peaks.size(); i++) {
    out.pb_mus_precal.push_back(pre.bkg_channel.peaks[i].mu);
    out.pb_mu_errs_precal.push_back(GuardedMuError(
        pair_tag + " Pb-Ka", pre.bkg_channel.peaks[i], pre.bkg_channel));
  }

  // Pb K-beta doublet (Kb1,3 group ~84.77 and Kb2 group ~87.32): same
  // FitDoublePeak method that seeds the K-alpha, run on the bkg spectrum in the
  // K-beta window. Fit as a doublet to cleanly separate the peaks, but anchor
  // only on Kb1 -- three Pb points (Ka2, Ka1, Kb1) -> a wide lever arm per
  // pair. With USE_AM_TRANSFER off these three ARE the entire calibration for
  // every pair (no Am, nothing shared between runs); with it on they are joined
  // by Am on the reference pair and serve as the gain handle on all the others.
  // link_sigma = kTRUE: the Pb K-beta1/K-beta2 groups share one detector
  // resolution, so tie their widths. Stops the weaker Kb2 from floating its
  // sigma and trading against the continuum, which would skew the Kb1 anchor.
  FitResult kb = FitCalPeak(bkg_events, bkg_run, "Pb_Kbeta", 81, 91, kTRUE,
                            E_PB_KB1, E_PB_KB2, interactive, kTRUE);

  // Kb1 joins the stored Pb centroids so the affine transfer has a THIRD point.
  // Two Ka lines determine an offset+gain map exactly, leaving no residual to
  // estimate an uncertainty from; a third gives 1 degree of freedom.
  if (kb.valid && USE_KB1_ANCHOR && !kb.peaks.empty()) {
    out.pb_mus_precal.push_back(kb.peaks.at(0).mu);
    out.pb_mu_errs_precal.push_back(
        GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb));
  }

  CalibrationData cal_data;
  AddCalPoint(cal_data, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
              GuardedMuError(pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0),
                             pre.bkg_channel),
              E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
  AddCalPoint(cal_data, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
              GuardedMuError(pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1),
                             pre.bkg_channel),
              E_PB_KA1, -1, dE_PB_KA);
  // Anchor on Kb1 only: Kb2 (weak KN line) has a ~0.1 keV centroid error that
  // tilts the unweighted cal fit and drags the Ge extrapolation low.
  if (kb.valid && USE_KB1_ANCHOR)
    AddCalPoint(cal_data, pair_tag + " Pb-Kb1", kb.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb),
                E_PB_KB1, kb.reduced_chi2, dE_PB_KB1);
  else
    std::cerr << "WARNING: Pb K-beta fit failed for " << bkg_run
              << "; calibration falls back to K-alpha only" << std::endl;

  // Three schemes, selected by USE_AM_TRANSFER and by what this pair has:
  //   transfer  - drift pair, rescale the day's Am+Pb reference line by gain
  //   am_anchor - reference pair, fit Am(59.5) + Ka1 + Ka2 + Kb1 -> 68.75 is
  //               bracketed below and above rather than extrapolated
  //   published - Pb-only self-cal (the scheme behind 68.7606)
  Bool_t transfer_mode =
      (USE_AM_TRANSFER && ref_cal != nullptr && ref_pb_mus != nullptr &&
       ref_pb_mus->size() >= 2 && pre.bkg_channel.peaks.size() >= 2);
  Bool_t am_anchor_mode = (USE_AM_TRANSFER && !transfer_mode && have_am);

  CalFit calfit;
  GainResult gain;
  AffineResult affine;

  if (transfer_mode) {
    std::vector<Float_t> mu_run = {(Float_t)pre.bkg_channel.peaks.at(0).mu,
                                   (Float_t)pre.bkg_channel.peaks.at(1).mu};
    std::vector<Float_t> err_run = {
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka2",
                                pre.bkg_channel.peaks.at(0), pre.bkg_channel),
        (Float_t)GuardedMuError(pair_tag + " Pb-Ka1",
                                pre.bkg_channel.peaks.at(1), pre.bkg_channel)};
    if (kb.valid && USE_KB1_ANCHOR && !kb.peaks.empty()) {
      mu_run.push_back((Float_t)kb.peaks.at(0).mu);
      err_run.push_back(
          (Float_t)GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb));
    }
    std::vector<Float_t> er =
        ref_pb_errs ? *ref_pb_errs
                    : std::vector<Float_t>(ref_pb_mus->size(), 0.f);

    if (USE_AFFINE_TRANSFER)
      affine = PbAffineVsReference(mu_run, err_run, *ref_pb_mus, er);
    if (affine.valid) {
      calfit = AffineScaleCal(*ref_cal, affine, cal_label);
      std::cout << "AFFINE-TRANSFER " << cal_label << ": a = " << std::fixed
                << std::setprecision(6) << affine.a << " +/- "
                << std::sqrt(std::max(0.0, affine.var_a))
                << "   b = " << affine.b << " +/- "
                << std::sqrt(std::max(0.0, affine.var_b))
                << "  (n = " << affine.n
                << ", chi2/ndf = " << std::setprecision(2) << affine.chi2_ndf
                << ")" << std::endl;
    } else {
      // Fewer than three usable Pb lines: fall back to the pure-gain model,
      // which needs only two but cannot represent an offset.
      gain = PbGainVsReferenceWithError(mu_run, err_run, *ref_pb_mus, er);
      calfit = GainScaleCal(*ref_cal, gain.g, cal_label);
      std::cout << "GAIN-TRANSFER (fallback) " << cal_label
                << ": g = " << std::fixed << std::setprecision(6) << gain.g
                << " +/- " << gain.sigma_g << "  (n = " << gain.n
                << ", chi2/ndf = " << std::setprecision(2) << gain.chi2_ndf
                << ")" << std::endl;
    }
  } else if (am_anchor_mode) {
    CalibrationData cal_am;
    AddCalPoint(cal_am, pair_tag + " Am-241 59.5", am.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Am-241", am.peaks.at(0), am),
                E_AM241, am.reduced_chi2, dE_AM241);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Ka2",
                               pre.bkg_channel.peaks.at(0), pre.bkg_channel),
                E_PB_KA2, pre.bkg_channel.reduced_chi2, dE_PB_KA);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
                GuardedMuError(pair_tag + " Pb-Ka1",
                               pre.bkg_channel.peaks.at(1), pre.bkg_channel),
                E_PB_KA1, -1, dE_PB_KA);
    if (kb.valid && USE_KB1_ANCHOR)
      AddCalPoint(cal_am, pair_tag + " Pb-Kb1", kb.peaks.at(0).mu,
                  GuardedMuError(pair_tag + " Pb-Kb1", kb.peaks.at(0), kb),
                  E_PB_KB1, kb.reduced_chi2, dE_PB_KB1);
    PrintCalSummary(cal_am, cal_label);
    // pol2: with Am at 59.5 plus the three Pb lines the response curves, and a
    // straight line cannot serve both ends. 68.75 is then INTERPOLATED between
    // Am below and Pb above rather than extrapolated below every anchor.
    calfit = CreateAndSavePol1Cal(cal_am, cal_label, kFALSE, 0, REF_CAL_DEGREE);
  } else {
    PrintCalSummary(cal_data, cal_label);
    calfit = CreateAndSavePol1Cal(cal_data, cal_label, kFALSE, 0);
  }
  TF1 *cal = calfit.func;

  // Pb-only line (Ka1, Ka2, Kb1 -- Am excluded by construction). In the
  // published scheme this IS the nominal; under Am anchoring it is rebuilt here
  // solely so the LINEARITY check below stays non-circular.
  CalFit pb_probe = calfit;
  if (USE_AM_TRANSFER && have_am)
    pb_probe = CreateAndSavePol1Cal(cal_data, cal_label + "_PbOnly", kFALSE, 0);

  // On the Am reference pair, keep the master line so the day's other pairs
  // (and the 14th, which has no Am of its own) can gain-transfer off it.
  if (USE_AM_TRANSFER && have_am) {
    out.ref_calfit = calfit;
    out.has_ref_cal = kTRUE;
  }
  out.cal_func = cal;
  // Calibration-term uncertainty propagated to the Ge energy (cal covariance at
  // the precal Ge centroid). Ge fit-error term is added later when combined.
  if (!pre.sig_channel.peaks.empty()) {
    Double_t mu_ge = pre.sig_channel.peaks.back().mu;
    out.ge_cal_err = CalEnergyError(calfit, mu_ge);

    if (transfer_mode) {
      // A gain error is a pure multiplicative scale error, i.e. slope-like, so
      // it folds into the SLOPE leg of the pivot split. This keeps the split's
      // defining property intact: offset^2 + slope^2 == total^2.
      out.ge_gain_err =
          affine.valid
              ? AffineTransferEnergyError(*ref_cal, mu_ge, affine)
              : GainTransferEnergyError(*ref_cal, mu_ge, gain.g, gain.sigma_g);
      out.ge_cal_err = std::hypot(out.ge_cal_err, out.ge_gain_err);

      std::cout << "GAIN-CAL " << cal_label << ": transfer term " << std::fixed
                << std::setprecision(1) << out.ge_gain_err * 1000.0
                << " eV   total cal " << out.ge_cal_err * 1000.0 << " eV"
                << std::endl;
    }
  }

  // RAW GE/PB RATIO -- calibration-independent by construction. Any spread in
  // this across runs is the FITS disagreeing, not the calibrations; a run whose
  // final energy moves while this stays put has a calibration problem, and vice
  // versa. This is the handle for telling the two apart.
  if (pre.bkg_channel.peaks.size() >= 2 && !pre.sig_channel.peaks.empty()) {
    Double_t mu_pb1 = pre.bkg_channel.peaks[0].mu;
    Double_t mu_pb2 = pre.bkg_channel.peaks[1].mu;
    Double_t denom = mu_pb2 - mu_pb1;
    if (std::fabs(denom) > 0) {
      Double_t ge_raw = pre.sig_channel.peaks.back().mu;
      Double_t r_ge_pb = (ge_raw - mu_pb1) / denom;
      std::cout << "R_GEPB " << cal_label << ": raw Ka2 " << std::fixed
                << std::setprecision(4) << mu_pb1 << "  raw Ka1 " << mu_pb2
                << "  raw Ge " << ge_raw << "  R = " << std::setprecision(6)
                << r_ge_pb << std::endl;
    }
  }

  // Cross-checks on the pair that also has Am. All lines are evaluated at the
  // SAME raw Ge centroid, so the fit cancels and only the calibration differs
  // -- this is the one place the two effects are cleanly separable.
  //
  // Under USE_AM_TRANSFER the nominal line already contains Am, so the Pb-only
  // line is rebuilt separately above (pb_probe) and the LINEARITY test below is
  // scored against THAT: Am must play no part in the line whose Am prediction
  // is being tested, or the check is circular.
  if (have_am && !pre.sig_channel.peaks.empty()) {
    CalibrationData cal_am;
    AddCalPoint(cal_am, pair_tag + " Am-241 59.5", am.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Am-241", am.peaks.at(0), am),
                E_AM241, am.reduced_chi2, dE_AM241);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka2", pre.bkg_channel.peaks.at(0).mu,
                GuardedMuError(pair_tag + " Pb-Ka2",
                               pre.bkg_channel.peaks.at(0), pre.bkg_channel),
                E_PB_KA2, -1, dE_PB_KA);
    AddCalPoint(cal_am, pair_tag + " Pb-Ka1", pre.bkg_channel.peaks.at(1).mu,
                GuardedMuError(pair_tag + " Pb-Ka1",
                               pre.bkg_channel.peaks.at(1), pre.bkg_channel),
                E_PB_KA1, -1, dE_PB_KA);
    CalFit cal_check =
        CreateAndSavePol1Cal(cal_am, cal_label + "_AmCheck", kFALSE, 0);
    Double_t ge_pre = pre.sig_channel.peaks.back().mu;
    // Three lines evaluated at the SAME raw Ge centroid, so the fit cancels and
    // only the calibration differs. nominal-vs-Pb-only is the quantity this
    // whole exercise turns on: what anchoring on Am does to the Ge energy.
    Double_t e_nominal = cal->Eval(ge_pre);
    Double_t e_pbonly = pb_probe.func->Eval(ge_pre);
    Double_t e_amka = cal_check.func->Eval(ge_pre);
    std::cout << "VALIDATION " << pair_tag << " (Ge precal " << std::fixed
              << std::setprecision(4) << ge_pre << "):  nominal -> "
              << e_nominal << " keV   Pb-only -> " << e_pbonly
              << " keV   Am+Ka -> " << e_amka << " keV   (nominal - Pb-only) "
              << std::setprecision(1) << (e_nominal - e_pbonly) * 1000.0
              << " eV   (nominal - Am+Ka) " << (e_nominal - e_amka) * 1000.0
              << " eV" << std::endl;

    // LINEARITY CHECK (for the paper). The Pb-only calibration (Ka1, Ka2, Kb1
    // -- Am NOT used) is extrapolated DOWN to the Am-241 raw centroid and
    // compared to the reference Am energy. The Am point played no role in
    // building this line, so agreement is a non-circular demonstration of
    // linear response across 59.5-85 keV -> the Ge at 68.75 is bracketed (Am
    // below, Pb above), NOT extrapolated. sigma combines the cal covariance at
    // the Am position (dominant) with the Am centroid fit error.
    Double_t am_raw = am.peaks.at(0).mu;
    Double_t am_pred = pb_probe.func->Eval(am_raw);
    Double_t sig_cal = CalEnergyError(pb_probe, am_raw);
    Double_t sig_amfit = pb_probe.p1 * am.peaks.at(0).mu_error;
    Double_t sig_pred = std::sqrt(sig_cal * sig_cal + sig_amfit * sig_amfit);
    Double_t resid = am_pred - E_AM241;
    std::cout << "LINEARITY " << pair_tag
              << ": Pb-only line predicts Am = " << std::fixed
              << std::setprecision(4) << am_pred << " +/- " << sig_pred
              << " keV   (truth " << E_AM241 << ", residual "
              << std::setprecision(1) << resid * 1000.0
              << " eV = " << std::setprecision(2)
              << (sig_pred > 0 ? resid / sig_pred : 0) << " sigma)"
              << std::endl;
  }

  std::vector<TString> pair_runs = {bkg_run, sig_run};
  if (have_am)
    pair_runs.push_back(am_run);
  ApplyPol1Cal(pair_runs, cal, cal_label);

  std::vector<Double_t> bc = LoadCalibrated(bkg_run);
  std::vector<Double_t> sc = LoadCalibrated(sig_run);

  Double_t cal_p0 = cal->GetParameter(0);
  Double_t cal_p1 = cal->GetParameter(1);
  FitResult bkg_seed_post = pre.bkg_channel;
  FitResult sig_seed_post = pre.sig_channel;
  if (cal_p1 > 0) {
    Double_t s_b = pre.bkg_channel.lin_bkg_slope;
    Double_t s_s = pre.sig_channel.lin_bkg_slope;
    bkg_seed_post.lin_bkg_slope = s_b / (cal_p1 - s_b * cal_p0);
    sig_seed_post.lin_bkg_slope = s_s / (cal_p1 - s_s * cal_p0);
    bkg_seed_post.bkg_constant = pre.bkg_channel.bkg_constant / cal_p1;
    sig_seed_post.bkg_constant = pre.sig_channel.bkg_constant / cal_p1;
  }
  // Map the seed peak mus through the calibration so the now-unlocked Pb (and
  // Ge) peaks start at their calibrated positions, not the precal seed values.
  // Map mu through the calibration AND convert the shape into calibrated units.
  //
  // Previously only mu was mapped: sigma and the tail parameters were carried
  // over as fitted in RAW units, and with LOCK_POSTCAL_STATE they were then
  // held FIXED at those raw values while fitting CALIBRATED data. The local
  // calibration derivative is f'(mu) ~ 0.989-0.995, so every width and tail
  // length was ~1% too large (sigma too wide by 7-10 eV) and, with the shape
  // frozen, that mismatch could only be absorbed by mu. Freeing the shape
  // instead moved the Ge centroids up by 24-68 eV, confirming the bias -- but
  // it also destabilised the fits and lost runs, so the fix is to correct the
  // units and KEEP the lock, not to unlock.
  //
  // Scaling under E = f(x) with s = f'(mu), evaluated per peak because the
  // reference calibration is a pol2 and the derivative varies across the range:
  //   sigma, *TailRatio (lengths in x)  -> multiply by s
  //   LowLinTailSlope   (per unit x)    -> divide by s
  //   *TailAmplitude, step (ratios)     -> unchanged
  auto convert_peaks = [&cal](const FitResult &src, FitResult &dst) {
    for (size_t i = 0; i < dst.peaks.size() && i < src.peaks.size(); i++) {
      Double_t mu_raw = src.peaks[i].mu;
      Double_t s = cal->Derivative(mu_raw);
      if (!(s > 0))
        s = 1.0;
      dst.peaks[i].mu = cal->Eval(mu_raw);
      dst.peaks[i].sigma = src.peaks[i].sigma * s;
      dst.peaks[i].low_exp_tail_ratio = src.peaks[i].low_exp_tail_ratio * s;
      dst.peaks[i].high_exp_tail_ratio = src.peaks[i].high_exp_tail_ratio * s;
      dst.peaks[i].low_lin_tail_slope = src.peaks[i].low_lin_tail_slope / s;
    }
  };
  convert_peaks(pre.bkg_channel, bkg_seed_post);
  convert_peaks(pre.sig_channel, sig_seed_post);

  // UNIT-CONVERSION DIAGNOSTIC.
  //
  // mu is mapped through the calibration above, and the background constant and
  // slope are rescaled by cal_p1 -- but sigma and the tail parameters are NOT.
  // They were fitted in RAW units and, with LOCK_POSTCAL_STATE, are then held
  // FIXED at those raw values while fitting CALIBRATED data. Anything carrying
  // units of energy is therefore wrong by the local calibration derivative.
  //
  // Scaling rules under E = f(x), with s = f'(mu) the local derivative (which
  // varies with position now that the reference cal is a pol2, so it must be
  // evaluated per peak rather than using p1 globally):
  //   sigma, *_tail_ratio (decay lengths in x)  ->  multiply by s
  //   low_lin_tail_slope  (per unit x)          ->  divide by s
  //   *_amplitude         (dimensionless ratio) ->  unchanged
  //
  // This block only REPORTS the mismatch; it does not correct it yet.
  if (!pre.sig_channel.peaks.empty()) {
    const PeakFitResult &ge_pre = pre.sig_channel.peaks.back();
    Double_t s = cal->Derivative(ge_pre.mu);
    std::cout << "UNIT-CHECK " << cal_label
              << " (Ge peak, f'(mu) = " << std::fixed << std::setprecision(6)
              << s << "):" << std::endl;
    std::cout << "  sigma            raw " << std::setprecision(5)
              << ge_pre.sigma << "   used " << ge_pre.sigma << "   should be "
              << ge_pre.sigma * s << "   (delta "
              << (ge_pre.sigma - ge_pre.sigma * s) * 1000.0 << " eV)"
              << std::endl;
    std::cout << "  lowExpTailRatio  raw " << ge_pre.low_exp_tail_ratio
              << "   used " << ge_pre.low_exp_tail_ratio << "   should be "
              << ge_pre.low_exp_tail_ratio * s << std::endl;
    std::cout << "  lowLinTailSlope  raw " << ge_pre.low_lin_tail_slope
              << "   used " << ge_pre.low_lin_tail_slope << "   should be "
              << (s != 0 ? ge_pre.low_lin_tail_slope / s : 0.0) << std::endl;
    std::cout << "  highExpTailRatio raw " << ge_pre.high_exp_tail_ratio
              << "   used " << ge_pre.high_exp_tail_ratio << "   should be "
              << ge_pre.high_exp_tail_ratio * s << std::endl;
  }

  // Pb peaks left FREE post-cal (lock_bkg_peaks = kFALSE): where the calibrated
  // Pb doublet actually lands is a blind check at the energy of interest. The
  // whole chain (precal -> cal -> postcal) uses ONE tail model (TAIL_SHARED),
  // so this run is internally consistent; the alternative tail model is a
  // separate macro and the lineshape systematic is taken in the combiner.
  // force_sim_fit = kTRUE, same as the pre-cal. With POSTCAL_INTERACTIVE the
  // saved .simroofits is loaded -- which is what we want, since it carries the
  // hand-tuned starting point the automated cold start cannot reach -- but
  // WITHOUT this flag the loader returns those parameters AS THE RESULT and no
  // minimisation runs at all. That is how three post-cal rows came back holding
  // day-one constants with chi2 of 7-9: old parameters scored against current
  // data. With the flag, the saved state SEEDS a real fit.
  out.postcal = RunBkgGeSim(
      bc, sc, pb_mus, bkg_run + "_postcal_" + pair_tag,
      sig_run + "_postcal_" + pair_tag, "PostCal_" + pair_tag,
      POSTCAL_INTERACTIVE, LOCK_POSTCAL_PB, &bkg_seed_post, &sig_seed_post,
      LOCK_POSTCAL_STATE, TAIL_SHARED, CAL_RANGE_LOW, CAL_RANGE_HIGH, kTRUE,
      CONSTRAIN_PB_KA_SEPARATION ? (Double_t)(E_PB_KA1 - E_PB_KA2) : -1.0,
      dE_PB_KA_SEP);

  if (out.postcal.valid) {
    std::cout << "Post-cal FREE Pb mu for " << cal_label << ":" << std::endl;
    for (size_t i = 0;
         i < out.postcal.bkg_channel.peaks.size() && i < pb_mus.size(); i++)
      std::cout << "  Pb peak " << i << ": " << std::fixed
                << std::setprecision(4) << out.postcal.bkg_channel.peaks[i].mu
                << " keV (truth " << pb_mus[i] << ", shift "
                << (out.postcal.bkg_channel.peaks[i].mu - pb_mus[i]) << ")"
                << std::endl;

    // BACKGROUND / CONTINUUM systematic = | cal(precal Ge mu) - postcal Ge mu
    // |. Both are already-converged fits (chi2/ndf ~ 0.3-1.5), so no new
    // minimization to destabilize. The precal fit determines its continuum in
    // raw space with a fully free peak shape; the postcal fit re-determines the
    // continuum independently in calibrated space (background FREE). Under a
    // linear cal these two Ge energies would coincide if nothing differed, so
    // their difference measures the Ge centroid's sensitivity to how the
    // continuum was estimated. (It also folds in the precal->postcal
    // shape-lock, making this a slightly CONSERVATIVE bound -- the safe
    // direction.)
    if (!pre.sig_channel.peaks.empty()) {
      Double_t ge_precal_mapped = cal->Eval(pre.sig_channel.peaks.back().mu);
      Double_t ge_postcal = out.postcal.sig_channel.peaks.back().mu;
      out.ge_bkg_err = std::fabs(ge_precal_mapped - ge_postcal);
      std::cout << "BACKGROUND " << cal_label
                << ": cal(precal Ge)=" << std::fixed << std::setprecision(4)
                << ge_precal_mapped << "  postcal Ge=" << ge_postcal
                << "  |diff|=" << std::setprecision(1)
                << out.ge_bkg_err * 1000.0 << " eV" << std::endl;
    }
  }
  return out;
}

void AppendPairToDay(DayResult &result, const PairCalResult &p,
                     const TString &date_label, const TString &ge_label) {
  if (p.cal_func) {
    result.cal_labels.push_back(date_label + " " + p.pair_tag);
    result.cal_funcs.push_back(p.cal_func);
  }
  AddGeResult(result, p.postcal, ge_label, p.ge_cal_err, p.ge_gain_err,
              p.ge_bkg_err);
}

void PrintPrecalPbComparison(const TString &date_label,
                             const std::vector<PairCalResult *> &pairs) {
  std::cout << "Pre-cal Pb mu comparison for " << date_label << std::endl;
  std::cout << std::left << std::setw(14) << "Pair" << std::right
            << std::setw(14) << "Pb-Ka2 mu" << std::setw(14) << "shift vs true"
            << std::setw(14) << "Pb-Ka1 mu" << std::setw(14) << "shift vs true"
            << std::endl;
  for (size_t i = 0; i < pairs.size(); i++) {
    const PairCalResult *p = pairs[i];
    if (!p || p->pb_mus_precal.size() < 2)
      continue;
    Float_t mu1 = p->pb_mus_precal[0];
    Float_t mu2 = p->pb_mus_precal[1];
    std::cout << std::left << std::setw(14) << p->pair_tag << std::right
              << std::fixed << std::setprecision(4) << std::setw(14) << mu1
              << std::setw(14) << (mu1 - E_PB_KA2) << std::setw(14) << mu2
              << std::setw(14) << (mu2 - E_PB_KA1) << std::endl;
  }
}

DayResult ProcessLineCalDay(const LineCalConfig &cfg, Bool_t interactive) {
  CalibrationData cal_data;

  std::vector<Double_t> am_events = LoadFiltered(cfg.am_run);
  FitResult am = FitCalPeak(am_events, cfg.am_run, "Am_59.5keV", cfg.am_lo,
                            cfg.am_hi, kFALSE, 0, 0, interactive);
  if (am.valid)
    AddCalPoint(cal_data, "Am-241 59.5", am.peaks.at(0).mu,
                am.peaks.at(0).mu_error, E_AM241, am.reduced_chi2, dE_AM241);

  // Ba-133 53.16 is the worst-behaved anchor on both line-cal days: chi2/ndf
  // 1.73 (01/15) and 2.14 (01/16), the poorest of any anchor, with a 42 eV
  // centroid error on the 15th against Am's 9 eV. Excluded by default. Am
  // (59.54) and Ba 80.98 still bracket the Ge line at 68.75, so it remains
  // INTERPOLATED; the cost is that two anchors exactly determine a pol1, giving
  // zero degrees of freedom and therefore no calibration chi2 to inspect.
  const Bool_t USE_BA133_53 = kFALSE;

  std::vector<Double_t> ba_events = LoadFiltered(cfg.ba_run);
  if (USE_BA133_53) {
    FitResult ba53 = FitCalPeak(ba_events, cfg.ba_run, "Ba_53.16keV", 48, 58,
                                kFALSE, 0, 0, interactive);
    if (ba53.valid)
      AddCalPoint(cal_data, "Ba-133 53.16", ba53.peaks.at(0).mu,
                  ba53.peaks.at(0).mu_error, E_BA133_53, ba53.reduced_chi2,
                  dE_LINE);
  }

  // Ba-133 79.61 / 80.98: fitted as a doublet so the weak satellite is MODELLED
  // (it carries ~2.65% intensity against the main line's ~34%, and sits 1.384
  // keV below it -- about -1.6 sigma, close enough to drag a bare single-peak
  // centroid low by O(100 eV)), but only the strong 80.98 line is used as a
  // calibration anchor. Exactly the Pb K-beta recipe elsewhere in this file.
  //
  // link_sigma = kTRUE is the fix for the 20260115 blow-up: unlinked, the weak
  // 79.61 component inflates its own width and trades against the continuum,
  // which compressed the measured doublet separation to 0.998 keV against a
  // true 1.384 keV and drove the line cal to chi2/ndf = 213.
  FitResult ba81 =
      FitCalPeak(ba_events, cfg.ba_run, "Ba_80.98keV", 75, 90, kTRUE,
                 E_BA133_79, E_BA133_81, interactive, kTRUE);
  if (ba81.valid)
    AddCalPoint(cal_data, "Ba-133 80.98", ba81.peaks.at(1).mu,
                ba81.peaks.at(1).mu_error, E_BA133_81, ba81.reduced_chi2,
                dE_LINE);

  PrintCalSummary(cal_data, cfg.date_label);
  CalFit calfit = CreateAndSavePol1Cal(cal_data, cfg.date_label);
  TF1 *cal = calfit.func;
  // Calibration-term uncertainty at the Ge energy. The line cal is in true-keV
  // already, so evaluate the covariance at the precalibrated mu that maps to
  // ~68.75 keV: mu_at_Ge = (E_GE - p0)/p1.
  Double_t mu_at_ge =
      (calfit.p1 != 0) ? (E_GE_73M - calfit.p0) / calfit.p1 : E_GE_73M;
  Double_t ge_cal_err = CalEnergyError(calfit, mu_at_ge);
  ApplyPol1Cal(cfg.day_datasets, cal, cfg.date_label);

  DayResult result;
  result.cal_labels.push_back(cfg.date_label);
  result.cal_funcs.push_back(cal);

  // The 15th/16th are measured by live-time RATE SUBTRACTION only -- no in-situ
  // simultaneous fit on these days. (The 13th/14th use the simultaneous fit;
  // the cross-check is the 13/14 in-situ result vs the 15/16 rate-sub result.)
  // NOT interactive. FittingUtils' interactive path blocks on a GUI it cannot
  // render headless: the macro sleeps at 0% CPU indefinitely after emitting
  // "gPad has at least one zero dimension" warnings. Unlike the RooFitUtils
  // line fits above -- which replay saved .roofits state -- there is nothing to
  // replay here, since the residual histogram is rebuilt from scratch on every
  // run, so this fit must genuinely run each time.
  RateSubResult rs =
      FitRateSubtractedGe(cfg.postcal_sig, cfg.postcal_bkg, cal, cfg.date_label,
                          kFALSE, cfg.rs_lo, cfg.rs_hi, cfg.rs_linear_bkg);
  if (rs.valid) {
    result.rs_labels.push_back(cfg.ge_label + " [rate-sub]");
    result.rs_mus.push_back(rs.mu);
    result.rs_errs.push_back(rs.mu_err);
    result.rs_chi2s.push_back(rs.chi2);
    result.rs_cal_errs.push_back(ge_cal_err);
    result.rs_gain_errs.push_back(0.0f);
    std::cout << "Rate-sub Ge mu for " << cfg.date_label << ": " << std::fixed
              << std::setprecision(4) << rs.mu << " +/- " << rs.mu_err
              << " keV (chi2/ndf = " << std::setprecision(3) << rs.chi2 << ")"
              << std::endl;
  }
  return result;
}

struct ShieldPair {
  TString bkg_run;
  TString sig_run;
  TString tag;
  TString ge_label;
};

// Under USE_AM_TRANSFER the 13th and 14th share one reference line: Am-241
// (59.5) + the Pb points of the 13th dataset closest in time to the Am run.
// Every other pair -- and the 14th, which has no Am run of its own -- rescales
// that line by its own measured Pb-Ka gain.
//
// With the toggle OFF nothing is shared: each pair self-cals on Pb alone. That
// is how the published 68.7606 was produced, notwithstanding what the older
// version of this comment claimed. ref_p0_out is still returned, but only for
// signature stability -- see the note on ref_p0 in ProcessPbAnchoredPair.
DayResult ProcessDay_20260113(Bool_t interactive, Float_t &ref_p0_out,
                              PairCalResult &ref_pair_out) {
  DayResult result;

  std::vector<ShieldPair> pairs;
  ShieldPair cu;
  cu.bkg_run = Constants::CUSHIELDBACKGROUND_10PERCENT_20260113;
  cu.sig_run = Constants::CUSHIELDSIGNAL_10PERCENT_20260113;
  cu.tag = "CuShield10";
  cu.ge_label = "Cu Shield Signal 10% (01/13)";
  pairs.push_back(cu);
  ShieldPair cd10;
  cd10.bkg_run = Constants::CDSHIELDBACKGROUND_10PERCENT_20260113;
  cd10.sig_run = Constants::CDSHIELDSIGNAL_10PERCENT_20260113;
  cd10.tag = "CdShield10";
  cd10.ge_label = "Cd Shield Signal 10% (01/13)";
  pairs.push_back(cd10);
  ShieldPair cd25;
  cd25.bkg_run = Constants::CDSHIELDBACKGROUND_25PERCENT_20260113;
  cd25.sig_run = Constants::CDSHIELDSIGNAL_25PERCENT_20260113;
  cd25.tag = "CdShield25";
  cd25.ge_label = "Cd Shield Signal 25% (01/13)";
  pairs.push_back(cd25);

  TString am_run = Constants::POSTREACTOR_AM241_20260113;
  Float_t am_time = LoadRunStartTime(am_run);
  Int_t ref_idx = 0;
  Float_t best_dt = 1e30f;
  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    Float_t dt = TMath::Abs(LoadRunStartTime(pairs[i].bkg_run) - am_time);
    if (dt < best_dt) {
      best_dt = dt;
      ref_idx = i;
    }
  }
  std::cout << "Reference pair for 13/14 drift: " << pairs[ref_idx].tag
            << " (Am from " << am_run << ")" << std::endl;

  std::vector<PairCalResult> results(pairs.size());

  results[ref_idx] = ProcessPbAnchoredPair(
      pairs[ref_idx].bkg_run, pairs[ref_idx].sig_run, "20260113",
      pairs[ref_idx].tag, interactive, am_run, kTRUE, 0);
  ref_p0_out = results[ref_idx].cal_func
                   ? results[ref_idx].cal_func->GetParameter(0)
                   : 0;

  // The day's reference line plus the Pb-Ka centroids measured in that same
  // run. Under USE_AM_TRANSFER every other pair rescales this line by its own
  // Pb-Ka gain; when the toggle is off these stay null and each pair self-cals.
  const CalFit *ref_cal =
      results[ref_idx].has_ref_cal ? &results[ref_idx].ref_calfit : nullptr;
  const std::vector<Float_t> *ref_mus = &results[ref_idx].pb_mus_precal;
  const std::vector<Float_t> *ref_errs = &results[ref_idx].pb_mu_errs_precal;

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
    if (i == ref_idx)
      continue;
    results[i] = ProcessPbAnchoredPair(
        pairs[i].bkg_run, pairs[i].sig_run, "20260113", pairs[i].tag,
        interactive, "", kFALSE, ref_p0_out, ref_cal, ref_mus, ref_errs);
  }

  for (Int_t i = 0; i < (Int_t)pairs.size(); i++)
    AppendPairToDay(result, results[i], "20260113", pairs[i].ge_label);

  // Rate-subtracted measurement of the same three pairs; see RATESUB_ALL_DAYS.
  if (RATESUB_ALL_DAYS)
    for (Int_t i = 0; i < (Int_t)pairs.size(); i++) {
      if (!results[i].cal_func)
        continue;
      RateSubResult rs = FitRateSubtractedGe(
          pairs[i].sig_run, pairs[i].bkg_run, results[i].cal_func,
          "20260113_" + pairs[i].tag, kFALSE, 62.0, 71.5, kTRUE);
      if (!rs.valid)
        continue;
      result.rs_labels.push_back(pairs[i].ge_label + " [rate-sub]");
      result.rs_mus.push_back(rs.mu);
      result.rs_errs.push_back(rs.mu_err);
      result.rs_chi2s.push_back(rs.chi2);
      result.rs_cal_errs.push_back(results[i].ge_cal_err);
      result.rs_gain_errs.push_back(results[i].ge_gain_err);
      std::cout << "Rate-sub Ge mu for 20260113_" << pairs[i].tag << ": "
                << std::fixed << std::setprecision(4) << rs.mu << " +/- "
                << rs.mu_err << " keV (chi2/ndf = " << std::setprecision(3)
                << rs.chi2 << ")" << std::endl;
    }

  std::vector<PairCalResult *> ptrs;
  for (Int_t i = 0; i < (Int_t)pairs.size(); i++)
    ptrs.push_back(&results[i]);
  PrintPrecalPbComparison("20260113", ptrs);
  // Hand the reference pair to the 14th, which has no Am run of its own.
  ref_pair_out = results[ref_idx];
  return result;
}

// NOTE: there is no Am-241 measurement on the 14th (Constants has only
// POSTREACTOR_AM241_20260113), so this pair can never be Am-anchored directly.
// It transfers the 13th's reference line ACROSS A DAY BOUNDARY, where the
// multiplicative-drift assumption is weakest. GainTransferEnergyError inflates
// this run's calibration error accordingly rather than the run being dropped.
DayResult ProcessDay_20260114(Bool_t interactive, Float_t ref_p0,
                              const PairCalResult &ref_pair) {
  DayResult result;
  const CalFit *ref_cal = ref_pair.has_ref_cal ? &ref_pair.ref_calfit : nullptr;
  PairCalResult cu = ProcessPbAnchoredPair(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, "20260114", "CuShield10",
      interactive, "", kFALSE, ref_p0, ref_cal, &ref_pair.pb_mus_precal,
      &ref_pair.pb_mu_errs_precal);
  AppendPairToDay(result, cu, "20260114", "Cu Shield Signal 10% (01/14)");

  // Second, independent measurement of the same run. The in-situ fit ties the
  // Ge centroid to the Pb doublet through the linked Ka shape parameters, and
  // on this run that doublet comes out ~47 eV wide against a truth known to
  // <1 eV -- so a Pb inconsistency propagates straight into Mu3. Live-time
  // rate subtraction breaks the link: the background run is subtracted off and
  // the Ge line is fit alone, with no Pb constraint on its shape at all.
  //
  // The 15th/16th rate-subs are statistics-starved and cannot constrain a tail
  // (hence the FIXED RATESUB_TAIL_* above). This run is the largest of the
  // campaign, so its residual is the one place that assumption can actually be
  // tested. Window is the 64-73 default, NOT the 16th's tighter 66-71.5: the
  // Pb-Ka2 edge at 72.8 needs no avoiding here, because Pb sits in the
  // background run too and the subtraction removes it. A 5.5 keV window around
  // a sigma ~0.8 peak leaves only ~3.4 sigma of baseline per side, far too
  // little lever arm to separate a linear background from a low tail -- which
  // is the degeneracy that railed tailFrac at both bounds in the scan below.
  if (cu.cal_func) {
    // DIAGNOSTIC ONLY -- prints, does not push a result row.
    //
    // WINDOW scan first, because the window turned out to matter more than the
    // tail. At 66-71.5 this fit has chi2/ndf 1.59; widened to 64-73 it goes to
    // 11.89 with sigma collapsing 0.80 -> 0.655. So the residual is NOT a clean
    // peak on a linear background across the wider range -- something out there
    // is unmodelled, and the prime suspect is Pb-Ka2 at 72.804 failing to
    // cancel because the signal and background runs do not have identical Pb
    // rates. Extending LOW costs nothing (no line down there) and buys the
    // baseline lever arm that a 5.5 keV window around a sigma ~0.8 peak lacks;
    // extending HIGH walks into the Pb residue. This scan separates the two
    // sides so the choice is measured rather than assumed.
    if (RATESUB_TAIL_SCAN) {
      // Background-scale cross-check on the CLEAN window. k = 1 is the
      // live-time normalisation; if the minimum sits elsewhere the live-time
      // ratio is off. This cannot and does not address the induced Pb
      // fluorescence, which is additive and confined to the excluded region.
      const Double_t scan_k[] = {0.94, 0.97, 1.00, 1.03, 1.06};
      for (size_t k = 0; k < sizeof(scan_k) / sizeof(scan_k[0]); k++)
        FitRateSubtractedGe(
            Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, cu.cal_func,
            TString::Format("20260114_k%03d", (Int_t)(scan_k[k] * 100)), kFALSE,
            62.0, 71.5, kTRUE, kFALSE, 0.0, 1.0, scan_k[k]);
      // Then the tail, on the widest window that stayed clean on the low side.
      const Double_t scan_fr[] = {0.0, 0.05, 0.10, 0.20, 0.30};
      for (size_t k = 0; k < sizeof(scan_fr) / sizeof(scan_fr[0]); k++)
        FitRateSubtractedGe(
            Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, cu.cal_func,
            TString::Format("20260114_fix%03d", (Int_t)(scan_fr[k] * 100)),
            kFALSE, 62.0, 71.5, kTRUE, kFALSE, scan_fr[k], 1.0);
      FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                          Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                          cu.cal_func, "20260114_freetail_flat", kFALSE, 62.0,
                          71.5, kFALSE, kTRUE);
      FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                          Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                          cu.cal_func, "20260114_freetail_lin", kFALSE, 62.0,
                          71.5, kTRUE, kTRUE);
    }
    RateSubResult rs =
        FitRateSubtractedGe(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                            Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                            cu.cal_func, "20260114", kFALSE, 62.0, 71.5, kTRUE);
    if (rs.valid) {
      result.rs_labels.push_back("Cu Shield Signal 10% (01/14) [rate-sub]");
      result.rs_mus.push_back(rs.mu);
      result.rs_errs.push_back(rs.mu_err);
      result.rs_chi2s.push_back(rs.chi2);
      result.rs_cal_errs.push_back(cu.ge_cal_err);
      result.rs_gain_errs.push_back(cu.ge_gain_err);
      std::cout << "Rate-sub Ge mu for 20260114: " << std::fixed
                << std::setprecision(4) << rs.mu << " +/- " << rs.mu_err
                << " keV (chi2/ndf = " << std::setprecision(3) << rs.chi2 << ")"
                << std::endl;
    }
  }
  return result;
}

DayResult ProcessDay_20260115(Bool_t interactive) {
  LineCalConfig cfg;
  cfg.date_label = "20260115";
  cfg.am_run = Constants::POSTREACTOR_AM241_20260115;
  cfg.am_lo = 51;
  cfg.am_hi = 71;
  cfg.ba_run = Constants::POSTREACTOR_BA133_20260115;
  cfg.day_datasets = {Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::POSTREACTOR_AM241_20260115,
                      Constants::POSTREACTOR_BA133_20260115,
                      Constants::SHUTTERCLOSED_20260115};
  cfg.postcal_bkg = Constants::NOSHIELDBACKGROUND_5PERCENT_20260115;
  cfg.postcal_sig = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;
  cfg.ge_label = "No Shield Signal 5% (01/15)";
  return ProcessLineCalDay(cfg, interactive);
}

DayResult ProcessDay_20260116(Bool_t interactive) {
  LineCalConfig cfg;
  cfg.date_label = "20260116";
  cfg.am_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.am_lo = 55;
  cfg.am_hi = 70;
  cfg.ba_run = Constants::POSTREACTOR_AM241_BA133_20260116;
  cfg.day_datasets = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  cfg.postcal_bkg =
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116;
  cfg.postcal_sig = Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116;
  cfg.ge_label = "Graphite Castle Signal 10% (01/16)";
  // Tightened from the default 64-73. This day's subtraction leaves an
  // over-subtracted shoulder below ~66 keV and unsettled structure above ~71.5;
  // 66-71.5 still spans the Ge peak at 68.75 to beyond +3 sigma (sigma ~0.85),
  // so nothing of the line is lost.
  cfg.rs_lo = 62.0;
  cfg.rs_hi = 71.5;
  return ProcessLineCalDay(cfg, interactive);
}

void PrintDayCals(const DayResult &d) {
  for (size_t i = 0; i < d.cal_funcs.size(); i++) {
    TF1 *c = d.cal_funcs[i];
    if (!c)
      continue;
    std::cout << "Cal pol1 " << d.cal_labels[i] << ": p0 = " << std::fixed
              << std::setprecision(6) << c->GetParameter(0) << " +/- "
              << c->GetParError(0) << ", p1 = " << c->GetParameter(1) << " +/- "
              << c->GetParError(1) << std::endl;
  }
}

void CalibrationLow() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/calibration_low_" +
                                    RESULT_TAG,
                                project_root + "/root_files");
  Bool_t interactive = kTRUE;

  Float_t ref_p0 = 0;
  PairCalResult ref_pair;
  DayResult d13 = ProcessDay_20260113(interactive, ref_p0, ref_pair);
  DayResult d14 = ProcessDay_20260114(interactive, ref_p0, ref_pair);
  DayResult d15 = ProcessDay_20260115(interactive);
  DayResult d16 = ProcessDay_20260116(interactive);

  PrintDayCals(d13);
  PrintDayCals(d14);
  PrintDayCals(d15);
  PrintDayCals(d16);

  // 15th/16th: line-cal (Am + Ba-133) days measured by live-time RATE
  // SUBTRACTION. Independent of the 13th/14th along both axes -- different
  // calibration sources and a different background treatment -- so they are the
  // cross-check that makes the in-situ result defensible, and they revive the
  // in-situ-vs-rate-sub method systematic in CombineGeResult.
  // Emit one .result row per run. CombineGeResult.cpp reads this file (+ the
  // Ka1-only one) for the budget. This macro does NO budget math -- only rows.
  // 01/14 is EXCLUDED from the combination. It is the only run that cannot be
  // Am-anchored (Constants has no Am-241 measurement on the 14th), so it is the
  // only cross-day gain transfer, where the multiplicative-drift assumption is
  // weakest. It has carried the worst chi2, cal_err and bkg_err of any run
  // throughout, and it is the only run whose Pb Ka doublet rejects a shared
  // width -- linking the sigmas moved its Ka2 centroid 116 eV and inflated its
  // doublet separation 2-3% above every other run, consistent with rate-
  // dependent pileup broadening at its 13.8 M events.
  //
  // Note it pulled the combined value TOWARD the literature number (68.8016 vs
  // the 01/13 cluster at 68.694-68.709). That is not a reason to keep it; the
  // three 01/13 runs agree to 15 eV among themselves and 01/14 sits 100 eV
  // above all of them.
  //
  // Flip INCLUDE_20260114 to restore it as a cross-check. It is also the
  // slowest fit in the chain, so leaving it off saves real time.
  const Bool_t INCLUDE_20260114 = kTRUE;

  std::vector<DayResult *> day_vec = {&d13};
  if (INCLUDE_20260114)
    day_vec.push_back(&d14);
  day_vec.push_back(&d15);
  day_vec.push_back(&d16);
  DayResult **days = day_vec.data();
  const Int_t n_days = (Int_t)day_vec.size();

  const TString results_dir = project_root + "/results";
  gSystem->mkdir(results_dir, kTRUE);
  const TString result_path =
      results_dir + "/ge_" + RESULT_TAG + RESULT_VARIANT + ".result";
  std::ofstream out(result_path.Data());
  out << "# Ge-73m gamma energy results   tag=" << RESULT_TAG << std::endl;
  out << "# method  ge_mu  fit_err  cal_err  bkg_err  chi2  gain_err  label"
      << std::endl;
  out << std::fixed << std::setprecision(6);

  auto sanitize = [](TString s) {
    s.ReplaceAll(" ", "_");
    return s;
  };

  std::cout << "\nMethod 1 -- in-situ simultaneous fit (Ge Peak mu, post-cal):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->ge_labels.size(); i++) {
      Float_t ce =
          i < days[d]->ge_cal_errs.size() ? days[d]->ge_cal_errs[i] : 0;
      Float_t be =
          i < days[d]->ge_bkg_errs.size() ? days[d]->ge_bkg_errs[i] : 0;
      Float_t gz =
          i < days[d]->ge_gain_errs.size() ? days[d]->ge_gain_errs[i] : 0;
      out << "insitu  " << days[d]->ge_mus[i] << "  " << days[d]->ge_errs[i]
          << "  " << ce << "  " << be << "  " << days[d]->ge_chi2s[i] << "  "
          << gz << "  " << sanitize(days[d]->ge_labels[i]) << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->ge_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->ge_mus[i]
                << " +/- " << days[d]->ge_errs[i] << " (fit) +/- " << ce
                << " (cal) keV (chi2/ndf = " << std::setprecision(3)
                << days[d]->ge_chi2s[i] << ")" << std::endl;
    }
  }

  std::cout << "\nMethod 2 -- live-time rate subtraction (Ge Peak mu):"
            << std::endl;
  for (Int_t d = 0; d < n_days; d++) {
    for (size_t i = 0; i < days[d]->rs_labels.size(); i++) {
      Float_t ce =
          i < days[d]->rs_cal_errs.size() ? days[d]->rs_cal_errs[i] : 0;
      // Same 5-column layout as the insitu rows. There is no |cal(pre)-post|
      // for the rate-sub path (no pre/post pair), so bkg_err is 0 here.
      out << "ratesub  " << days[d]->rs_mus[i] << "  " << days[d]->rs_errs[i]
          << "  " << ce << "  " << 0.0 << "  " << days[d]->rs_chi2s[i] << "  "
          << (i < days[d]->rs_gain_errs.size() ? days[d]->rs_gain_errs[i]
                                               : 0.0f)
          << "  " << sanitize(days[d]->rs_labels[i]) << std::endl;
      std::cout << std::left << std::setw(50) << days[d]->rs_labels[i] << ": "
                << std::fixed << std::setprecision(4) << days[d]->rs_mus[i]
                << " +/- " << days[d]->rs_errs[i] << " (fit) +/- " << ce
                << " (cal) keV" << std::endl;
    }
  }
  out.close();
  std::cout << "\nWrote " << result_path << std::endl;
  std::cout << "Run CombineGeResult.cpp for the budget." << std::endl;
}

// Export the spectrum + fit + residuals of one shielded pair to CSV for the
// collaborator's plotting. Reuses the EXACT postcal fit path (loads the cached
// converged fit non-interactively, so it reproduces the published curve), then
// the in-process CSV hook dumps <prefix>_{bkg,sig}.csv. Driven from the thin
// ExportFitCSV.cpp entry point. Each CSV:
//   block 1: smooth fit curve   energy_keV, fit_total, fit_background
//   block 2: binned data + pull bin_center_keV, data_counts, fit_total,
//   residual_pull
void RunFitCSVExport() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/calibration_low_" +
                                    RESULT_TAG,
                                project_root + "/root_files");
  gSystem->mkdir(project_root + "/results", kTRUE);

  // The sig channel (Pb-Ka2, Pb-Ka1, Ge) is the spectrum of interest. Process
  // the reference pair INTERACTIVELY (kTRUE): with the saved
  // .roofits/.simroofits caches present, interactive mode LOADS the converged
  // params and skips the editor -- reproducing the published Boggs-Pike fit
  // (tails and all). Passing kFALSE instead would trigger the automated
  // group-and-prune path, which re-fits from scratch and strips the photopeaks
  // to bare Gaussians. The CSV hook then fires inside RunBkgGeSim's postcal
  // call.
  G_CSV_PREFIX = project_root + "/results/spectrum_CuShield10_20260113";
  Float_t ref_p0 = 0;
  ProcessPbAnchoredPair(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                        Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                        "20260113", "CuShield10", kTRUE,
                        Constants::POSTREACTOR_AM241_20260113, kTRUE, ref_p0);
  G_CSV_PREFIX = "";
  std::cout << "\nCSV export done. The sig file is the Pb-Ka + Ge spectrum; "
               "the bkg file is the Pb-only background channel."
            << std::endl;
}
