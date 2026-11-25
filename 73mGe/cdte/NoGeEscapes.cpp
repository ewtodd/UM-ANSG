// The 95.902 keV line's K X-ray escapes in the PSI CdTe beam runs without a
// Ge target in front of the CdTe, which fix the escape tie of the Ge budget
// (GeBudget.cpp) by measurement instead of Geant4. The in-run model
// (CdTeAnalysis.hpp) without the Ge line, the four escape families free: Cd
// Ka, Cd Kb, Te Ka, Te Kb, each with its sub-lines at fixed relative
// intensities, at the measured capture line minus the X-ray energies. Same
// rise-time slices as the budget.
// Run groups: 10/10 Beam_NOGe; 10/11-12 CdTe_NoGe_Si_YesGe (a Ge target in
// front of the Si, so the CdTe may still see some 68.75 keV); 10/13 NoPb with
// nothing in front of either detector.
// The Cd Ka escapes (72.73, 72.92 keV) sit under Pb Ka2 (72.80 keV), so with
// the Pb Ka2/Ka1 ratio free they are not measured; the Kb/Ka ratio comes from
// fits with Pb Ka2/Ka1 fixed at 0.595 (DDEP) and from the NoPb run.
//
// Usage: root -b -q 'NoGeEscapes.cpp+'
#include "CdTeAnalysis.hpp"

#include <TCanvas.h>
#include <TGaxis.h>
#include <TLegend.h>

namespace CdTeNoGe {

using namespace CdTe;

// Escape family yields (counts); the Cd Kb yield is Q_ESC_A.
enum { R_CD_KA = N_Q, R_TE_KA, R_TE_KB, R_GE_FIXED_E_A, N_R };

// The parameter holding each escape line's family yield, in
// ReferenceLines::ESCAPE_XRAY_KEV's order.
const Int_t FAMILY_PAR[ReferenceLines::N_ESCAPE_LINES] = {
    R_CD_KA, R_CD_KA, Q_ESC_A, Q_ESC_A, R_TE_KA, R_TE_KA, R_TE_KB, R_TE_KB};
const Double_t PB_KA2_OVER_KA1_DDEP = 0.595;

Bool_t g_fix_pb_ratio = kFALSE;
// A 68.752 keV line of free yield, as a contamination check.
Bool_t g_free_ge = kFALSE;
// Without Pb lines the scale cannot be fitted: hold its gain at the Pb runs'
// value and its offset at the capture line.
Bool_t g_anchor_capture = kFALSE;
Double_t g_b_fixed = 1.0;

Double_t EscapeFamily(Double_t x, const Double_t *p, Int_t first, Int_t last) {
  Double_t e_capture = ScaleEnergy(g_capture_meas, p);
  Double_t v = 0;
  for (Int_t k = first; k <= last; k++)
    v += p[FAMILY_PAR[k]] * GEANT4_SERIES_SHARE[k] *
         Shape(x, e_capture - ReferenceLines::ESCAPE_XRAY_KEV[k], p);
  return v;
}

Double_t FixedEnergyGe(Double_t x, const Double_t *p) {
  return p[R_GE_FIXED_E_A] * Shape(x, ReferenceLines::GE73M_ENSDF, p);
}

Double_t NoGeModel(Double_t *x, Double_t *p) {
  return Bkg(x[0], p) + PbLines(x[0], p) +
         EscapeFamily(x[0], p, 0, ReferenceLines::N_ESCAPE_LINES - 1) +
         FixedEnergyGe(x[0], p);
}

Double_t DrawCd(Double_t x, const Double_t *p) {
  return Bkg(x, p) + EscapeFamily(x, p, 0, 3);
}
Double_t DrawTe(Double_t x, const Double_t *p) {
  return Bkg(x, p) + EscapeFamily(x, p, 4, 7);
}
Double_t DrawGeFixed(Double_t x, const Double_t *p) {
  return Bkg(x, p) + FixedEnergyGe(x, p);
}

struct EscapeResult {
  Double_t n96 = 0, chi2 = 0;
  Double_t y[4] = {0, 0, 0, 0},
           e[4] = {0, 0, 0, 0}; // Cd Ka, Cd Kb, Te Ka, Te Kb
  Double_t te_over_cdkb = 0, te_over_cdkb_err = 0;
  Double_t kb_over_ka = 0, kb_over_ka_err = 0;
  Double_t te_over_cdka = 0, te_over_cdka_err = 0;
  Double_t tekb_over_teka = 0, tekb_over_teka_err = 0;
  Double_t pb_ratio = 0, ge_a = 0, ge_a_err = 0;
  Double_t a = 0, b = 1;
  Int_t status = -1;
  Bool_t ok = kFALSE;
};

Double_t RatioError(Double_t a, Double_t b, Double_t vaa, Double_t vbb,
                    Double_t vab) {
  if (a == 0 || b == 0)
    return 0;
  Double_t r = a / b;
  Double_t rel2 = vaa / (a * a) + vbb / (b * b) - 2.0 * vab / (a * b);
  return TMath::Abs(r) * TMath::Sqrt(TMath::Max(rel2, 0.0));
}

// Staged as FitInRun: the scale fixed at the Pb Ka1 peak, then released, then
// the ratios and escapes released.
TF1 *FitNoGe(TH1D *h, const TString &tag, Double_t n96, EscapeResult &res) {
  TF1 *f = new TF1("noge_" + tag, NoGeModel, g_lo, g_hi, N_R);
  f->SetNpx(3000);
  h->GetXaxis()->SetRangeUser(74.2, 75.6);
  Double_t ka1_mode = h->GetBinCenter(h->GetMaximumBin());
  Double_t ka1_height = h->GetMaximum();
  h->GetXaxis()->SetRange(0, 0);
  Double_t level = h->GetBinContent(h->FindBin(65.0));
  Double_t a0 = ka1_mode - ReferenceLines::PB_KA1 + 0.3;
  Double_t b0 = 1.0;
  Double_t ka1 = TMath::Max((ka1_height - level) * 0.6 * 2.5 / BIN, 10.0);
  if (g_anchor_capture) {
    b0 = g_b_fixed;
    a0 = g_capture_meas - ReferenceLines::CD113_CAPTURE * b0;
    ka1 = 0.01 * n96;
  }

  f->SetParameter(Q_SIG0, 0.6);
  f->SetParLimits(Q_SIG0, 0.1, 1.5);
  f->SetParameter(Q_TAU, 1.0);
  f->SetParLimits(Q_TAU, 0.02, 5.0);
  f->SetParameter(Q_C0, level);
  f->SetParameter(Q_C1, 0);
  f->FixParameter(Q_GE_E, ReferenceLines::GE73M_ENSDF);
  f->FixParameter(Q_GE_A, 0.0);
  f->SetParameter(Q_KA1_A, ka1);
  f->SetParLimits(Q_KA1_A, 0, 1e9);
  f->SetParameter(Q_KB1_A, 0.25 * ka1);
  f->SetParLimits(Q_KB1_A, 0, 1e9);
  f->FixParameter(Q_A, a0);
  f->FixParameter(Q_B, b0);
  f->FixParameter(Q_SCALE2, 0.0);
  f->FixParameter(Q_KA2_R, PB_KA2_OVER_KA1_DDEP);
  f->FixParameter(Q_KB2_R, 0.3);
  f->FixParameter(Q_ESC_A, 0.0);
  f->FixParameter(R_CD_KA, 0.0);
  f->FixParameter(R_TE_KA, 0.0);
  f->FixParameter(R_TE_KB, 0.0);
  f->FixParameter(R_GE_FIXED_E_A, 0.0);
  f->FixParameter(Q_TAU_GE, 1.0);
  f->FixParameter(Q_SIG1, 0.0);
  h->Fit(f, "QLRN0");

  if (!g_anchor_capture) {
    f->ReleaseParameter(Q_A);
    f->ReleaseParameter(Q_B);
    f->SetParLimits(Q_A, a0 - 1.0, a0 + 1.0);
    f->SetParLimits(Q_B, 0.97, 1.03);
    h->Fit(f, "QLRN0");
  }

  if (!g_fix_pb_ratio) {
    f->ReleaseParameter(Q_KA2_R);
    f->SetParLimits(Q_KA2_R, 0.3, 1.2);
  }
  f->ReleaseParameter(Q_KB2_R);
  f->SetParLimits(Q_KB2_R, 0.1, 0.6);
  // Yields may go negative so a family absent from the data is not railed at
  // zero with a meaningless error.
  const Int_t family[4] = {R_CD_KA, Q_ESC_A, R_TE_KA, R_TE_KB};
  const Double_t start[4] = {0.05, 0.015, 0.003, 0.001};
  for (Int_t i = 0; i < 4; i++) {
    f->ReleaseParameter(family[i]);
    f->SetParameter(family[i], start[i] * n96);
    f->SetParLimits(family[i], -0.1 * n96, 0.5 * n96);
  }
  if (g_free_ge) {
    f->ReleaseParameter(R_GE_FIXED_E_A);
    f->SetParameter(R_GE_FIXED_E_A, 0.001 * n96);
    f->SetParLimits(R_GE_FIXED_E_A, -0.1 * n96, 0.5 * n96);
  }
  f->ReleaseParameter(Q_SIG1);
  f->SetParLimits(Q_SIG1, -0.05, 0.05);
  h->Fit(f, "QLRN0");
  TFitResultPtr fr = h->Fit(f, "SQLRN0M");

  res.status = fr.Get() ? (Int_t)fr->Status() : -1;
  res.ok = fr.Get() && fr->IsValid();
  res.a = f->GetParameter(Q_A);
  res.b = f->GetParameter(Q_B);
  res.chi2 = f->GetChisquare() / f->GetNDF();
  for (Int_t i = 0; i < 4; i++) {
    res.y[i] = f->GetParameter(family[i]);
    res.e[i] = f->GetParError(family[i]);
  }
  res.pb_ratio = f->GetParameter(Q_KA2_R);
  res.ge_a = f->GetParameter(R_GE_FIXED_E_A);
  res.ge_a_err = f->GetParError(R_GE_FIXED_E_A);
  if (fr.Get()) {
    res.te_over_cdkb = res.y[2] / res.y[1];
    res.te_over_cdkb_err = RatioError(
        res.y[2], res.y[1], fr->CovMatrix(R_TE_KA, R_TE_KA),
        fr->CovMatrix(Q_ESC_A, Q_ESC_A), fr->CovMatrix(R_TE_KA, Q_ESC_A));
    res.kb_over_ka = res.y[1] / res.y[0];
    res.kb_over_ka_err = RatioError(
        res.y[1], res.y[0], fr->CovMatrix(Q_ESC_A, Q_ESC_A),
        fr->CovMatrix(R_CD_KA, R_CD_KA), fr->CovMatrix(Q_ESC_A, R_CD_KA));
    res.te_over_cdka = res.y[2] / res.y[0];
    res.te_over_cdka_err = RatioError(
        res.y[2], res.y[0], fr->CovMatrix(R_TE_KA, R_TE_KA),
        fr->CovMatrix(R_CD_KA, R_CD_KA), fr->CovMatrix(R_TE_KA, R_CD_KA));
    res.tekb_over_teka = res.y[3] / res.y[2];
    res.tekb_over_teka_err = RatioError(
        res.y[3], res.y[2], fr->CovMatrix(R_TE_KB, R_TE_KB),
        fr->CovMatrix(R_TE_KA, R_TE_KA), fr->CovMatrix(R_TE_KB, R_TE_KA));
  }
  return f;
}

void DrawNoGe(TH1D *h, TF1 *f, const TString &name, const TString &label) {
  const Double_t *p = f->GetParameters();
  std::vector<TGraph *> components;
  components.push_back(Component(g_lo, g_hi, kGray + 2, 2, DrawBkg, p));
  components.push_back(Component(g_lo, g_hi, kOrange + 7, 1, DrawPb, p));
  components.push_back(Component(g_lo, g_hi, kBlue + 1, 1, DrawCd, p));
  components.push_back(Component(g_lo, g_hi, kMagenta + 1, 1, DrawTe, p));
  if (g_free_ge)
    components.push_back(Component(g_lo, g_hi, kGreen + 2, 1, DrawGeFixed, p));
  Plot(h, f, components, g_lo, g_hi, name, label);
}

void Print(const TString &label, const EscapeResult &r) {
  std::cout << "NOGE " << label << ": N96 " << std::fixed
            << std::setprecision(0) << r.n96
            << "  fractions of 95.9 [%]: Cd Ka " << std::setprecision(2)
            << 100 * r.y[0] / r.n96 << " +/- " << 100 * r.e[0] / r.n96
            << ", Cd Kb " << 100 * r.y[1] / r.n96 << " +/- "
            << 100 * r.e[1] / r.n96 << ", Te Ka " << 100 * r.y[2] / r.n96
            << " +/- " << 100 * r.e[2] / r.n96 << ", Te Kb "
            << 100 * r.y[3] / r.n96 << " +/- " << 100 * r.e[3] / r.n96
            << std::endl;
  std::cout << "NOGE " << label << ":   Te Ka/Cd Kb " << std::setprecision(3)
            << r.te_over_cdkb << " +/- " << r.te_over_cdkb_err << "  Cd Kb/Ka "
            << r.kb_over_ka << " +/- " << r.kb_over_ka_err << "  Te Ka/Cd Ka "
            << r.te_over_cdka << " +/- " << r.te_over_cdka_err << "  Te Kb/Ka "
            << r.tekb_over_teka << " +/- " << r.tekb_over_teka_err
            << "  Pb Ka2/Ka1 " << r.pb_ratio << "  chi2/ndf "
            << std::setprecision(2) << r.chi2 << "  scale a "
            << std::setprecision(4) << r.a << " b " << r.b << "  status "
            << r.status << (r.ok ? "" : " [invalid]");
  if (g_free_ge)
    std::cout << "  68.752 line " << std::setprecision(2)
              << 100 * r.ge_a / r.n96 << " +/- " << 100 * r.ge_a_err / r.n96
              << " % of 95.9";
  std::cout << std::endl;
}

EscapeResult RunOne(TH1D *h, const TString &tag, const TString &label,
                    Bool_t draw) {
  CaptureLine capture = FitCaptureLine(h, "noge_" + tag);
  g_capture_meas = capture.mu;
  EscapeResult r;
  r.n96 = capture.area;
  TF1 *f = FitNoGe(h, tag, capture.area, r);
  Print(label, r);
  if (draw)
    DrawNoGe(h, f, "psi_cdte_noge_escapes_" + tag, label);
  delete f;
  return r;
}

// The NoPb run against a Pb run, scaled to the same counts in 90-99 keV: what
// the Pb lines leave behind.
void DrawNoPbAgainstPb(const std::vector<TString> &nopb,
                       const std::vector<TString> &pb) {
  TH1D *h_fast = FillRise(nopb, "nopb_fast", 0, 280);
  TH1D *h_slower = FillRise(nopb, "nopb_280_300", 280, 300);
  TH1D *h_pb = FillRise(pb, "pb_fast", 0, 280);
  Double_t reference = h_pb->Integral(h_pb->FindBin(90.0), h_pb->FindBin(99.0));
  h_fast->Scale(reference /
                h_fast->Integral(h_fast->FindBin(90.0), h_fast->FindBin(99.0)));
  h_slower->Scale(reference / h_slower->Integral(h_slower->FindBin(90.0),
                                                 h_slower->FindBin(99.0)));
  h_fast->Rebin(4);
  h_slower->Rebin(4);
  h_pb->Rebin(4);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  canvas->SetLogy();
  h_fast->GetXaxis()->SetRangeUser(55, 100);
  h_fast->GetYaxis()->SetTitle("Counts / 200 eV");
  h_fast->SetLineColor(kBlue + 1);
  h_slower->SetLineColor(kRed + 1);
  h_pb->SetLineColor(kGray + 2);
  h_fast->Draw("hist");
  h_slower->Draw("hist same");
  h_pb->Draw("hist same");
  TLegend *leg = PlottingUtils::AddLegend(0.16, 0.68, 0.62, 0.9);
  leg->AddEntry(h_fast, "10/13 no Pb, rise < 280 ns", "l");
  leg->AddEntry(h_slower, "10/13 no Pb, rise 280-300 ns", "l");
  leg->AddEntry(h_pb, "10/10 Pb, rise < 280 ns", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(canvas, "psi_cdte_noge_nopb_vs_pb_look", "",
                            PlotSaveOptions::kLOG);
}

} // namespace CdTeNoGe

void NoGeEscapes() {
  using namespace CdTeNoGe;
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root +
                                    "/plots/studies/psi_cdte/noge_escapes",
                                project_root + "/root_files");
  TGaxis::SetMaxDigits(3);
  g_lo = 62.0;
  g_hi = 89.0;
  ResetKbGroups();

  std::vector<TString> runs_pb = NoGePbRuns();
  std::vector<TString> runs_1010(runs_pb.begin(), runs_pb.begin() + 2);
  std::vector<TString> runs_1011(runs_pb.begin() + 2, runs_pb.end());
  std::vector<TString> runs_nopb;
  runs_nopb.push_back("PSI_Beam_NoPb_CdTe_Nothing_Si_Nothing_20251013_0652_1");

  // The budget's rise slices, on the runs with Pb (the Ge runs' setup).
  Double_t b_sum = 0;
  Int_t b_n = 0;
  for (Int_t s = 0; s < N_RISE_SLICES; s++) {
    TH1D *h = FillRise(runs_pb, Form("noge_slice_%d", s), RISE_SLICE_LO[s],
                       RISE_SLICE_HI[s]);
    for (Int_t variant = 0; variant < 2; variant++) {
      g_fix_pb_ratio = (variant == 1);
      TString tag = TString::Format("pb_slice%.0f_%.0f", RISE_SLICE_LO[s],
                                    RISE_SLICE_HI[s]) +
                    (g_fix_pb_ratio ? "_pbfixed" : "");
      TString label = TString::Format("No Ge, rise %.0f-%.0f ns",
                                      RISE_SLICE_LO[s], RISE_SLICE_HI[s]) +
                      (g_fix_pb_ratio ? ", Pb Ka2/Ka1 0.595" : "");
      EscapeResult r = RunOne(h, tag, "SLICE " + label, kTRUE);
      if (variant == 0 && RISE_SLICE_HI[s] <= 300.0) {
        b_sum += r.b;
        b_n++;
      }
    }
    g_fix_pb_ratio = kFALSE;
    delete h;
  }
  g_b_fixed = b_sum / b_n;
  std::cout << "NOGE scale gain from the Pb runs' slices <= 300 ns: "
            << std::setprecision(5) << g_b_fixed << std::endl;

  // Coarse selections per run group; the NoPb run anchored at its capture
  // line.
  const Int_t n_groups = 4;
  const char *group_tags[n_groups] = {"1010", "1011", "pb", "1013nopb"};
  const char *group_names[n_groups] = {
      "10/10 NOGe", "10/11-12 CdTe NoGe, Si Ge", "10/10-12 with Pb",
      "10/13 NoPb, nothing"};
  std::vector<TString> group_runs[n_groups] = {runs_1010, runs_1011, runs_pb,
                                               runs_nopb};
  const Int_t n_sel = 3;
  const Double_t sel_hi[n_sel] = {280, 300, 450};
  for (Int_t g = 0; g < n_groups; g++) {
    g_anchor_capture = (g == 3);
    for (Int_t s = 0; s < n_sel; s++) {
      TH1D *h =
          FillRise(group_runs[g], Form("noge_h_%d_%d", g, s), 0.0, sel_hi[s]);
      TString base =
          TString::Format("%s, rise < %.0f ns", group_names[g], sel_hi[s]);
      TString tag = TString::Format("%s_rise%.0f", group_tags[g], sel_hi[s]);
      for (Int_t variant = 0; variant < 3; variant++) {
        g_fix_pb_ratio = (variant == 1);
        g_free_ge = (variant == 2);
        if (g_anchor_capture && g_fix_pb_ratio)
          continue;
        TString vtag = (variant == 0)   ? ""
                       : (variant == 1) ? "_pbfixed"
                                        : "_gefree";
        TString vname = (variant == 0)   ? ", Pb Ka2/Ka1 free"
                        : (variant == 1) ? ", Pb Ka2/Ka1 0.595"
                                         : ", plus a 68.752 line";
        Bool_t draw = (s == 1) || (g_anchor_capture && variant == 0);
        RunOne(h, tag + vtag, base + vname, draw);
      }
      g_fix_pb_ratio = kFALSE;
      g_free_ge = kFALSE;
      delete h;
    }
  }
  g_anchor_capture = kFALSE;

  std::vector<TString> pb_1010_run;
  pb_1010_run.push_back("PSI_Beam_NOGe_20251010_1106_3");
  DrawNoPbAgainstPb(runs_nopb, pb_1010_run);
}
