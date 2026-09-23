#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <utility>
#include <vector>

const Float_t E_CD558 = 558.456;

const Float_t FIT_INIT_LO = 525.0;
const Float_t FIT_INIT_HI = 590.0;
const Float_t SIGMA_INIT_KEV = 2.0;

const Int_t MIN_FIT_EVENTS = 100;
const Int_t BAD_MIN_EVENTS = 1000;
const Float_t SIGMA_HARD_MAX_KEV = 5.0;
const Float_t SIGMA_MAD_CUT = 4.5;
const Float_t MU_ERR_MAX_KEV = 1.0;

const TString RUN_NAME =
    Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116;
const TString CAL_FILE = "pixel_calibration.root";

const std::set<std::pair<Int_t, Int_t>> INTERACTIVE_PIXELS = {
    // {ix, iy} entries here override the default auto fit for that pixel.
};

const Float_t HALFWIDTH = Constants::PIXEL_ACCEPT_HALFWIDTH_MM;

struct PixelEvents {
  std::vector<Float_t> x;
  std::vector<Float_t> y;
  std::vector<Float_t> e;
};

struct PixelRow {
  Int_t crystal = -1;
  Int_t ix = -1;
  Int_t iy = -1;
  Int_t n_events = 0;
  Float_t mu = 0;
  Float_t mu_err = 0;
  Float_t sigma = 0;
  Bool_t valid = kFALSE;
  Bool_t bad = kFALSE;
  TString reason;
};

Int_t GetCrystalIndex(Float_t x, Float_t y) {
  if (x < 0 && y < 0)
    return 0;
  if (x > 0 && y < 0)
    return 1;
  if (x < 0 && y > 0)
    return 2;
  if (x > 0 && y > 0)
    return 3;
  return -1;
}

PixelEvents LoadBefEvents(const TString &run_name) {
  PixelEvents out;
  TFile *file = IO::OpenForReading("filtered/" + run_name + ".root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open filtered/" << run_name << ".root"
              << std::endl;
    return out;
  }
  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  if (!tree) {
    std::cerr << "ERROR: bef_tree not found in filtered/" << run_name << ".root"
              << std::endl;
    file->Close();
    delete file;
    return out;
  }

  Float_t energy = 0, x = 0, y = 0;
  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);

  Int_t n = tree->GetEntries();
  out.x.reserve(n);
  out.y.reserve(n);
  out.e.reserve(n);
  for (Int_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    out.x.push_back(x);
    out.y.push_back(y);
    out.e.push_back(energy);
  }
  file->Close();
  delete file;
  return out;
}

TH1F *BuildPixelHist(const std::vector<Float_t> &energies,
                     const TString &name) {
  TH1F *hist = new TH1F(
      name, Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
      Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
  hist->SetDirectory(0);
  for (size_t i = 0; i < energies.size(); i++)
    hist->Fill(energies[i]);
  return hist;
}

static Float_t Median(std::vector<Float_t> v) {
  if (v.empty())
    return std::nanf("");
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

static Float_t Mad(std::vector<Float_t> v) {
  if (v.empty())
    return std::nanf("");
  std::sort(v.begin(), v.end());
  Float_t med = v[v.size() / 2];
  for (size_t i = 0; i < v.size(); i++)
    v[i] = TMath::Abs(v[i] - med);
  std::sort(v.begin(), v.end());
  return 1.4826f * v[v.size() / 2];
}

void ClassifyBadPixels(std::vector<PixelRow> &rows) {
  std::vector<std::vector<TString>> reasons(rows.size());

  for (size_t i = 0; i < rows.size(); i++) {
    PixelRow &r = rows[i];
    if (r.n_events < BAD_MIN_EVENTS)
      reasons[i].push_back(Form("low_stats_%d", r.n_events));
    if (!r.valid)
      reasons[i].push_back("invalid_fit");
    if (std::isfinite(r.sigma) && r.sigma > SIGMA_HARD_MAX_KEV)
      reasons[i].push_back(
          Form("sigma=%.2f>%.2f", r.sigma, SIGMA_HARD_MAX_KEV));
    if (std::isfinite(r.mu_err) && r.mu_err > MU_ERR_MAX_KEV)
      reasons[i].push_back(Form("mu_err=%.3f>%.2f", r.mu_err, MU_ERR_MAX_KEV));
  }

  std::map<Int_t, std::vector<Float_t>> sig_by_crystal;
  std::map<Int_t, std::vector<size_t>> idx_by_crystal;
  for (size_t i = 0; i < rows.size(); i++) {
    const PixelRow &r = rows[i];
    if (r.valid && std::isfinite(r.sigma) && r.sigma <= SIGMA_HARD_MAX_KEV) {
      sig_by_crystal[r.crystal].push_back(r.sigma);
      idx_by_crystal[r.crystal].push_back(i);
    }
  }
  for (std::map<Int_t, std::vector<Float_t>>::iterator it =
           sig_by_crystal.begin();
       it != sig_by_crystal.end(); ++it) {
    if (it->second.size() < 5)
      continue;
    Float_t med = Median(it->second);
    Float_t m = Mad(it->second);
    if (m == 0.0f || !std::isfinite(m))
      continue;
    Float_t cut = SIGMA_MAD_CUT * m;
    const std::vector<size_t> &idxs = idx_by_crystal[it->first];
    for (size_t k = 0; k < idxs.size(); k++) {
      size_t i = idxs[k];
      Float_t s = rows[i].sigma;
      if (TMath::Abs(s - med) > cut)
        reasons[i].push_back(
            Form("sigma_outlier_%.3f_vs_%.3f_mad_%.3f", s, med, m));
    }
  }

  for (size_t i = 0; i < rows.size(); i++) {
    if (reasons[i].empty())
      continue;
    rows[i].bad = kTRUE;
    TString joined;
    for (size_t k = 0; k < reasons[i].size(); k++) {
      if (k)
        joined += "|";
      joined += reasons[i][k];
    }
    rows[i].reason = joined;
  }
}

std::vector<Double_t> PixelEdges(const std::vector<Float_t> &centers) {
  Int_t n = (Int_t)centers.size();
  std::vector<Double_t> edges(n + 1);
  for (Int_t i = 1; i < n; i++)
    edges[i] = 0.5 * (centers[i - 1] + centers[i]);
  edges[0] = centers[0] - 0.5 * (centers[1] - centers[0]);
  edges[n] = centers[n - 1] + 0.5 * (centers[n - 1] - centers[n - 2]);
  return edges;
}

TH2D *BuildGainHist(const std::vector<PixelRow> &rows) {
  std::vector<Double_t> x_edges = PixelEdges(Constants::PIXEL_CENTERS_X_MM);
  std::vector<Double_t> y_edges = PixelEdges(Constants::PIXEL_CENTERS_Y_MM);
  Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();

  TH2D *gain = new TH2D(
      "pixel_gain", Form("; x [mm]; y [mm]; gain (%.3f keV / #mu)", E_CD558),
      nx, x_edges.data(), ny, y_edges.data());
  gain->SetDirectory(0);

  for (size_t i = 0; i < rows.size(); i++) {
    const PixelRow &r = rows[i];
    if (r.bad || !std::isfinite(r.mu) || r.mu <= 0)
      continue;
    gain->SetBinContent(r.ix + 1, r.iy + 1, E_CD558 / r.mu);
  }
  return gain;
}

TH2D *BuildMuHist(const std::vector<PixelRow> &rows) {
  std::vector<Double_t> x_edges = PixelEdges(Constants::PIXEL_CENTERS_X_MM);
  std::vector<Double_t> y_edges = PixelEdges(Constants::PIXEL_CENTERS_Y_MM);
  Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();

  TH2D *mu =
      new TH2D("pixel_mu", Form("; x [mm]; y [mm]; #mu [keV] (%.3f)", E_CD558),
               nx, x_edges.data(), ny, y_edges.data());
  mu->SetDirectory(0);
  for (size_t i = 0; i < rows.size(); i++) {
    const PixelRow &r = rows[i];
    if (!r.valid)
      continue;
    mu->SetBinContent(r.ix + 1, r.iy + 1, r.mu);
  }
  return mu;
}

void PlotMuHeatmap(TH2D *mu_hist, const std::vector<PixelRow> &rows) {
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(mu_hist, canvas);
  Float_t mu_min = 1e9, mu_max = -1e9;
  for (size_t i = 0; i < rows.size(); i++) {
    if (!rows[i].valid)
      continue;
    if (rows[i].mu < mu_min)
      mu_min = rows[i].mu;
    if (rows[i].mu > mu_max)
      mu_max = rows[i].mu;
  }
  if (mu_max > mu_min) {
    Float_t pad = 0.05 * (mu_max - mu_min);
    mu_hist->GetZaxis()->SetRangeUser(mu_min - pad, mu_max + pad);
  }
  PlottingUtils::SaveFigure(canvas, "pixel_mu_cd558", "calibration",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

void PrintSummary(const std::vector<PixelRow> &rows) {
  Int_t n_total = (Int_t)rows.size();
  Int_t n_valid = 0, n_bad = 0;
  for (size_t i = 0; i < rows.size(); i++) {
    if (rows[i].valid)
      n_valid++;
    if (rows[i].bad)
      n_bad++;
  }
  std::cout << "total pixels: " << n_total << std::endl;
  std::cout << "fits valid: " << n_valid << std::endl;
  std::cout << "bad pixels: " << n_bad << " / " << n_total << " (" << std::fixed
            << std::setprecision(1) << (100.0 * n_bad / n_total) << "%)"
            << std::endl;

  std::map<TString, Int_t> reason_counts;
  for (size_t i = 0; i < rows.size(); i++) {
    if (!rows[i].bad)
      continue;
    TObjArray *toks = rows[i].reason.Tokenize("|");
    for (Int_t t = 0; t < toks->GetEntries(); t++) {
      TString tag = ((TObjString *)toks->At(t))->GetString();
      TString key;
      if (tag.BeginsWith("sigma_outlier"))
        key = "sigma_outlier";
      else if (tag.BeginsWith("sigma"))
        key = "sigma_hard_cut";
      else if (tag.BeginsWith("mu_err"))
        key = "mu_err";
      else if (tag.BeginsWith("low_stats"))
        key = "low_stats";
      else if (tag.BeginsWith("invalid_fit"))
        key = "invalid_fit";
      else
        key = tag;
      reason_counts[key]++;
    }
    delete toks;
  }
  if (!reason_counts.empty()) {
    std::cout << "Breakdown (a pixel may have multiple reasons):" << std::endl;
    std::vector<std::pair<Int_t, TString>> sorted;
    for (std::map<TString, Int_t>::iterator it = reason_counts.begin();
         it != reason_counts.end(); ++it)
      sorted.push_back(std::make_pair(it->second, it->first));
    std::sort(
        sorted.begin(), sorted.end(),
        [](const std::pair<Int_t, TString> &a,
           const std::pair<Int_t, TString> &b) { return a.first > b.first; });
    for (size_t i = 0; i < sorted.size(); i++)
      std::cout << "  " << std::left << std::setw(20) << sorted[i].second
                << ": " << sorted[i].first << std::endl;
  }

  std::vector<Float_t> good_gains;
  for (size_t i = 0; i < rows.size(); i++) {
    const PixelRow &r = rows[i];
    if (r.bad || !std::isfinite(r.mu) || r.mu <= 0)
      continue;
    good_gains.push_back(E_CD558 / r.mu);
  }
  if (!good_gains.empty()) {
    Float_t gmin = *std::min_element(good_gains.begin(), good_gains.end());
    Float_t gmax = *std::max_element(good_gains.begin(), good_gains.end());
    Double_t mean = 0;
    for (size_t i = 0; i < good_gains.size(); i++)
      mean += good_gains[i];
    mean /= good_gains.size();
    Double_t var = 0;
    for (size_t i = 0; i < good_gains.size(); i++)
      var += (good_gains[i] - mean) * (good_gains[i] - mean);
    var /= (good_gains.size() - 1);
    std::cout << "Per-pixel gain across " << good_gains.size()
              << " good pixels:" << std::endl;
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  median = " << Median(good_gains) << std::endl;
    std::cout << "  min / max = " << gmin << " / " << gmax << std::endl;
    std::cout << "  spread (stddev) = " << TMath::Sqrt(var) << std::endl;
  }

  std::cout << "Per-crystal bad count:" << std::endl;
  std::map<Int_t, Int_t> bad_per_crystal, total_per_crystal;
  for (size_t i = 0; i < rows.size(); i++) {
    total_per_crystal[rows[i].crystal]++;
    if (rows[i].bad)
      bad_per_crystal[rows[i].crystal]++;
  }
  for (std::map<Int_t, Int_t>::iterator it = total_per_crystal.begin();
       it != total_per_crystal.end(); ++it) {
    std::cout << "  crystal " << it->first << ": " << bad_per_crystal[it->first]
              << " / " << it->second << std::endl;
  }
}

void PixelCalibration() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots/pixel_calibration",
                                project_root + "/root_files");

  std::cout << "Loading BEF events from " << RUN_NAME << "..." << std::endl;
  PixelEvents ev = LoadBefEvents(RUN_NAME);
  std::cout << "  " << ev.e.size() << " interactions loaded" << std::endl;
  if (ev.e.empty())
    return;

  std::vector<PixelRow> rows;
  Int_t n_fits = 0, n_low = 0, n_failed = 0;

  Int_t nx = (Int_t)Constants::PIXEL_CENTERS_X_MM.size();
  Int_t ny = (Int_t)Constants::PIXEL_CENTERS_Y_MM.size();

  for (Int_t ix = 0; ix < nx; ix++) {
    Float_t x_c = Constants::PIXEL_CENTERS_X_MM[ix];
    for (Int_t iy = 0; iy < ny; iy++) {
      Float_t y_c = Constants::PIXEL_CENTERS_Y_MM[iy];

      std::vector<Float_t> pix_e;
      pix_e.reserve(1024);
      Int_t n_in_window = 0;
      for (size_t k = 0; k < ev.e.size(); k++) {
        if (TMath::Abs(ev.x[k] - x_c) < HALFWIDTH &&
            TMath::Abs(ev.y[k] - y_c) < HALFWIDTH) {
          pix_e.push_back(ev.e[k]);
          if (ev.e[k] >= FIT_INIT_LO && ev.e[k] <= FIT_INIT_HI)
            n_in_window++;
        }
      }

      PixelRow row;
      row.crystal = GetCrystalIndex(x_c, y_c);
      row.ix = ix;
      row.iy = iy;
      row.n_events = n_in_window;

      if (row.n_events < MIN_FIT_EVENTS) {
        n_low++;
        rows.push_back(row);
        continue;
      }

      n_fits++;
      TString peak_name = Form("Cd558_ix%02d_iy%02d", ix, iy);
      std::cout << "Pixel (" << ix << "," << iy << ") x=" << x_c << " y=" << y_c
                << " N=" << row.n_events << " " << peak_name << std::endl;

      TH1F *hist = BuildPixelHist(pix_e, "h_" + peak_name);

      FittingUtils *fitter = new FittingUtils(hist, FIT_INIT_LO, FIT_INIT_HI,
                                              /*use_flat_background=*/kTRUE,
                                              /*use_step=*/kFALSE,
                                              /*use_low_exp_tail=*/kTRUE,
                                              /*use_low_lin_tail=*/kTRUE,
                                              /*use_high_exp_tail=*/kFALSE);
      fitter->GetFitFunction()->SetParameter(0, E_CD558);
      fitter->GetFitFunction()->SetParameter(1, SIGMA_INIT_KEV);
      if (INTERACTIVE_PIXELS.count(std::make_pair(ix, iy)))
        fitter->SetInteractive();

      FitResult result;
      try {
        result = fitter->FitSinglePeak(RUN_NAME, peak_name);
      } catch (...) {
        n_failed++;
        delete fitter;
        delete hist;
        rows.push_back(row);
        continue;
      }

      if (result.valid && !result.peaks.empty()) {
        const PeakFitResult &p = result.peaks.at(0);
        row.mu = p.mu;
        row.mu_err = p.mu_error;
        row.sigma = p.sigma;
        row.valid = kTRUE;
      }

      delete fitter;
      delete hist;
      rows.push_back(row);
    }
  }

  std::cout << "fits attempted: " << n_fits << std::endl;
  std::cout << "fits failed (exception): " << n_failed << std::endl;
  std::cout << "pixels with low stats: " << n_low << std::endl;

  ClassifyBadPixels(rows);

  TH2D *gain = BuildGainHist(rows);
  TH2D *mu = BuildMuHist(rows);

  TFile *out = IO::OpenForWriting(CAL_FILE);
  gain->Write("pixel_gain", TObject::kOverwrite);
  mu->Write("pixel_mu", TObject::kOverwrite);
  out->Close();
  delete out;
  std::cout << "Wrote root_files/" << CAL_FILE << " (pixel_gain, pixel_mu)"
            << std::endl;

  PlotMuHeatmap(mu, rows);

  PrintSummary(rows);

  delete gain;
  delete mu;
}
