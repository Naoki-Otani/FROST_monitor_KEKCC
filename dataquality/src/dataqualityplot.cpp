#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TError.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TDatime.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TPad.h>

#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include "/home/nu/notani/FROST_monitor/config/config.hpp"

// ------------------------
// Global configuration
// ------------------------

static const char* TREE_NAME = "tree";

// Timing hist settings
static const Int_t    NBINS_TDC = FrostmonConfig::NBINS_TDC;
static const Double_t XMIN_TDC  = FrostmonConfig::XMIN_TDC;
static const Double_t XMAX_TDC  = FrostmonConfig::XMAX_TDC;

static const Int_t    NBINS_ADC = FrostmonConfig::NBINS_ADC;
static const Double_t XMIN_ADC  = FrostmonConfig::XMIN_ADC;
static const Double_t XMAX_ADC  = FrostmonConfig::XMAX_ADC;

// Lightyield settings
static const Int_t    NOUT       = FrostmonConfig::NOUT;
static const Int_t    NBUNCH     = FrostmonConfig::NBUNCH;

static const Int_t    NBINS_LYAVG = FrostmonConfig::NBINS_LYAVG;
static const Double_t XMIN_LYAVG  = FrostmonConfig::XMIN_LYAVG;
static const Double_t XMAX_LYAVG  = FrostmonConfig::XMAX_LYAVG;

// xg–yg settings
static const Int_t    XBINS_XY  = FrostmonConfig::XBINS_XY;
static const Int_t    YBINS_XY  = FrostmonConfig::YBINS_XY;
static const Double_t XMIN_XY   = FrostmonConfig::XMIN_XY;
static const Double_t XMAX_XY   = FrostmonConfig::XMAX_XY;
static const Double_t YMIN_XY   = FrostmonConfig::YMIN_XY;
static const Double_t YMAX_XY   = FrostmonConfig::YMAX_XY;

static const Double_t XG_WEIGHT = FrostmonConfig::XG_WEIGHT;        // weight exponent
static const Double_t LIGHTMAX_MIN = FrostmonConfig::LIGHTMAX_MIN;    // threshold for xg–yg selection

// 6-hour binning for LY history
static const std::string PATH_LY_PROC  = FrostmonConfig::OUTPUT_DIR + "/dataquality/lightyield/processed_files.tsv";
static const std::string PATH_LY_BINS  = FrostmonConfig::OUTPUT_DIR + "/dataquality/lightyield/chavg_bins.tsv";
static const Double_t BINW_SEC = FrostmonConfig::BINW_SEC;
static const Int_t MIN_COUNTS_PER_BIN = FrostmonConfig::MIN_COUNTS_PER_BIN; // Minimum total number of ly>=10 p.e. entries per 6-hour bin to be shown

// I/O paths
static const std::string PATH_ROOT     = FrostmonConfig::OUTPUT_DIR + "/rootfile_aftercalib/";
static const std::string PATH_OUT_TDC  = FrostmonConfig::OUTPUT_DIR + "/dataquality/tdc/";
static const std::string PATH_OUT_LY   = FrostmonConfig::OUTPUT_DIR + "/dataquality/lightyield/";
static const std::string PATH_OUT_XY   = FrostmonConfig::OUTPUT_DIR + "/dataquality/xgyg/";
static const std::string PATH_OUT_UNIXTIME = FrostmonConfig::OUTPUT_DIR + "/dataquality/unixtime/";
static const std::string PATH_OUT_SPILLNUM = FrostmonConfig::OUTPUT_DIR + "/dataquality/spillnum/";
static const std::string PATH_OUT_LY_EACHCH = FrostmonConfig::OUTPUT_DIR + "/dataquality/lightyield_each_ch/";

// ------------------------
// Utilities
// ------------------------

// Convert a local date-time to Unix time (seconds since 1970-01-01).
static double unixTimeLocal(int year, int month, int day,
                            int hour, int min, int sec)
{
  std::tm tm{};
  tm.tm_year = year  - 1900;  // years since 1900
  tm.tm_mon  = month - 1;     // 0–11
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min  = min;
  tm.tm_sec  = sec;
  tm.tm_isdst = -1;           // let libc determine DST
  return static_cast<double>(std::mktime(&tm));
}

struct FileInfoLite {
  std::string path;
  Long_t size = 0;
  Long_t mtime = 0;
};

static bool GetPathInfo(const std::string& p, Long_t& id, Long_t& size, Long_t& flags, Long_t& mtime) {
  return gSystem->GetPathInfo(p.c_str(), &id, &size, &flags, &mtime) == 0;
}

static bool IsLikelyStableFile(const std::string& path, int wait_ms = 500) {
  Long_t id1=0, sz1=0, fl1=0, mt1=0;
  if (!GetPathInfo(path, id1, sz1, fl1, mt1)) return false;
  gSystem->Sleep(wait_ms);
  Long_t id2=0, sz2=0, fl2=0, mt2=0;
  if (!GetPathInfo(path, id2, sz2, fl2, mt2)) return false;
  return (sz1 == sz2);
}

static std::string ExtractRunTag(const std::string& fname) {
  // e.g. .../run00103_0_9999_lightyield.root -> run00103
  std::regex re("^((?:.*/)?)(run\\d{5})_\\d+_\\d+_lightyield\\.root$");
  std::smatch m;
  if (std::regex_match(fname, m, re)) return m[2].str();
  return "";
}

static std::vector<FileInfoLite> ListLightyieldFiles(const std::string& dir) {
  std::vector<FileInfoLite> out;
  TSystemDirectory d("rootdir", dir.c_str());
  TList* files = d.GetListOfFiles();
  if (!files) return out;

  TIter it(files);
  while (TSystemFile* f = (TSystemFile*)it()) {
    TString name = f->GetName();
    if (f->IsDirectory()) continue;
    if (!name.EndsWith("_lightyield.root")) continue;

    std::string path = dir + name.Data();
    Long_t id=0, sz=0, fl=0, mt=0;
    if (GetPathInfo(path, id, sz, fl, mt)) out.push_back({path, sz, mt});
  }
  std::sort(out.begin(), out.end(), [](const FileInfoLite& a, const FileInfoLite& b){
    return a.mtime > b.mtime;
  });
  return out;
}

static std::vector<FileInfoLite> FilterByRun(const std::vector<FileInfoLite>& v, const std::string& runTag) {
  std::vector<FileInfoLite> out;
  const std::string needle = "/" + runTag + "_";
  for (const auto& fi : v) if (fi.path.find(needle) != std::string::npos) out.push_back(fi);
  return out;
}

// Lightyield ROOT file readiness condition:
//   - size >= 10 KB
//   - file size stable (short check via IsLikelyStableFile)
//   - last modification time older than 60 seconds
static bool IsReadyLightyieldFile(const FileInfoLite& fi,
                                  Long_t minSizeBytes = 10L * 1024,
                                  Long_t stableSec    = 60)
{
  // (1) minimum size
  if (fi.size < minSizeBytes) return false;

  // (2) quick stability check (~0.5 s)
  if (!IsLikelyStableFile(fi.path, 500)) return false;

  // (3) mtime older than stableSec
  const Long_t nowSec = static_cast<Long_t>(time(nullptr));
  if (nowSec - fi.mtime < stableSec) return false;

  return true;
}

// ------------------------
// Timing histograms
// ------------------------

static TH1D* CreateHistTiming(const char* hname, bool isTDC, const char* titlePrefix, const char* xname) {
  if (isTDC) {
    TH1D* h = new TH1D(hname, Form("%s;%s [0.833 ns];Number of events", titlePrefix, xname), NBINS_TDC, XMIN_TDC, XMAX_TDC);
    h->SetLineWidth(2);
    return h;
  } else {
    TH1D* h = new TH1D(hname, Form("%s;%s [13.3 ns];Number of events", titlePrefix, xname), NBINS_ADC, XMIN_ADC, XMAX_ADC);
    h->SetLineWidth(2);
    return h;
  }
}

static void FillHistsTimingFromFile(TH1D* h_lead, TH1D* h_trail, TH1D* h_lead_adc, TH1D* h_trail_adc,
                                    const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* tin = (TTree*)fin.Get(TREE_NAME);
  if (!tin) { fin.Close(); return; }

  tin->SetBranchStatus("*", 0);

  std::vector<std::vector<double>> *leading = nullptr;
  std::vector<std::vector<double>> *trailing = nullptr;
  std::vector<std::vector<double>> *leading_fromadc = nullptr;
  std::vector<std::vector<double>> *trailing_fromadc = nullptr;

  if (tin->GetBranch("leading")) { tin->SetBranchStatus("leading", 1); tin->SetBranchAddress("leading", &leading); }
  if (tin->GetBranch("trailing")) { tin->SetBranchStatus("trailing", 1); tin->SetBranchAddress("trailing", &trailing); }
  if (tin->GetBranch("leading_fromadc")) { tin->SetBranchStatus("leading_fromadc", 1); tin->SetBranchAddress("leading_fromadc", &leading_fromadc); }
  if (tin->GetBranch("trailing_fromadc")) { tin->SetBranchStatus("trailing_fromadc", 1); tin->SetBranchAddress("trailing_fromadc", &trailing_fromadc); }

  const Long64_t nent = tin->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tin->GetEntry(ie);

    if (leading && h_lead) {
      const int nch = (int)leading->size();
      for (int ch = 0; ch < nch; ++ch) for (const double x : leading->at(ch))  h_lead->Fill(x);
    }
    if (trailing && h_trail) {
      const int nch = (int)trailing->size();
      for (int ch = 0; ch < nch; ++ch) for (const double x : trailing->at(ch)) h_trail->Fill(x);
    }
    if (leading_fromadc && h_lead_adc) {
      const int nch = (int)leading_fromadc->size();
      for (int ch = 0; ch < nch; ++ch) for (const double x : leading_fromadc->at(ch))  h_lead_adc->Fill(x);
    }
    if (trailing_fromadc && h_trail_adc) {
      const int nch = (int)trailing_fromadc->size();
      for (int ch = 0; ch < nch; ++ch) for (const double x : trailing_fromadc->at(ch)) h_trail_adc->Fill(x);
    }
  }
  fin.Close();
}

static void DrawAndSaveTiming(TH1D* h, const std::string& runTag, const char* nameSuffix) {
  TString outpdf = TString::Format("%s%s_%s_hist.pdf", PATH_OUT_TDC.c_str(), runTag.c_str(), nameSuffix);

  TCanvas* c = (TCanvas*)gROOT->FindObject("c_realtime");
  if (!c) c = new TCanvas("c_realtime", "timing hist (realtime)", 1000, 700);
  c->Clear();
  c->SetGrid();
  c->SetLogy();

  gStyle->SetOptStat(0);
  h->SetStats(0);
  h->SetMinimum(0.5);
  h->Draw("HIST");

  TDatime now;
  TPaveText note(0.75, 0.85, 0.99, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("Run: %s", runTag.c_str()));
  note.AddText(Form("Entries: %.0f", h->GetEntries()));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  c->SaveAs(outpdf);
}

// ------------------------
// Average LY per channel
// ------------------------

static void AccumulateLyAvgPerChannel(std::vector<double>& sum, std::vector<int>& cnt, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  Double_t ly[NOUT][NBUNCH];
  t->SetBranchStatus("lightyield", 1);
  t->SetBranchAddress("lightyield", ly);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    for (int ch = 0; ch < NOUT; ++ch) {
      for (int b = 0; b < NBUNCH; ++b) {
        const double v = ly[ch][b];
        if (std::isfinite(v) && v >= 10.0) { sum[ch] += v; cnt[ch] += 1; }
      }
    }
  }
  fin.Close();
}

static void BuildAndSaveLyAvgHist(const std::vector<double>& sum, const std::vector<int>& cnt, const std::string& runTag) {
  std::unique_ptr<TH1D> h(new TH1D("h_chavg_ly", "Average light yield per channel (>=10 p.e.);Light yield [p.e.];Number of channels",
                                   NBINS_LYAVG, XMIN_LYAVG, XMAX_LYAVG));
  h->SetLineWidth(2);
  // h->SetStats(1);
  for (int ch = 0; ch < NOUT; ++ch) {
    const double avg = (cnt[ch] > 0) ? (sum[ch] / (double)cnt[ch]) : 0.0;
    h->Fill(avg);
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  TCanvas c("c_chavg_ly", "chavg ly", 1000, 700);
  c.SetGrid();
  h->SetMinimum(0.5);

  // gStyle->SetOptStat(1110);
  // gStyle->SetStatW(0.22);
  // gStyle->SetStatH(0.16);
  // gStyle->SetStatY(0.78);
  h->Draw("HIST");

  TDatime now;
  TPaveText note(0.75, 0.72, 0.99, 0.9, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("Run: %s", runTag.c_str()));
  note.AddText("Entries: 272");
  note.AddText(Form("Mean: %.2f", h->GetMean()));
  note.AddText(Form("Std Dev: %.2f", h->GetStdDev()));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_chavg_lightyield_hist.pdf", PATH_OUT_LY.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// Lightyield hist per channel (RAYRAW planes)
// ------------------------

static void BuildAndSaveLyEachChannel(const std::string& runTag,
                                      const std::vector<FileInfoLite>& stableFiles)
{
  // Ensure histograms are not owned by any TFile
  const Bool_t oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // Configuration: RAYRAW planes and channels per plane
  static const int NPLANES = FrostmonConfig::N_RAYRAW;   // RAYRAW = 1..11
  static const int NCHPL   = FrostmonConfig::N_CH_PER_PLANE;   // channels per plane

  if (stableFiles.empty()) {
    TH1::AddDirectory(oldAddDir);
    return;
  }

  gSystem->mkdir(PATH_OUT_LY_EACHCH.c_str(), true);

  // (plane, local ch) -> cablenum mapping
  int cable_map[NPLANES][NCHPL];
  for (int p = 0; p < NPLANES; ++p) {
    for (int lc = 0; lc < NCHPL; ++lc) {
      cable_map[p][lc] = -1;
    }
  }

  bool mapping_built = false;

  // Histograms: [plane][local channel]
  std::vector<std::vector<TH1D*>> hists(NPLANES, std::vector<TH1D*>(NCHPL, (TH1D*)nullptr));

  const int nbins = FrostmonConfig::NBINS_LYCH;
  const double xmin = FrostmonConfig::XMIN_LYCH;
  const double xmax = FrostmonConfig::XMAX_LYCH;

  // Loop over all ready ROOT files of the latest run
  for (const auto& fi : stableFiles) {
    TFile f(fi.path.c_str(), "READ");
    if (f.IsZombie()) continue;
    TTree* t = (TTree*)f.Get(TREE_NAME);
    if (!t) { f.Close(); continue; }

    t->SetBranchStatus("*", 0);

    // Branch arrays
    Int_t   cablenum[NOUT];
    Int_t   rayraw[NOUT];
    Int_t   rayrawch[NOUT];
    Double_t ly[NOUT][NBUNCH];

    // Require all needed branches
    if (!t->GetBranch("cablenum") ||
        !t->GetBranch("rayraw")   ||
        !t->GetBranch("rayrawch") ||
        !t->GetBranch("lightyield")) {
      f.Close();
      continue;
    }

    t->SetBranchStatus("cablenum",   1);
    t->SetBranchStatus("rayraw",     1);
    t->SetBranchStatus("rayrawch",   1);
    t->SetBranchStatus("lightyield", 1);
    t->SetBranchAddress("cablenum",   cablenum);
    t->SetBranchAddress("rayraw",     rayraw);
    t->SetBranchAddress("rayrawch",   rayrawch);
    t->SetBranchAddress("lightyield", ly);

    const Long64_t nent = t->GetEntries();
    if (nent <= 0) { f.Close(); continue; }

    // Build mapping from (plane, local ch) to cablenum only once
    if (!mapping_built) {
      t->GetEntry(0);
      for (int i = 0; i < NOUT; ++i) {
        const int rr = rayraw[i];   // expected 1..NPLANES
        const int lc = rayrawch[i]; // expected 0..NCHPL-1
        if (rr < 1 || rr > NPLANES) continue;
        if (lc < 0 || lc >= NCHPL)  continue;
        cable_map[rr-1][lc] = cablenum[i];
      }

      // Create histograms for all (plane, ch) in the global directory
      gROOT->cd();  // <-- add this line

      // Create histograms for all (plane, ch)
      for (int p = 0; p < NPLANES; ++p) {
        for (int lc = 0; lc < NCHPL; ++lc) {
          TString hname  = Form("h_rayraw%02d_ch%02d", p+1, lc);
          const int cable = cable_map[p][lc];

          TString htitle;
          if (cable > 0) {
            htitle = Form("RAYRAW#%d ch%02d (cable %03d);Light yield [p.e.];Number of events",
                          p+1, lc, cable);
          } else {
            htitle = Form("RAYRAW#%d ch%02d;Light yield [p.e.];Number of events",
                          p+1, lc);
          }

          TH1D* h = new TH1D(hname, htitle, nbins, xmin, xmax);
          h->SetLineWidth(1);
          h->SetStats(1);
          hists[p][lc] = h;
        }
      }

      mapping_built = true;
    }

    // Fill histograms for all events in this file
    for (Long64_t ie = 0; ie < nent; ++ie) {
      t->GetEntry(ie);

      for (int i = 0; i < NOUT; ++i) {
        const int rr = rayraw[i];   // 1..NPLANES
        const int lc = rayrawch[i]; // 0..NCHPL-1
        if (rr < 1 || rr > NPLANES) continue;
        if (lc < 0 || lc >= NCHPL)  continue;

        TH1D* h = hists[rr-1][lc];
        if (!h) continue;

        // Fill from all bunches for this channel
        for (int b = 0; b < NBUNCH; ++b) {
          const double v = ly[i][b];
          if (!std::isfinite(v)) continue;
          h->Fill(v);
        }
      }
    }

    f.Close();
  }

  // If mapping was never built (no suitable file/branch), nothing to do
  if (!mapping_built) return;

  // Drawing and saving: one PDF per RAYRAW plane, 8x4 pads, log-y
  gStyle->SetOptStat(1110);
  gStyle->SetStatW(0.22);
  gStyle->SetStatH(0.16);

  for (int p = 0; p < NPLANES; ++p) {
    TCanvas* c = new TCanvas(Form("c_rayraw%02d", p+1),
                             Form("RAYRAW#%d light yield", p+1),
                             2000, 1400);
    c->Divide(8, 4);

    for (int lc = 0; lc < NCHPL; ++lc) {
      c->cd(lc + 1);
      gPad->SetLogy(1);
      gPad->SetLeftMargin(0.15);

      TH1D* h = hists[p][lc];
      if (!h || h->GetEntries() <= 0) continue;

      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.9);
      h->GetYaxis()->SetTitleOffset(1.4);

      // Set a small positive minimum for log scale if needed
      if (h->GetMinimum() <= 0.0) h->SetMinimum(0.5);

      h->Draw();
      gPad->Modified();
      gPad->Update();
    }

    c->Modified();
    c->Update();
    gSystem->ProcessEvents();

    TString out = TString::Format("%s%s_RAYRAW#%02d.pdf",
                                  PATH_OUT_LY_EACHCH.c_str(),
                                  runTag.c_str(), p+1);
    c->SaveAs(out.Data());
  }
  // Restore previous AddDirectory setting
  TH1::AddDirectory(oldAddDir);
}

// ------------------------
// xg–yg 2D barycenter
// ------------------------

static inline bool CableToPosition(int cablenum, double& x, double& y) {
  x = std::numeric_limits<double>::quiet_NaN();
  y = std::numeric_limits<double>::quiet_NaN();
  if (1 <= cablenum && cablenum <= 140) { y = 695.0 - 10.0 * (cablenum - 1); return true; }
  if (201 <= cablenum && cablenum <= 332) { x = 655.0 - 10.0 * (cablenum - 201); return true; }
  return false;
}

static void AccumulateXY(TH2D* hXY, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  Int_t cablenum[NOUT];
  Double_t ly[NOUT][NBUNCH];

  t->SetBranchStatus("cablenum", 1);
  t->SetBranchStatus("lightyield", 1);
  t->SetBranchAddress("cablenum", cablenum);
  t->SetBranchAddress("lightyield", ly);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);

    for (int b = 0; b < NBUNCH; ++b) {
      double x_num = 0.0, x_den = 0.0;
      double y_num = 0.0, y_den = 0.0;
      double lightmax_x = 0.0, lightmax_y = 0.0;

      for (int i = 0; i < NOUT; ++i) {
        const double ly_b = ly[i][b];
        if (!(std::isfinite(ly_b) && ly_b > 0.0)) continue;

        double x_i, y_i;
        if (!CableToPosition(cablenum[i], x_i, y_i)) continue;

        if (std::isfinite(x_i)) {
          const double w = std::pow(ly_b, XG_WEIGHT);
          x_num += w * x_i;
          x_den += w;
          if (ly_b > lightmax_x) lightmax_x = ly_b;
        }
        if (std::isfinite(y_i)) {
          const double w = std::pow(ly_b, XG_WEIGHT);
          y_num += w * y_i;
          y_den += w;
          if (ly_b > lightmax_y) lightmax_y = ly_b;
        }
      }

      const bool ok_den = (x_den > 0.0) || (y_den > 0.0);
      if (ok_den && (lightmax_x > LIGHTMAX_MIN) && (lightmax_y > LIGHTMAX_MIN)) {
        const double xg = (x_den > 0.0) ? (x_num / x_den) : std::numeric_limits<double>::quiet_NaN();
        const double yg = (y_den > 0.0) ? (y_num / y_den) : std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(xg) && std::isfinite(yg)) hXY->Fill(xg, yg);
      }
    }
  }

  fin.Close();
}

static void BuildAndSaveXY2D(const std::string& runTag, const std::vector<FileInfoLite>& stableFiles) {
  gSystem->mkdir(PATH_OUT_XY.c_str(), true);
  std::unique_ptr<TH2D> hXY(new TH2D("hXY", ";x_{g} [mm];y_{g} [mm]",
                                     XBINS_XY, XMIN_XY, XMAX_XY,
                                     YBINS_XY, YMIN_XY, YMAX_XY));

  for (const auto& fi : stableFiles) AccumulateXY(hXY.get(), fi.path);

  gStyle->SetOptStat(0);
  TCanvas c("c_xy", "xg-yg", 1000, 900);
  c.SetGrid();
  c.SetTopMargin(0.15);
  hXY->SetContour(99);
  hXY->Draw("COLZ");
  hXY->GetYaxis()->SetTitleOffset(1.1);

  TDatime now;
  TPaveText note(0.7, 0.85, 0.95, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("Run: %s", runTag.c_str()));
  note.AddText(Form("Entries: %.0f", hXY->GetEntries()));
  note.AddText(Form("Mean x: %.2f, Std Dev x: %.2f", hXY->GetMean(1), hXY->GetStdDev(1)));
  note.AddText(Form("Mean y: %.2f, Std Dev y: %.2f", hXY->GetMean(2), hXY->GetStdDev(2)));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_xgyg.pdf", PATH_OUT_XY.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// LY profiles vs position (per run)
// ------------------------

static void BuildAndSaveLyVsPosition(
    const std::string& runTag,
    const std::vector<FileInfoLite>& stableFiles,
    const std::vector<double>& sum,   // per-channel sum of ly (>=10 p.e.)
    const std::vector<int>&    cnt)   // per-channel count of entries (>=10 p.e.)
{
  if (stableFiles.empty()) return;

  // (1) Map channel index -> cablenum from one ROOT file of this run
  std::vector<int> cab_of_ch(NOUT, -1);
  {
    TFile f(stableFiles.front().path.c_str(), "READ");
    if (f.IsZombie()) return;
    TTree* t = (TTree*)f.Get(TREE_NAME);
    if (!t) { f.Close(); return; }

    t->SetBranchStatus("*", 0);
    if (!t->GetBranch("cablenum")) { f.Close(); return; }

    Int_t cablenum_arr[NOUT];
    t->SetBranchStatus("cablenum", 1);
    t->SetBranchAddress("cablenum", cablenum_arr);

    if (t->GetEntries() <= 0) { f.Close(); return; }

    // Assume cablenum is fixed per channel and does not change with event
    t->GetEntry(0);
    for (int ch = 0; ch < NOUT; ++ch) {
      cab_of_ch[ch] = cablenum_arr[ch];
    }
    f.Close();
  }

  // (2) Define position bins
  const int    NX       = FrostmonConfig::XBINS_XY;
  const double XMIN_POS = FrostmonConfig::XMIN_XY;
  const double XMAX_POS = FrostmonConfig::XMAX_XY;
  const int    NY       = FrostmonConfig::YBINS_XY;
  const double YMIN_POS = FrostmonConfig::YMIN_XY;
  const double YMAX_POS = FrostmonConfig::YMAX_XY;

  // Average light yield vs position
  std::unique_ptr<TH1D> hXavg(new TH1D(
      "h_lyavg_xpos",
      "Average light yield vs fiber position (x);Fiber position (x) [mm];Average light yield [p.e.]",
      NX, XMIN_POS, XMAX_POS));

  std::unique_ptr<TH1D> hYavg(new TH1D(
      "h_lyavg_ypos",
      "Average light yield vs fiber position (y);Fiber position (y) [mm];Average light yield [p.e.]",
      NY, YMIN_POS, YMAX_POS));

  // Number of entries with ly >= 10 p.e. vs position
  // For example: if cablenum=1 has 80 entries with ly>=10 in this run,
  // the corresponding bin height will be 80.
  std::unique_ptr<TH1D> hXcount(new TH1D(
      "h_nentries_over10_xpos",
      "Number of events with light yield >= 10 p.e. vs fiber position (x);Fiber position (x) [mm];Number of events (light yeild >= 10 p.e.)",
      NX, XMIN_POS, XMAX_POS));

  std::unique_ptr<TH1D> hYcount(new TH1D(
      "h_nentries_over10_ypos",
      "Number of events with light yield >= 10 p.e. vs fiber position (y);Fiber position (y) [mm];Number of events (light yield >= 10 p.e.)",
      NY, YMIN_POS, YMAX_POS));

  hXavg->SetStats(0);
  hYavg->SetStats(0);
  hXcount->SetStats(0);
  hYcount->SetStats(0);
  hXavg->SetMinimum(0);
  hYavg->SetMinimum(0);
  hXcount->SetMinimum(0);
  hYcount->SetMinimum(0);

  // (3) Loop over channels and fill profiles
  for (int ch = 0; ch < NOUT; ++ch) {
    const int cab = cab_of_ch[ch];
    if (cab <= 0) continue;

    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    if (!CableToPosition(cab, x, y)) continue;

    const int    nEntries = cnt[ch];  // number of entries with ly>=10 for this channel in this run
    const double avg_ly   = (nEntries > 0)
                            ? (sum[ch] / static_cast<double>(nEntries))
                            : 0.0;

    // X side (cablenum 201–332)
    if (std::isfinite(x)) {
      const int binx = hXavg->FindBin(x);

      if (nEntries > 0 && avg_ly > 0.0) {
        // Set the average light yield for this position.
        // If multiple channels somehow map to the same x-bin,
        // we average them.
        const double prev_avg   = hXavg->GetBinContent(binx);
        const double prev_nchan = hXcount->GetBinContent(binx); // reuse as a "channel counter" for averaging
        if (prev_nchan <= 0.0) {
          hXavg->SetBinContent(binx, avg_ly);
        } else {
          const double new_avg = (prev_avg * prev_nchan + avg_ly) / (prev_nchan + 1.0);
          hXavg->SetBinContent(binx, new_avg);
        }
      }

      // For the count histogram, the bin height is the total number of entries (ly>=10)
      // for channels in this position.
      if (nEntries > 0) {
        hXcount->SetBinContent(binx, hXcount->GetBinContent(binx) + nEntries);
      }
    }

    // Y side (cablenum 1–140)
    if (std::isfinite(y)) {
      const int biny = hYavg->FindBin(y);

      if (nEntries > 0 && avg_ly > 0.0) {
        const double prev_avg   = hYavg->GetBinContent(biny);
        const double prev_nchan = hYcount->GetBinContent(biny); // reuse as a "channel counter" for averaging
        if (prev_nchan <= 0.0) {
          hYavg->SetBinContent(biny, avg_ly);
        } else {
          const double new_avg = (prev_avg * prev_nchan + avg_ly) / (prev_nchan + 1.0);
          hYavg->SetBinContent(biny, new_avg);
        }
      }

      if (nEntries > 0) {
        hYcount->SetBinContent(biny, hYcount->GetBinContent(biny) + nEntries);
      }
    }
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  gStyle->SetOptStat(0);

  // (4) Canvas 1: average lightyield vs position (top: X, bottom: Y)
  {
    TCanvas c1(Form("c_ly_pos_avg_%s", runTag.c_str()),
               "Average LY vs position", 1200, 900);
    c1.Divide(1, 2);

    TDatime now;

    c1.cd(1);
    gPad->SetGrid();
    hXavg->Draw("HIST");

    TPaveText note(0.8, 0.9, 0.99, 0.99, "NDC");
    note.SetFillColor(0);
    note.SetTextAlign(12);
    note.AddText(Form("Run: %s", runTag.c_str()));
    note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                      now.GetYear(), now.GetMonth(), now.GetDay(),
                      now.GetHour(), now.GetMinute(), now.GetSecond()));
    note.Draw();

    c1.cd(2);
    gPad->SetGrid();
    hYavg->Draw("HIST");

    TString outpdf = TString::Format(
        "%s%s_chavg_lightyield_profile_xy.pdf",
        PATH_OUT_LY.c_str(), runTag.c_str());
    c1.SaveAs(outpdf);
  }

  // (5) Canvas 2: number of entries (ly>=10) vs position (top: X, bottom: Y)
  {
    TCanvas c2(Form("c_ly_pos_nentries_%s", runTag.c_str()),
               "Entries with LY >= 10 p.e. vs position", 1200, 900);
    c2.Divide(1, 2);

    TDatime now;

    c2.cd(1);
    gPad->SetGrid();
    hXcount->Draw("HIST");

    TPaveText note(0.8, 0.9, 0.99, 0.99, "NDC");
    note.SetFillColor(0);
    note.SetTextAlign(12);
    note.AddText(Form("Run: %s", runTag.c_str()));
    note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                      now.GetYear(), now.GetMonth(), now.GetDay(),
                      now.GetHour(), now.GetMinute(), now.GetSecond()));
    note.Draw();

    c2.cd(2);
    gPad->SetGrid();
    hYcount->Draw("HIST");

    TString outpdf = TString::Format(
        "%s%s_nevents_over10_profile_xy.pdf",
        PATH_OUT_LY.c_str(), runTag.c_str());
    c2.SaveAs(outpdf);
  }
}

// ------------------------
// evnum vs unixtime
// ------------------------

static void AccumulateEvTime(std::vector<double>& v_ev, std::vector<double>& v_ut, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  Int_t   evnum = 0;
  Double_t unixtime_arr[NOUT];
  t->SetBranchStatus("evnum", 1);
  t->SetBranchStatus("unixtime", 1);
  t->SetBranchAddress("evnum", &evnum);
  t->SetBranchAddress("unixtime", unixtime_arr);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    const double ut = unixtime_arr[0];
    if (std::isfinite(ut)) { v_ev.push_back((double)evnum); v_ut.push_back(ut); }
  }
  fin.Close();
}

static void BuildAndSaveEvTime(const std::string& runTag, const std::vector<FileInfoLite>& stableFiles) {
  gSystem->mkdir(PATH_OUT_UNIXTIME.c_str(), true);

  std::vector<double> v_ev; v_ev.reserve(100000);
  std::vector<double> v_ut; v_ut.reserve(100000);

  for (const auto& fi : stableFiles) AccumulateEvTime(v_ev, v_ut, fi.path);
  if (v_ev.empty()) return;

  TGraph gr((int)v_ev.size(), v_ev.data(), v_ut.data());
  gr.SetTitle("evnum vs unixtime;evnum;unixtime [s]");
  gr.SetMarkerStyle(20);
  gr.SetMarkerSize(0.5);
  gr.GetYaxis()->SetTitleOffset(1.8);

  TCanvas c("c_evtime", "evnum vs unixtime", 1000, 700);
  c.SetGrid();
  c.SetLeftMargin(0.15);
  gr.Draw("AP");

  TDatime now;
  TPaveText note(0.7, 0.9, 0.95, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("Run: %s", runTag.c_str()));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_evnum_vs_unixtime.pdf", PATH_OUT_UNIXTIME.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// evnum vs spillnum
// ------------------------

static void AccumulateEvSpillnum(std::vector<double>& v_ev, std::vector<double>& v_sn, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  Int_t evnum = 0;
  Int_t spillnum = -1;
  t->SetBranchStatus("evnum", 1);
  t->SetBranchStatus("spillnum", 1);
  t->SetBranchAddress("evnum", &evnum);
  t->SetBranchAddress("spillnum", &spillnum);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    const int sn = spillnum;
    // spillnum is scalar; no need for isfinite()
    v_ev.push_back((double)evnum);
    v_sn.push_back((double)sn);
  }
  fin.Close();
}

static void BuildAndSaveSpillnum(const std::string& runTag, const std::vector<FileInfoLite>& stableFiles) {
  gSystem->mkdir(PATH_OUT_SPILLNUM.c_str(), true);

  std::vector<double> v_ev; v_ev.reserve(100000);
  std::vector<double> v_sn; v_sn.reserve(100000);

  for (const auto& fi : stableFiles) AccumulateEvSpillnum(v_ev, v_sn, fi.path);
  if (v_ev.empty()) return;

  TGraph gr((int)v_ev.size(), v_ev.data(), v_sn.data());
  gr.SetTitle("evnum vs spillnum;evnum;spill number");
  gr.SetMarkerStyle(20);
  gr.SetMarkerSize(0.5);
  gr.SetMinimum(0);
  gr.SetMaximum(35000);
  gr.GetYaxis()->SetTitleOffset(1.5);

  TCanvas c("c_spill", "evnum vs spillnum", 1000, 700);
  c.SetGrid();
  gr.Draw("AP");

  TDatime now;
  TPaveText note(0.7, 0.9, 0.95, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("Run: %s", runTag.c_str()));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_evnum_vs_spillnum.pdf", PATH_OUT_SPILLNUM.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// 6-hour binned LY history
// ------------------------

struct SumCnt { double sum=0.0; long long cnt=0; };

static inline std::string key_bin_ch(long long bin, int ch) {
  return std::to_string(bin) + "\t" + std::to_string(ch);
}

static void LoadChAvgBins(std::unordered_map<std::string, SumCnt>& bins,
                          long long &bin_min, long long &bin_max)
{
  bins.clear();
  bin_min = std::numeric_limits<long long>::max();
  bin_max = std::numeric_limits<long long>::min();

  std::ifstream fin(PATH_LY_BINS);
  if (!fin) return;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    long long bin; int ch; double sum; long long cnt;
    if (!(iss >> bin >> ch >> sum >> cnt)) continue;
    bins[key_bin_ch(bin,ch)] = {sum, cnt};
    if (bin < bin_min) bin_min = bin;
    if (bin > bin_max) bin_max = bin;
  }
  if (bin_min==std::numeric_limits<long long>::max()) { bin_min=0; bin_max=-1; }
}

static void SaveChAvgBins(const std::unordered_map<std::string, SumCnt>& bins)
{
  std::ofstream fout(PATH_LY_BINS, std::ios::trunc);
  fout << "#bin\tch\tsum\tcnt\n";
  for (const auto& kv : bins) {
    std::istringstream iss(kv.first);
    long long bin; int ch;
    iss >> bin >> ch;
    fout << bin << "\t" << ch << "\t"
         << std::setprecision(10) << kv.second.sum << "\t" << kv.second.cnt << "\n";
  }
}

static void LoadProcessedFiles(std::unordered_map<std::string, Long_t>& seen)
{
  seen.clear();
  std::ifstream fin(PATH_LY_PROC);
  if (!fin) return;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    std::string path; Long_t mt=0;
    if (!(iss >> std::quoted(path) >> mt)) continue;
    seen[path] = mt;
  }
}

static void SaveProcessedFiles(const std::unordered_map<std::string, Long_t>& seen)
{
  std::ofstream fout(PATH_LY_PROC, std::ios::trunc);
  fout << "#path\tmtime\n";
  for (const auto& kv : seen) {
    fout << std::quoted(kv.first) << "\t" << kv.second << "\n";
  }
}

static void AccumulateBinsFromFile(std::unordered_map<std::string, SumCnt>& bins,
                                   const std::string& path)
{
  TFile f(path.c_str(), "READ");
  if (f.IsZombie()) return;
  TTree* t = (TTree*)f.Get(TREE_NAME);
  if (!t) { f.Close(); return; }

  t->SetBranchStatus("*", 0);
  Double_t ly[NOUT][NBUNCH];
  Double_t unixtime_arr[NOUT];
  t->SetBranchStatus("lightyield", 1);
  t->SetBranchStatus("unixtime",   1);
  t->SetBranchAddress("lightyield", ly);
  t->SetBranchAddress("unixtime",   unixtime_arr);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie=0; ie<nent; ++ie) {
    t->GetEntry(ie);
    const double ut = unixtime_arr[0];
    if (!std::isfinite(ut)) continue;
    const long long bin = (long long)std::floor(ut / BINW_SEC);

    for (int ch=0; ch<NOUT; ++ch) {
      for (int b=0; b<NBUNCH; ++b) {
        const double v = ly[ch][b];
        if (std::isfinite(v) && v >= 10.0) {
          auto &sc = bins[key_bin_ch(bin, ch)];
          sc.sum += v;
          sc.cnt += 1;
        }
      }
    }
  }
  f.Close();
}

static void UpdateLyAvgBinsIncremental()
{
  auto all = ListLightyieldFiles(PATH_ROOT);

  std::unordered_map<std::string, SumCnt> bins;
  long long bin_min=0, bin_max=-1;
  LoadChAvgBins(bins, bin_min, bin_max);

  std::unordered_map<std::string, Long_t> seen;
  LoadProcessedFiles(seen);

  bool changed_bins=false, changed_seen=false;

  for (const auto& fi : all) {
    if (!IsReadyLightyieldFile(fi)) continue;
    auto it = seen.find(fi.path);
    if (it != seen.end() && it->second == fi.mtime) continue;

    AccumulateBinsFromFile(bins, fi.path);
    seen[fi.path] = fi.mtime;
    changed_bins = true;
    changed_seen = true;
  }

  if (changed_bins) SaveChAvgBins(bins);
  if (changed_seen) SaveProcessedFiles(seen);
}

static void BuildAndSaveLyAvgHistory2D_Binned()
{
  UpdateLyAvgBinsIncremental();

  std::unordered_map<std::string, SumCnt> bins;
  long long bin_min = 0, bin_max = -1;
  LoadChAvgBins(bins, bin_min, bin_max);
  if (bin_max < bin_min) return;

  const double tmin = static_cast<double>(bin_min) * BINW_SEC;
  const double tmax = (static_cast<double>(bin_max) + 1.0) * BINW_SEC;
  const int nbx = std::max(1, static_cast<int>(bin_max - bin_min + 1));

  std::unique_ptr<TH2D> h(new TH2D(
      "h_lyavg_hist2d",
      "Average light yield history;time;Average light yield per channel (>=10 p.e.) [p.e.]",
      nbx, tmin, tmax,
      NBINS_LYAVG, XMIN_LYAVG, XMAX_LYAVG));

  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeOffset(0, "local");
  h->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelOffset(0.02);
  h->GetXaxis()->SetTitleOffset(1.3);

  long long npoints = 0;
  for (long long bin = bin_min; bin <= bin_max; ++bin) {
    long long total_cnt_in_bin = 0;
    for (int ch = 0; ch < NOUT; ++ch) {
      const auto it = bins.find(key_bin_ch(bin, ch));
      if (it != bins.end()) {
        total_cnt_in_bin += it->second.cnt;
      }
    }

    // Skip bins where the total number of ly>=10 p.e. entries
    // in this 6-hour window is too small (i.e. beam-off periods).
    if (total_cnt_in_bin < MIN_COUNTS_PER_BIN) continue;

    const double t = static_cast<double>(bin) * BINW_SEC;
    for (int ch = 0; ch < NOUT; ++ch) {
      const auto it = bins.find(key_bin_ch(bin, ch));
      double avg = 0.0;
      if (it != bins.end() && it->second.cnt > 0) {
        avg = it->second.sum / static_cast<double>(it->second.cnt);
      }
      h->Fill(t, avg);
      ++npoints;
    }
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  gStyle->SetOptStat(0);

  TCanvas c("c_lyavg_hist2d", "LY avg history (6h bins)", 2000, 800);
  c.SetGrid();
  h->SetContour(99);
  h->GetXaxis()->SetNdivisions(520, kTRUE);

  // Limit visible X range to >= 2025-11-29 18:00 (local time).
  {
    const double t_cut  = unixTimeLocal(2025, 11, 29, 18, 0, 0);
    const double x_min_plot = std::max(tmin, t_cut);
    const double x_max_plot = tmax;
    if (x_min_plot < x_max_plot) {
      h->GetXaxis()->SetRangeUser(x_min_plot, x_max_plot);
    }
  }

  h->Draw("COLZ");
  h->GetYaxis()->SetTitleOffset(1.1);

  TDatime now;
  TPaveText note(0.75, 0.85, 0.99, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText("Run: ALL (6h bins)");
  note.AddText(Form("Entries: %lld", npoints));
  note.AddText(Form("Updated: %04d/%02d/%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%sALL_chavg_lightyield_history_2D_6h.pdf", PATH_OUT_LY.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// Helpers for per-run PDF management
// ------------------------

// Check if all expected PDFs for a given run already exist.
// If any of them is missing, this returns false.
static bool RunPdfsExist(const std::string& runTag)
{
  auto fileExists = [](const std::string& p) {
    // AccessPathName returns 0 if the file exists.
    return (gSystem->AccessPathName(p.c_str()) == kFALSE);
  };

  // Timing histograms (TDC / ADC)
  const std::string baseT = PATH_OUT_TDC + runTag;
  if (!fileExists(baseT + "_leading_hist.pdf"))          return false;
  if (!fileExists(baseT + "_trailing_hist.pdf"))         return false;
  if (!fileExists(baseT + "_leading_fromadc_hist.pdf"))  return false;
  if (!fileExists(baseT + "_trailing_fromadc_hist.pdf")) return false;

  // Average lightyield per channel
  if (!fileExists(PATH_OUT_LY + runTag + "_chavg_lightyield_hist.pdf")) return false;

  // Average lightyield vs position (X/Y)
  if (!fileExists(PATH_OUT_LY + runTag + "_chavg_lightyield_profile_xy.pdf")) return false;

  // Number of eventss (ly>=10 p.e.) vs position (X/Y)
  if (!fileExists(PATH_OUT_LY + runTag + "_nevents_over10_profile_xy.pdf")) return false;

  // xg–yg 2D
  if (!fileExists(PATH_OUT_XY + runTag + "_xgyg.pdf")) return false;

  // evnum vs unixtime
  if (!fileExists(PATH_OUT_UNIXTIME + runTag + "_evnum_vs_unixtime.pdf")) return false;

  // evnum vs spillnum
  if (!fileExists(PATH_OUT_SPILLNUM + runTag + "_evnum_vs_spillnum.pdf")) return false;

  // RAYRAW-wise per-channel lightyield PDFs
  for (int p = 1; p <= FrostmonConfig::N_RAYRAW; ++p) {
    TString name = TString::Format("%s%s_RAYRAW#%02d.pdf",
                                   PATH_OUT_LY_EACHCH.c_str(),
                                   runTag.c_str(), p);
    if (!fileExists(name.Data())) return false;
  }

  return true;
}

// Collect "ready" lightyield files for a given runTag from the full file list.
// Returns true if at least one ready file is found.
static bool CollectStableFilesForRun(const std::string& runTag,
                                     const std::vector<FileInfoLite>& allFiles,
                                     std::vector<FileInfoLite>& stableFilesOut)
{
  stableFilesOut.clear();
  auto runFiles = FilterByRun(allFiles, runTag);
  for (const auto& fi : runFiles) {
    if (IsReadyLightyieldFile(fi)) {
      stableFilesOut.push_back(fi);
    }
  }
  return !stableFilesOut.empty();
}

// Produce all per-run PDFs (timing, average LY, per-channel LY, xg-yg,
// evnum vs unixtime, evnum vs spillnum) for a given run and its ready files.
static void MakeRunPdfs(const std::string& runTag,
                        const std::vector<FileInfoLite>& stableFiles)
{
  // Timing histograms
  std::unique_ptr<TH1D> h_lead(
      CreateHistTiming("h_leading", true,  "Leading","leading"));
  std::unique_ptr<TH1D> h_trail(
      CreateHistTiming("h_trailing", true, "Trailing","trailing"));
  std::unique_ptr<TH1D> h_lead_adc(
      CreateHistTiming("h_leading_adc", false,  "Leading_fromADC","leading_fromADC"));
  std::unique_ptr<TH1D> h_trail_adc(
      CreateHistTiming("h_trailing_adc", false, "Trailing_fromADC","trailing_fromADC"));

  for (const auto& fi : stableFiles) {
    FillHistsTimingFromFile(h_lead.get(), h_trail.get(),
                            h_lead_adc.get(), h_trail_adc.get(),
                            fi.path);
  }

  DrawAndSaveTiming(h_lead.get(),      runTag, "leading");
  DrawAndSaveTiming(h_trail.get(),     runTag, "trailing");
  DrawAndSaveTiming(h_lead_adc.get(),  runTag, "leading_fromadc");
  DrawAndSaveTiming(h_trail_adc.get(), runTag, "trailing_fromadc");

  // Average lightyield per channel (>=10 p.e.)
  std::vector<double> sum(NOUT, 0.0);
  std::vector<int>    cnt(NOUT, 0);
  for (const auto& fi : stableFiles) {
    AccumulateLyAvgPerChannel(sum, cnt, fi.path);
  }
  BuildAndSaveLyAvgHist(sum, cnt, runTag);

  // LY profiles vs position (average LY and #events with ly>=10 p.e.)
  BuildAndSaveLyVsPosition(runTag, stableFiles, sum, cnt);

  // Per-channel lightyield histograms (RAYRAW planes)
  BuildAndSaveLyEachChannel(runTag, stableFiles);

  // xg–yg 2D barycenter
  BuildAndSaveXY2D(runTag, stableFiles);

  // evnum vs unixtime
  BuildAndSaveEvTime(runTag, stableFiles);

  // evnum vs spillnum
  BuildAndSaveSpillnum(runTag, stableFiles);
}

// ------------------------
// Main loop
// ------------------------

static void OneIteration(int refresh_ms) {
  // Ensure output directories exist.
  gSystem->mkdir(PATH_OUT_TDC.c_str(),       true);
  gSystem->mkdir(PATH_OUT_LY.c_str(),        true);
  gSystem->mkdir(PATH_OUT_XY.c_str(),        true);
  gSystem->mkdir(PATH_OUT_UNIXTIME.c_str(),  true);
  gSystem->mkdir(PATH_OUT_SPILLNUM.c_str(),  true);
  gSystem->mkdir(PATH_OUT_LY_EACHCH.c_str(), true);

  // List all *_lightyield.root files (sorted by mtime in descending order).
  auto all = ListLightyieldFiles(PATH_ROOT);
  if (all.empty()) {
    ::Info("dataqualityplot", "No *_lightyield.root yet in: %s", PATH_ROOT.c_str());
    gSystem->Sleep(refresh_ms);
    return;
  }

  // Determine the latest runTag from the newest file that matches "runXXXXX".
  std::string latestRunTag = "";
  for (const auto& fi : all) {
    latestRunTag = ExtractRunTag(fi.path);
    if (!latestRunTag.empty()) break;
  }
  if (latestRunTag.empty()) {
    gSystem->Sleep(refresh_ms);
    return;
  }

  // 1. Always update PDFs for the latest run (online monitoring).
  {
    std::vector<FileInfoLite> stableFilesLatest;
    if (CollectStableFilesForRun(latestRunTag, all, stableFilesLatest)) {
      MakeRunPdfs(latestRunTag, stableFilesLatest);
    } else {
      ::Info("dataqualityplot",
             "No ready *_lightyield.root for latest run %s "
             "(require size >= 10kB and stable >= 60s).",
             latestRunTag.c_str());
    }
  }

  // 2. Backfill older runs: if some run does not have all PDFs yet, create them once.
  {
    // Collect unique run tags (excluding the latest run), in descending mtime order.
    std::vector<std::string> runTags;
    for (const auto& fi : all) {
      std::string tag = ExtractRunTag(fi.path);
      if (tag.empty() || tag == latestRunTag) continue;
      if (std::find(runTags.begin(), runTags.end(), tag) == runTags.end()) {
        runTags.push_back(tag);
      }
    }

    for (const auto& tag : runTags) {
      // Skip runs that already have all PDFs.
      if (RunPdfsExist(tag)) continue;

      // Skip runs that have no ready ROOT files.
      std::vector<FileInfoLite> stableFiles;
      if (!CollectStableFilesForRun(tag, all, stableFiles)) continue;

      ::Info("dataqualityplot",
             "Backfilling plots for run %s (some PDFs are missing).",
             tag.c_str());
      MakeRunPdfs(tag, stableFiles);

      // To limit CPU load, only backfill one run per iteration.
      break;
    }
  }

  // 3. Lightyield history (6-hour bins, all runs) is updated incrementally as before.
  BuildAndSaveLyAvgHistory2D_Binned();

  gSystem->Sleep(refresh_ms);
}

int main(int argc, char** argv) {
  // Default: infinite loop; refresh 60s
  int max_loops = -1;
  int refresh_ms = 60000;

  // Simple CLI parsing
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if ((a == "-n" || a == "--loops") && i+1 < argc) {
      max_loops = std::stoi(argv[++i]);
    } else if ((a == "-r" || a == "--refresh") && i+1 < argc) {
      refresh_ms = std::stoi(argv[++i]);
    } else if (a == "-h" || a == "--help") {
      std::cout << "Usage: " << argv[0] << " [-n loops] [-r refresh_ms]\n"
                << "  -n, --loops     Number of iterations; negative for infinite (default -1)\n"
                << "  -r, --refresh   Refresh interval in milliseconds (default 60000)\n";
      return 0;
    }
  }

  int iter = 0;
  while (max_loops < 0 || iter < max_loops) {
    ++iter;
    OneIteration(refresh_ms);
  }

  ::Info("dataqualityplot", "Stopped.");
  return 0;
}
