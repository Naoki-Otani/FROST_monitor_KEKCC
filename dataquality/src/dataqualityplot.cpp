// dataqualityplot.cpp
// Purpose:
//   Standalone C++ (compiled) program using ROOT to generate data-quality plots.
//   It mirrors and extends the previous ROOT macro version.
//
// Major outputs (PDF):
//   - Timing histograms (leading / trailing / leading_fromadc / trailing_fromadc) for latest run
//   - Average lightyield per channel (>=10 p.e.) histogram for latest run
//   - xg–yg 2D barycenter (latest run)
//   - evnum vs unixtime (latest run)
//   - evnum vs spillnum (latest run)
//   - Average lightyield history 2D over all runs (6-hour bins) with incremental caching
//
// Build example:
//   g++ -O2 -std=c++17 dataqualityplot.cpp -o dataqualityplot $(root-config --cflags --libs)
//
// Run examples:
//   ./dataqualityplot                 // infinite loop, refresh 60s
//   ./dataqualityplot -n 10          // 10 iterations
//   ./dataqualityplot -n -1 -r 5000  // infinite loop, refresh every 5s

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

// ------------------------
// Global configuration
// ------------------------

static const char* TREE_NAME = "tree";

// Timing hist settings
static const int    NBINS_TDC = 2048;
static const double XMIN_TDC  = 0.0;
static const double XMAX_TDC  = 8192.0;

static const int    NBINS_ADC = 400;
static const double XMIN_ADC  = 0.0;
static const double XMAX_ADC  = 400.0;

// Lightyield settings
static const int    NOUT       = 272;
static const int    NBUNCH     = 8;

static const int    NBINS_LYAVG = 50;
static const double XMIN_LYAVG  = 0.0;
static const double XMAX_LYAVG  = 200.0;

// xg–yg settings
static const int    XBINS_XY  = 132;
static const int    YBINS_XY  = 140;
static const double XMIN_XY   = -660.0;
static const double XMAX_XY   =  660.0;
static const double YMIN_XY   = -700.0;
static const double YMAX_XY   =  700.0;

static const double XG_WEIGHT = 4.0;        // weight exponent
static const double LIGHTMAX_MIN = 10.0;    // threshold for xg–yg selection

// 6-hour binning for LY history
static const std::string PATH_LY_PROC  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/processed_files.tsv";
static const std::string PATH_LY_BINS  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/chavg_bins.tsv";
static const double BINW_SEC = 6.0 * 3600.0;

// I/O paths
static const std::string PATH_ROOT     = "/group/nu/ninja/work/otani/FROST_beamdata/test/rootfile_aftercalib/";
static const std::string PATH_OUT_TDC  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/";
static const std::string PATH_OUT_LY   = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/";
static const std::string PATH_OUT_XY   = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/xgyg/";
static const std::string PATH_OUT_UNIXTIME = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/unixtime/";
static const std::string PATH_OUT_SPILLNUM = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/spillnum/";

// ------------------------
// Utilities
// ------------------------

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
//   - size >= 1 MB
//   - file size stable (short check via IsLikelyStableFile)
//   - last modification time older than 60 seconds
static bool IsReadyLightyieldFile(const FileInfoLite& fi,
                                  Long_t minSizeBytes = 1L * 1024 * 1024,
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
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText(Form("entries: %.0f", h->GetEntries()));
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
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
  std::unique_ptr<TH1D> h(new TH1D("h_chavg_ly", "Average lightyield per channel (>=10 p.e.);Lightyield [p.e.];Channels",
                                   NBINS_LYAVG, XMIN_LYAVG, XMAX_LYAVG));
  h->SetLineWidth(2);
  h->SetStats(0);
  for (int ch = 0; ch < NOUT; ++ch) {
    const double avg = (cnt[ch] > 0) ? (sum[ch] / (double)cnt[ch]) : 0.0;
    h->Fill(avg);
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  TCanvas c("c_chavg_ly", "chavg ly", 1000, 700);
  c.SetGrid();
  h->SetMinimum(0.5);
  h->Draw("HIST");

  TDatime now;
  TPaveText note(0.75, 0.80, 0.99, 0.9, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText("entries: 272");
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_chavg_lightyield_hist.pdf", PATH_OUT_LY.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
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
  hXY->SetContour(99);
  hXY->Draw("COLZ");
  hXY->GetYaxis()->SetTitleOffset(1.1);

  TDatime now;
  TPaveText note(0.7, 0.9, 0.95, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText(Form("entries: %.0f", hXY->GetEntries()));
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%s%s_xgyg.pdf", PATH_OUT_XY.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
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
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
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
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
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
      "Average lightyield history;time;Average lightyield per channel (>=10 p.e.) [p.e.]",
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
    bool any_data = false;
    for (int ch = 0; ch < NOUT; ++ch) {
      const auto it = bins.find(key_bin_ch(bin, ch));
      if (it != bins.end() && it->second.cnt > 0) { any_data = true; break; }
    }
    if (!any_data) continue;

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
  h->Draw("COLZ");
  h->GetYaxis()->SetTitleOffset(1.1);

  TDatime now;
  TPaveText note(0.75, 0.85, 0.99, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText("run: ALL (6h bins)");
  note.AddText(Form("entries: %lld", npoints));
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%sALL_chavg_lightyield_history_2D_6h.pdf", PATH_OUT_LY.c_str());
  c.SaveAs(outpdf);
}

// ------------------------
// Main loop
// ------------------------

static void OneIteration(int refresh_ms) {
  // Make sure output directories exist
  gSystem->mkdir(PATH_OUT_TDC.c_str(),  true);
  gSystem->mkdir(PATH_OUT_LY.c_str(),   true);
  gSystem->mkdir(PATH_OUT_XY.c_str(),   true);
  gSystem->mkdir(PATH_OUT_UNIXTIME.c_str(), true);
  gSystem->mkdir(PATH_OUT_SPILLNUM.c_str(), true);

  // List all files and decide latest run
  auto all = ListLightyieldFiles(PATH_ROOT);
  if (all.empty()) {
    ::Info("dataqualityplot", "No *_lightyield.root yet in: %s", PATH_ROOT.c_str());
    gSystem->Sleep(refresh_ms);
    return;
  }

  std::string runTag = "";
  for (const auto& fi : all) { runTag = ExtractRunTag(fi.path); if (!runTag.empty()) break; }
  if (runTag.empty()) { gSystem->Sleep(refresh_ms); return; }

  // Collect "ready" files for the latest run
  auto runFiles = FilterByRun(all, runTag);
  std::vector<FileInfoLite> stableFiles;
  for (const auto& fi : runFiles) {
    if (IsReadyLightyieldFile(fi)) stableFiles.push_back(fi);
  }

  if (stableFiles.empty()) {
    ::Info("dataqualityplot",
           "No ready *_lightyield.root for run %s "
           "(require size >= 1MB and stable >= 60s).",
           runTag.c_str());
    gSystem->Sleep(refresh_ms);
    return;
  }

  // Timing histograms
  std::unique_ptr<TH1D> h_lead(CreateHistTiming("h_leading", true,  "Leading","leading"));
  std::unique_ptr<TH1D> h_trail(CreateHistTiming("h_trailing", true, "Trailing","trailing"));
  std::unique_ptr<TH1D> h_lead_adc(CreateHistTiming("h_leading_adc", false,  "Leading_fromADC","leading_fromADC"));
  std::unique_ptr<TH1D> h_trail_adc(CreateHistTiming("h_trailing_adc", false, "Trailing_fromADC","trailing_fromADC"));

  for (const auto& fi : stableFiles) {
    FillHistsTimingFromFile(h_lead.get(), h_trail.get(), h_lead_adc.get(), h_trail_adc.get(), fi.path);
  }

  DrawAndSaveTiming(h_lead.get(),      runTag, "leading");
  DrawAndSaveTiming(h_trail.get(),     runTag, "trailing");
  DrawAndSaveTiming(h_lead_adc.get(),  runTag, "leading_fromadc");
  DrawAndSaveTiming(h_trail_adc.get(), runTag, "trailing_fromadc");

  // Average LY per channel
  std::vector<double> sum(NOUT, 0.0);
  std::vector<int>    cnt(NOUT, 0);
  for (const auto& fi : stableFiles) AccumulateLyAvgPerChannel(sum, cnt, fi.path);
  BuildAndSaveLyAvgHist(sum, cnt, runTag);

  // xg–yg
  BuildAndSaveXY2D(runTag, stableFiles);

  // evnum–unixtime
  BuildAndSaveEvTime(runTag, stableFiles);

  // evnum–spillnum
  BuildAndSaveSpillnum(runTag, stableFiles);

  // LY avg history 2D (all runs, 6h bins, incremental)
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
