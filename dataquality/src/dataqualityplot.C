// dataqualityplot.C
// 目的:
//   既存: 最新ラン番号の *_lightyield.root 群から
//         leading / trailing / leading_fromadc / trailing_fromadc のヒストを更新（書き込み中は除外）
//   追加: ① ch毎の「10 p.e.以上のlightyield平均」ヒスト（272エントリ）
//        ② xg–yg の 2D ヒスト（バリセンタ、以前の実装に準拠）
//        ③ evnum–unixtime のグラフ
//
// 実行例:
//   root -l -q 'dataqualityplot.C()'
//   root -l -q 'dataqualityplot.C(-1, 60000)'  // 無限ループ, 60秒ごと更新
//
// 出力 (timing: 変更なし):
//   /group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/runXXXXX_<name>_hist.pdf
//   <name> = leading | trailing | leading_fromadc | trailing_fromadc
//
// 出力 (追加分):
//   ① ch平均LYヒスト: /group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/runXXXXX_chavg_lightyield_hist.pdf
//   ② xg–yg 2Dヒスト: /group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/xgyg/runXXXXX_xgyg.pdf
//   ③ evnum–unixtime:  /group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/unixtime/runXXXXX_evnum_vs_unixtime.pdf

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
#include <iomanip>

static const char* TREE_NAME = "tree";

// --- 既存 timing 用（変更なし）---
// TDC (leading, trailing)
static const int    NBINS_TDC = 2048;
static const double XMIN_TDC  = 0.0;
static const double XMAX_TDC  = 8192.0;
// ADC (leading_fromadc, trailing_fromadc)
static const int    NBINS_ADC = 400;
static const double XMIN_ADC  = 0.0;
static const double XMAX_ADC  = 400.0;

// --- 追加分の定数 ---
static const int    NOUT       = 272;
static const int    NBUNCH     = 8;
// ① ch平均 lightyield ヒスト（範囲は汎用的に 0〜400 p.e.）
static const int    NBINS_LYAVG = 50;
static const double XMIN_LYAVG  = 0.0;
static const double XMAX_LYAVG  = 200.0;
// ② xg–yg 2D ヒスト（以前の実装に準拠）
static const int    XBINS_XY  = 132;
static const int    YBINS_XY  = 140;
static const double XMIN_XY   = -660.0;
static const double XMAX_XY   =  660.0;
static const double YMIN_XY   = -700.0;
static const double YMAX_XY   =  700.0;
// xg 重み（以前の実装と同義。必要なら変更可）
static const double XG_WEIGHT = 4.0;
// xgyg ヒストに入れるイベントの最小 light 最大値
static const double LIGHTMAX_MIN = 10.0;

// ---- Lightyield average history (全ラン・6時間ビン) ----
static const std::string PATH_LY_PROC  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/processed_files.tsv";
static const std::string PATH_LY_BINS  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/chavg_bins.tsv";
static const double BINW_SEC = 6.0 * 3600.0;  // 6 hours

// --- 入出力パス ---
static const std::string PATH_ROOT     = "/group/nu/ninja/work/otani/FROST_beamdata/test/rootfile_aftercalib/";
static const std::string PATH_OUT_TDC  = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/";
static const std::string PATH_OUT_LY   = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/";
static const std::string PATH_OUT_XY   = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/xgyg/";
static const std::string PATH_OUT_TIME = "/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/unixtime/";

// --- ユーティリティ（既存＋追加）---

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
  // 例: .../run00103_0_9999_lightyield.root -> run00103
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
    return a.mtime > b.mtime; // 新しい順
  });
  return out;
}

static std::vector<FileInfoLite> FilterByRun(const std::vector<FileInfoLite>& v, const std::string& runTag) {
  std::vector<FileInfoLite> out;
  const std::string needle = "/" + runTag + "_";
  for (const auto& fi : v) if (fi.path.find(needle) != std::string::npos) out.push_back(fi);
  return out;
}

// --- 既存 timing ヒスト作成（変更なしのまま再掲） ---
static TH1D* CreateHistTiming(const char* hname, bool isTDC, const char* titlePrefix) {
  if (isTDC) {
    TH1D* h = new TH1D(hname, Form("%s;Time [0.833 ns];Number of events", titlePrefix), NBINS_TDC, XMIN_TDC, XMAX_TDC);
    h->SetLineWidth(2);
    return h;
  } else {
    TH1D* h = new TH1D(hname, Form("%s;Time [13.3 ns];Number of events", titlePrefix), NBINS_ADC, XMIN_ADC, XMAX_ADC);
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

  // ★ 既存の note レイアウトは変更しない
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

// --- 追加① ch毎の10 p.e.以上平均 lightyield ---
// 272 ch ぶんの平均値を算出（全イベント・全バンチのうち >=10 p.e. を平均。該当なしは0）
static void AccumulateLyAvgPerChannel(std::vector<double>& sum, std::vector<int>& cnt, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  // lightyield[272][8]
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
    h->Fill(avg); // 272 回だけ Fill
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  TCanvas c("c_chavg_ly", "chavg ly", 1000, 700);
  c.SetGrid();
  // c.SetLogy();          // 分布を見やすくする（必要なければ外してOK）
  h->SetMinimum(0.5);   // 0ビン回避
  h->Draw("HIST");

  // 小さめの注記
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

// --- 追加② xg–yg 2D バリセンタ ---
// cablenum -> (x,y) 変換（以前の実装準拠）
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

    // ★ 8バンチそれぞれで重心を計算
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

      // ★ しきい値: lightmax_x>10 かつ lightmax_y>10
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

  // 2Dヒストを積算
  for (const auto& fi : stableFiles) AccumulateXY(hXY.get(), fi.path);

  // 描画
  gStyle->SetOptStat(0);
  TCanvas c("c_xy", "xg-yg", 1000, 900);
  c.SetGrid();
  hXY->SetContour(99);
  hXY->Draw("COLZ");
  hXY->GetYaxis()->SetTitleOffset(1.1);

  // ★ note を他と同様に追加（位置・サイズも合わせる）
  TDatime now;
  TPaveText note(0.7, 0.9, 0.95, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText(Form("run: %s", runTag.c_str()));
  note.AddText(Form("entries: %.0f", hXY->GetEntries()));  // Fill回数（= xg,yg 点数）
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  // 保存
  TString outpdf = TString::Format("%s%s_xgyg.pdf", PATH_OUT_XY.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// --- 追加③ evnum–unixtime グラフ ---
static void AccumulateEvTime(std::vector<double>& v_ev, std::vector<double>& v_ut, const std::string& path) {
  TFile fin(path.c_str(), "READ");
  if (fin.IsZombie()) return;
  TTree* t = (TTree*)fin.Get(TREE_NAME);
  if (!t) { fin.Close(); return; }

  t->SetBranchStatus("*", 0);
  Int_t   evnum = 0;
  Double_t unixtime_arr[NOUT]; // 出力側は unixtime[272] で同値が複製されている想定
  t->SetBranchStatus("evnum", 1);
  t->SetBranchStatus("unixtime", 1);
  t->SetBranchAddress("evnum", &evnum);
  t->SetBranchAddress("unixtime", unixtime_arr);

  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    const double ut = unixtime_arr[0]; // 272個同値なので先頭を使用
    if (std::isfinite(ut)) { v_ev.push_back((double)evnum); v_ut.push_back(ut); }
  }
  fin.Close();
}

static void BuildAndSaveEvTime(const std::string& runTag, const std::vector<FileInfoLite>& stableFiles) {
  gSystem->mkdir(PATH_OUT_TIME.c_str(), true);

  std::vector<double> v_ev; v_ev.reserve(100000);
  std::vector<double> v_ut; v_ut.reserve(100000);

  for (const auto& fi : stableFiles) AccumulateEvTime(v_ev, v_ut, fi.path);
  if (v_ev.empty()) return;

  TGraph gr((int)v_ev.size(), v_ev.data(), v_ut.data());
  gr.SetTitle(";evnum;unixtime [s]");
  gr.SetMarkerStyle(20);
  gr.SetMarkerSize(0.5);

  TCanvas c("c_evtime", "evnum vs unixtime", 1000, 700);
  c.SetGrid();
  gr.Draw("AP");

  TString outpdf = TString::Format("%s%s_evnum_vs_unixtime.pdf", PATH_OUT_TIME.c_str(), runTag.c_str());
  c.SaveAs(outpdf);
}

// ---- 6hビン × ch の累積キャッシュ： (bin,ch) -> (sum, cnt) ----
struct SumCnt { double sum=0.0; long long cnt=0; };

// フラットキー（bin と ch を 1つの文字列に）
static inline std::string key_bin_ch(long long bin, int ch) {
  return std::to_string(bin) + "\t" + std::to_string(ch);
}

// 読み込み： chavg_bins.tsv
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

// 書き出し： chavg_bins.tsv
static void SaveChAvgBins(const std::unordered_map<std::string, SumCnt>& bins)
{
  std::ofstream fout(PATH_LY_BINS, std::ios::trunc);
  fout << "#bin\tch\tsum\tcnt\n";
  for (const auto& kv : bins) {
    // kv.first = "bin\tch"
    std::istringstream iss(kv.first);
    long long bin; int ch;
    iss >> bin >> ch;
    fout << bin << "\t" << ch << "\t"
         << std::setprecision(10) << kv.second.sum << "\t" << kv.second.cnt << "\n";
  }
}

// 読み込み： processed_files.tsv  (path \t mtime)
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
    if (!(iss >> std::quoted(path) >> mt)) continue; // パスはquoted対応
    seen[path] = mt;
  }
}

// 書き出し： processed_files.tsv
static void SaveProcessedFiles(const std::unordered_map<std::string, Long_t>& seen)
{
  std::ofstream fout(PATH_LY_PROC, std::ios::trunc);
  fout << "#path\tmtime\n";
  for (const auto& kv : seen) {
    fout << std::quoted(kv.first) << "\t" << kv.second << "\n";
  }
}

// 1ファイルを読み、6時間ビン×ch に >=10 p.e. の値を加算
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
    // イベントの代表時間（全ch同値の想定）として unixtime[0] を使用
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

// 増分更新：新規（未処理）かつ安定なファイルだけ bins に反映し、キャッシュ保存
static void UpdateLyAvgBinsIncremental()
{
  // ディレクトリ内の全 *_lightyield.root を列挙
  auto all = ListLightyieldFiles(PATH_ROOT);

  // 既存キャッシュ読み込み
  std::unordered_map<std::string, SumCnt> bins;
  long long bin_min=0, bin_max=-1;
  LoadChAvgBins(bins, bin_min, bin_max);

  std::unordered_map<std::string, Long_t> seen;
  LoadProcessedFiles(seen);

  bool changed_bins=false, changed_seen=false;

  for (const auto& fi : all) {
    // 書き込み中の可能性が少しでもあるならスキップ
    if (!IsLikelyStableFile(fi.path, 500)) continue;
    auto it = seen.find(fi.path);
    if (it != seen.end() && it->second == fi.mtime) continue; // 既に取り込み済み(同 mtime)

    // 新規または mtime 変化 → 取り込み（※同名再生成に対する差し戻しは実装簡略化のため未対応）
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
  // まず増分更新（未取り込みファイルがあれば bins に反映）
  UpdateLyAvgBinsIncremental();

  // bins を読み出し
  std::unordered_map<std::string, SumCnt> bins;
  long long bin_min=0, bin_max=-1;
  LoadChAvgBins(bins, bin_min, bin_max);
  if (bin_max < bin_min) return;

  // X軸範囲（6hビン開始時刻の最小・最大）を bins から決定
  const double tmin = (double)bin_min * BINW_SEC;
  const double tmax = ((double)bin_max + 1.0) * BINW_SEC; // 末尾ビンの終端
  const int nbx = std::max(1, (int)(bin_max - bin_min + 1));

  // Y軸は既存の ch-avg LY ヒストと同様（NBINS_LYAVG, XMIN_LYAVG, XMAX_LYAVG）
  std::unique_ptr<TH2D> h(new TH2D("h_lyavg_hist2d",
                                   "Average lightyield history;time;Average lightyield per channel (>=10 p.e.) [p.e.]",
                                   nbx, tmin, tmax,
                                   NBINS_LYAVG, XMIN_LYAVG, XMAX_LYAVG));

  h->GetXaxis()->SetTimeDisplay(1);
  // Unix epoch(1970-01-01)基準にする。ローカル時刻表示なら "gmt" は外してOK
  h->GetXaxis()->SetTimeOffset(0, "local");  // or: h->GetXaxis()->SetTimeOffset(TDatime(1970,1,1,0,0,0).Convert(), "gmt");
  h->GetXaxis()->SetTimeFormat("#splitline{%Y-%m-%d}{%H:%M}");
  // h->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
  // h->GetXaxis()->SetNdivisions(508, kTRUE); // 目盛りを適度に間引く（必要に応じて調整）
  h->GetXaxis()->SetLabelSize(0.03);        // ラベル小さめ
  h->GetXaxis()->SetLabelOffset(0.02);
  h->GetXaxis()->SetTitleOffset(1.3);

  // 各 (bin, ch) の平均を詰める
  long long npoints = 0;
  for (long long bin = bin_min; bin <= bin_max; ++bin) {
    const double t = (double)bin * BINW_SEC;
    for (int ch=0; ch<NOUT; ++ch) {
      const auto it = bins.find(key_bin_ch(bin, ch));
      double avg = 0.0;
      if (it != bins.end() && it->second.cnt > 0) avg = it->second.sum / (double)it->second.cnt;
      h->Fill(t, avg);
      ++npoints;
    }
  }

  gSystem->mkdir(PATH_OUT_LY.c_str(), true);
  gStyle->SetOptStat(0);

  TCanvas c("c_lyavg_hist2d", "LY avg history (6h bins)", 1200, 800);
  c.SetGrid();
  h->SetContour(99);
  h->Draw("COLZ");
  h->GetYaxis()->SetTitleOffset(1.1);

  // note（右上、既存の見た目に近づける）
  TDatime now;
  TPaveText note(0.75, 0.85, 0.99, 0.99, "NDC");
  note.SetFillColor(0);
  note.SetTextAlign(12);
  note.AddText("run: ALL (6h bins)");
  note.AddText(Form("entries: %lld", npoints)); // おおよそ = (#bins) * 272
  note.AddText(Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                    now.GetYear(), now.GetMonth(), now.GetDay(),
                    now.GetHour(), now.GetMinute(), now.GetSecond()));
  note.Draw();

  TString outpdf = TString::Format("%sALL_chavg_lightyield_history_2D_6h.pdf", PATH_OUT_LY.c_str());
  c.SaveAs(outpdf);
}

// --- メイン監視ループ（既存の timing に追加出力を付与） ---
void dataqualityplot(int max_loops = -1, int refresh_ms = 60000)
{
  // 出力ディレクトリ作成
  gSystem->mkdir(PATH_OUT_TDC.c_str(),  true);
  gSystem->mkdir(PATH_OUT_LY.c_str(),   true);
  gSystem->mkdir(PATH_OUT_XY.c_str(),   true);
  gSystem->mkdir(PATH_OUT_TIME.c_str(), true);

  int iter = 0;
  while (max_loops < 0 || iter < max_loops) {
    ++iter;

    // 最新ファイル群の取得
    auto all = ListLightyieldFiles(PATH_ROOT);
    if (all.empty()) {
      ::Info("make_timing_hists_realtime_plus", "No *_lightyield.root yet in: %s", PATH_ROOT.c_str());
      gSystem->Sleep(refresh_ms);
      continue;
    }

    // 最新 run タグ
    std::string runTag = "";
    for (const auto& fi : all) { runTag = ExtractRunTag(fi.path); if (!runTag.empty()) break; }
    if (runTag.empty()) { gSystem->Sleep(refresh_ms); continue; }

    // 安定ファイルのみ
    auto runFiles = FilterByRun(all, runTag);
    std::vector<FileInfoLite> stableFiles;
    for (const auto& fi : runFiles) if (IsLikelyStableFile(fi.path, 500)) stableFiles.push_back(fi);

    // --- 既存4ヒスト ---
    std::unique_ptr<TH1D> h_lead(CreateHistTiming("h_leading", true,  "Leading"));
    std::unique_ptr<TH1D> h_trail(CreateHistTiming("h_trailing", true, "Trailing"));
    std::unique_ptr<TH1D> h_lead_adc(CreateHistTiming("h_leading_adc", false,  "Leading_fromADC"));
    std::unique_ptr<TH1D> h_trail_adc(CreateHistTiming("h_trailing_adc", false, "Trailing_fromADC"));

    for (const auto& fi : stableFiles) {
      FillHistsTimingFromFile(h_lead.get(), h_trail.get(), h_lead_adc.get(), h_trail_adc.get(), fi.path);
    }

    DrawAndSaveTiming(h_lead.get(),      runTag, "leading");
    DrawAndSaveTiming(h_trail.get(),     runTag, "trailing");
    DrawAndSaveTiming(h_lead_adc.get(),  runTag, "leading_fromadc");
    DrawAndSaveTiming(h_trail_adc.get(), runTag, "trailing_fromadc");

    // --- 追加① ch平均 lightyield ---
    std::vector<double> sum(NOUT, 0.0);
    std::vector<int>    cnt(NOUT, 0);
    for (const auto& fi : stableFiles) AccumulateLyAvgPerChannel(sum, cnt, fi.path);
    BuildAndSaveLyAvgHist(sum, cnt, runTag);

    // --- 追加② xg–yg 2D ---
    BuildAndSaveXY2D(runTag, stableFiles);

    // --- 追加③ evnum–unixtime ---
    BuildAndSaveEvTime(runTag, stableFiles);

    // --- 追加④ lightyield average の 2Dヒスト（全ラン・6時間ビン、増分更新） ---
    BuildAndSaveLyAvgHistory2D_Binned();


    gSystem->Sleep(refresh_ms);
  }

  ::Info("make_timing_hists_realtime_plus", "Stopped.");
}
