// LightyieldHistRayraw.C
// 使い方：
//   root -l -b -q 'LightyieldHistRayraw.C+("run*****")'
//         eg) oot -l -b -q 'LightyieldHistRayraw.C+("run00044")'
// 出力： outprefix_RAYRAW###.pdf（各面 32pad, y軸ログ）

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

static const int NPLANES = 11;   // RAYRAW = 1..11
static const int NCHPL   = 32;   // ch per plane
static const int NOUT    = 272;  // 配列長

// NaN/Inf 安全チェック
static inline bool isfinite_d(double x){ return std::isfinite(x); }

void LightyieldHistRayraw_(const char* infile,
                               const char* outprefix,
                               int nbins = 100,
                               double xmin = -10,
                               double xmax = 200.0)
{
  // gROOT->SetBatch(kFALSE);
  TH1::AddDirectory(kFALSE);
  // 入力
  TFile f(infile, "READ");
  if (f.IsZombie()) { ::Error("LightyieldHistRayraw","Open failed: %s", infile); return; }
  TTree* t = (TTree*)f.Get("tree");
  if (!t) { ::Error("LightyieldHistRayraw","TTree 'tree' not found"); return; }

  // ブランチ配列
  int    cablenum[NOUT];
  int    rayraw  [NOUT];
  int    rayrawch[NOUT];
  double lightyield[NOUT];
  double unixtime[NOUT]; // 未使用だが読み込み可能に保持

  t->SetBranchAddress("cablenum",   cablenum);
  t->SetBranchAddress("rayraw",     rayraw);
  t->SetBranchAddress("rayrawch",   rayrawch);
  t->SetBranchAddress("lightyield", lightyield);
  if (t->GetBranch("unixtime")) t->SetBranchAddress("unixtime", unixtime);

  // 1) (plane, ch) -> cable の対応表を用意
  int cable_map[NPLANES][NCHPL];
  for (int p=0; p<NPLANES; ++p) for (int lc=0; lc<NCHPL; ++lc) cable_map[p][lc] = -1;

  // 2) 最初のエントリから対応を埋める（0/1始まり両対応の正規化を例示）
  t->GetEntry(0);
  for (int i=0; i<NOUT; ++i){
    int rr_raw = rayraw[i];
    int lc_raw = rayrawch[i];
    if (rr_raw < 0 || lc_raw < 0) continue;
    cable_map[rr_raw-1][lc_raw] = cablenum[i];  // ← 対応を記録
  }

  // ヒスト作成（面×ch）
  std::vector<std::vector<TH1D*>> h(NPLANES, std::vector<TH1D*>(NCHPL, (TH1D*)nullptr));
  for (int p = 0; p < NPLANES; ++p) {
    for (int lc = 0; lc < NCHPL; ++lc) {
      TString hname  = Form("h_p%02d_ch%02d", p+1, lc);
      TString htitle = Form("RAYRAW#%d ch%02d (cable %03d);Lightyield [p.e.];Number of events", p+1, lc,cable_map[p][lc]);
      h[p][lc] = new TH1D(hname, htitle, nbins, xmin, xmax);
      h[p][lc]->SetLineWidth(1);
    }
  }

  // ループして fill
  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);

    // 各イベント：配列 272 要素を面・ch に振り分ける
    for (int i = 0; i < NOUT; ++i) {
      const int rr = rayraw[i];         // 1..11 を想定
      const int lc = rayrawch[i];       // 0..31 を想定
      if (rr < 1 || rr > NPLANES) continue;
      if (lc < 0 || lc >= NCHPL) continue;

      const double ly = lightyield[i];
      if (!isfinite_d(ly)) continue;    // NaN/Inf スキップ

      h[rr-1][lc]->Fill(ly);
      // cout<<ly<<endl;
    }
  }

  // 描画・保存（面ごとに8×4のキャンバス、y軸ログ）
  gStyle->SetOptStat(1110);
  gStyle->SetStatW(0.22);
  gStyle->SetStatH(0.16);

  std::vector<TCanvas*> canv(NPLANES, (TCanvas*)nullptr);
  for (int p = 0; p < NPLANES; ++p) {
    canv[p] = new TCanvas(Form("c%d", p+1), Form("RAYRAW#%d lightyield", p+1), 2000, 1400);
    canv[p]->Divide(8, 4);
    for (int lc = 0; lc < NCHPL; ++lc) {
      canv[p]->cd(lc+1);
      gPad->SetLogy(1); // y軸ログ
      TH1D* hh = h[p][lc];
      // printf("p=%d lc=%d: Entries=%lld, Integral(in-range)=%g, MaxBin=%g, Under=%g, Over=%g\n",
      //  p+1, lc,
      //  (Long64_t)hh->GetEntries(),
      //  hh->Integral(1, hh->GetNbinsX()),    // 範囲内ビンの合計
      //  hh->GetMaximum(),
      //  hh->GetBinContent(0),
      //  hh->GetBinContent(hh->GetNbinsX()+1));

      if (!hh || hh->GetEntries() <= 0) continue;

      hh->GetXaxis()->SetLabelSize(0.05);
      hh->GetYaxis()->SetLabelSize(0.05);
      hh->GetXaxis()->SetTitleSize(0.05);
      hh->GetYaxis()->SetTitleSize(0.05);
      hh->GetXaxis()->SetTitleOffset(0.9);
      hh->GetYaxis()->SetTitleOffset(0.9);
      hh->Draw();
      gPad->Modified();           // ← 変更があったことを通知
      gPad->Update();             // ← 実際に再描画
      // gSystem->ProcessEvents();   // ← GUI イベントを処理（macOS 等で効く）
    }
    canv[p]->Modified();
    canv[p]->Update();
    gSystem->ProcessEvents();
    TString out = TString::Format("%s_RAYRAW#%02d.pdf", outprefix, p+1);
    canv[p]->SaveAs(out.Data());
  }


  f.Close();
}

// 呼び出し用ラッパ
void LightyieldHistRayraw(const char* infile,
                               int nbins = 100,
                               double xmin = -10,
                               double xmax = 200.0)
{
  std::string path_cosmictest = "/group/nu/ninja/work/otani/FROST_cosmictest_dataandoutput/";
  TString filename_inroot = TString(path_cosmictest) + "dataaftercalib/rayraw/" + infile + "_lightyield.root";
  TString filename_out = TString(path_cosmictest) + "LightyieldHisto/lightyield_" + infile;

  LightyieldHistRayraw_(filename_inroot.Data(), filename_out.Data(), nbins, xmin, xmax);
}
