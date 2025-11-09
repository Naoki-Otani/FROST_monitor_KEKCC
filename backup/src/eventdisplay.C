// eventdisplay.C
// 使い方：
// root -l
// .L eventdisplay.C+
// draw_event_rayraw_caen("rayraw.root", "caen.root", 0); // ev_start=0 から対話表示
//
// キー操作：Enter=次イベント、q=終了

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <string>
#include <limits>
#include <cmath>

static const int NCH_TOTAL = 272;
// ------- RAYRAW ヒスト設定（参照プログラムに準拠）-------
static const int    N_BINS_X   = 132;
static const double X_MIN_X    = -660.0;
static const double X_MAX_X    =  660.0;

static const int    N_BINS_Y   = 140;
static const double X_MIN_Y    = -700.0;
static const double X_MAX_Y    =  700.0;

// ------- CAEN ヒスト設定（-140..660 を 50 mm/bin で 16 ビン）-------
static const int    N_BINS_C   = 16;
static const double X_MIN_C    = -140.0;
static const double X_MAX_C    =  660.0;

// cab→ビン番号（1始まり）
// X: cab 332→bin1 (左=-660), 201→bin132 (右=+660)
static inline int cabToBinX(int cab) { return 333 - cab; }
// Y: cab 140→bin1 (左=-700), 1→bin140 (右=+700)
static inline int cabToBinY(int cab) { return 141 - cab; }

// CAEN: ch0→最右(+660)のビン, ch15→最左(-140)のビン
static inline int chToBinC(int ch) { return N_BINS_C - ch; } // ch=0→16, ch=15→1

void eventdisplay_(const char* rayraw_file,
                            const char* caen_file,
                            Long64_t ev_start = 0)
{
  // ----- 入力ファイルとツリー -----
  TFile fr(rayraw_file, "READ");
  if (fr.IsZombie()) { std::cerr << "Failed to open rayraw file: " << rayraw_file << "\n"; return; }
  TTree* tr = (TTree*)fr.Get("tree");
  if (!tr) { std::cerr << "TTree 'tree' not found in rayraw file.\n"; return; }

  TFile fc(caen_file, "READ");
  if (fc.IsZombie()) { std::cerr << "Failed to open caen file: " << caen_file << "\n"; return; }
  TTree* tc = (TTree*)fc.Get("t");
  if (!tc) { std::cerr << "TTree 't' not found in caen file.\n"; return; }

  // ----- ブランチ（必要最小限のみ）-----
  tr->SetBranchStatus("*", 0);
  tc->SetBranchStatus("*", 0);

  int    cablenum[NCH_TOTAL];
  // double adcint[NCH_TOTAL];
  double lightyield_RAYRAW[NCH_TOTAL];
  tr->SetBranchStatus("cablenum", 1);
  // tr->SetBranchStatus("adcint",   1);
  tr->SetBranchStatus("lightyield",   1);
  tr->SetBranchAddress("cablenum", cablenum);
  // tr->SetBranchAddress("adcint", adcint);
  tr->SetBranchAddress("lightyield", lightyield_RAYRAW);

  double lightyield_CAEN[16];
  tc->SetBranchStatus("lightyield", 1);
  tc->SetBranchAddress("lightyield", lightyield_CAEN);

  // ----- エントリ数（安全のため短い方に合わせる）-----
  const Long64_t nr = tr->GetEntries();
  const Long64_t nc = tc->GetEntries()-1;
  const Long64_t n  = std::min(nr, nc);
  if (ev_start < 0) ev_start = 0;
  if (ev_start >= n) {
    std::cerr << "ev_start >= nentries (" << ev_start << " >= " << n << ")\n";
    return;
  }
  if (nr != nc) {
    std::cerr << "[WARN] Entries differ (rayraw=" << nr << ", caen=" << nc
              << "). Using min=" << n << ".\n";
  }

  // ----- 画面 -----
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_rayraw_caen", "Event display: RAYRAW (X,Y) + CAEN light yield", 1400, 1000);
  c->Divide(1, 3); // 上: X, 中: Y, 下: CAEN

  std::string line;
  for (Long64_t i = ev_start; i < n; ++i) {
    tr->GetEntry(i);
    tc->GetEntry(i);

    // ヒスト作成（イベント毎にユニーク名 & メモリ管理）
    TH1D* hX = new TH1D(Form("hX_%lld", i),
                        Form("RAYRAW light yield X (Event %lld);Fiber position (x) [mm];Light yield [p.e.]", i),
                        N_BINS_X, X_MIN_X, X_MAX_X);
    TH1D* hY = new TH1D(Form("hY_%lld", i),
                        Form("RAYRAW light yield Y (Event %lld);Fiber position (y) [mm];Light yield [p.e.]", i),
                        N_BINS_Y, X_MIN_Y, X_MAX_Y);
    TH1D* hC = new TH1D(Form("hC_%lld", i),
                        Form("CAEN light yield (Event %lld);Fiber position (x) [mm];Light yield [p.e.]", i),
                        N_BINS_C, X_MIN_C, X_MAX_C);

    hX->SetDirectory(nullptr);
    hY->SetDirectory(nullptr);
    hC->SetDirectory(nullptr);

    // --- RAYRAW: cablenum と adcint から埋める ---
    for (int k = 0; k < NCH_TOTAL; ++k) {
      const int cab = cablenum[k];
      const double ly_RAYRAW = lightyield_RAYRAW[k];
      if (!std::isfinite(ly_RAYRAW)) continue;

      if (cab >= 201 && cab <= 332) {
        int bin = cabToBinX(cab);
        if (bin >= 1 && bin <= N_BINS_X) hX->SetBinContent(bin, ly_RAYRAW);
      } else if (cab >= 1 && cab <= 140) {
        int bin = cabToBinY(cab);
        if (bin >= 1 && bin <= N_BINS_Y) hY->SetBinContent(bin, ly_RAYRAW);
      }
      // それ以外の cab は無視
    }

    // --- CAEN: ch0→右端(+660mm)側、ch15→左端(-140mm)側 ---
    for (int ch = 0; ch < 16; ++ch) {
      const double ly_CAEN = lightyield_CAEN[ch];
      if (!std::isfinite(ly_CAEN)) continue;
      int bin = chToBinC(ch);
      if (bin >= 1 && bin <= N_BINS_C) hC->SetBinContent(bin, ly_CAEN);
    }

    // スタイル
    hX->SetLineColor(kBlue+1);  hX->SetLineWidth(2);
    hY->SetLineColor(kRed+1);   hY->SetLineWidth(2);
    hC->SetLineColor(kGreen+2); hC->SetLineWidth(2);

    // 描画
    c->cd(1); hX->Draw("HIST");
    c->cd(2); hY->Draw("HIST");
    c->cd(3); hC->Draw("HIST");

    c->Update();
    gSystem->ProcessEvents();

    std::cout << "Event " << i << " / " << (n-1)
              << "  —  <Enter>=next, 'q'=quit : " << std::flush;
    std::getline(std::cin, line);
    if (line == "q") break;
  }
}

// ROOT から関数形式で呼びやすいよう wrapper を用意
void eventdisplay(int runnum, Long64_t ev_start = 0)
{
  TString path_cosmictest = "/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/";
  TString filename_rayraw = TString::Format("%sdataaftercalib/rayraw/lightyield/run%05d_lightyield.root",path_cosmictest.Data(), runnum);
  // TString filename_rayraw = TString::Format("%sxgyg/xgyg_run%05d_weight4.0.root",path_cosmictest.Data(), runnum);
  TString filename_caen = TString::Format("%sdataaftercalib/caen/caenrun%02d_aftercalib.root",path_cosmictest.Data(), runnum);
  eventdisplay_(filename_rayraw, filename_caen, ev_start);
}
