// calibrayraw.C（0pe固定窓＋1pe右側ピークをガウスフィット／CSV出力・プロット保存）
// 使い方:
//   root -l -b -q 'calibrayraw.C+("run00042","chmap_run40")'
//
// ※ バッチ描画にしたい場合は -b を付けるか、関数内で gROOT->SetBatch(kTRUE) を使用

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TAxis.h>

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>

void calibrayraw_(TString filename, TString chmapPath, TString outputcsv, TString outputplot) {

  // バッチ描画（ウィンドウ非表示）
  gROOT->SetBatch(kTRUE);

  // ---------- 入力 ----------
  const Int_t N_RAYRAW = 11;          // RAYRAW 面数（1..11）
  const Int_t N_CH_PER_PLANE = 32;    // 各面の ch 数（0..31）
  const Int_t HIST_NBIN = 200;
  const Double_t HIST_XMIN = -100.;
  const Double_t HIST_XMAX =  300.;

  // ベースライン推定パラメータ
  const Int_t BL_WIN = 15;            // 平均窓幅
  const Int_t BL_START_MAX = 70;      // 開始点の上限（0,5,...,70）
  const Int_t BL_STEP = 5;            // スライド刻み

  // 積分窓（可変窓：ピーク中心 ±(20,25)）
  const Int_t INT_LEFT  = 20;
  const Int_t INT_RIGHT = 25;

  // 0pe用 固定窓（1 〜 INT_LEFT+INT_RIGHT）
  const Int_t FIX_START = 1;
  const Int_t FIX_RANGE = INT_LEFT + INT_RIGHT;

  // ---- TSpectrum（ピークサーチ：1pe側）設定 ----
  const Int_t   PEAKS_MAX     = 2;     // 常に2つまで
  const Double_t SPEC_SIGMA   = 3.0;
  const Double_t SPEC_THRESH  = 0.05;

  // ---------- chmap 読み込み ----------
  std::map<long long, int> cabMap; // key = ((long long)rr << 32) | (unsigned)lc
  {
    std::ifstream fin(chmapPath.Data());
    if (fin) {
      std::string line;
      bool first = true;
      while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') { first = false; continue; }
        if (first) { first = false; continue; } // 先頭行ヘッダ無視
        std::istringstream iss(line);
        int cab, rr, lc;
        if (!(iss >> cab >> rr >> lc)) continue;
        long long key = ( (long long)rr << 32 ) | (unsigned)lc;
        if (cabMap.find(key) == cabMap.end()) cabMap[key] = cab;
      }
    } else {
      std::cerr << "[WARN] chmap を開けませんでした: " << chmapPath << "\n"
                << "       cablenum なし/未定義chはヒストを作らず空欄表示します。\n";
    }
  }
  auto cabOf = [&](int rr, int lc)->int {
    long long key = ( (long long)rr << 32 ) | (unsigned)lc;
    auto it = cabMap.find(key);
    if (it == cabMap.end()) return -1;
    return it->second;
  };

  // ---------- 入力 ROOT ----------
  TFile *f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "File open failed: " << filename << std::endl;
    return;
  }
  TTree *t = (TTree*) f->Get("tree");
  if (!t) {
    std::cerr << "TTree 'tree' not found in file." << std::endl;
    return;
  }

  // ---------- ブランチ ----------
  std::vector<std::vector<double>> *waveform = nullptr;  // [ch][sample]
  std::vector<int> *plane = nullptr;                      // [ch] -> 0..10 (RAYRAW#-1)
  t->SetBranchAddress("waveform", &waveform);
  if (t->GetBranch("plane")) t->SetBranchAddress("plane", &plane);
  t->GetEntry(0);

  // ---------- ヒスト（面×ch） ----------
  // hInt1: 可変窓（ピーク中心 ±(20,25)）→ 1pe フィット
  // hInt0: 固定窓（0 〜 INT_LEFT+INT_RIGHT）→ 0pe フィット
  std::vector<std::vector<TH1D*>> hInt1(N_RAYRAW, std::vector<TH1D*>(N_CH_PER_PLANE, (TH1D*)nullptr));
  std::vector<std::vector<TH1D*>> hInt0(N_RAYRAW, std::vector<TH1D*>(N_CH_PER_PLANE, (TH1D*)nullptr));

  for (Int_t p = 0; p < N_RAYRAW; ++p) {
    for (Int_t lc = 0; lc < N_CH_PER_PLANE; ++lc) {
      int rr  = p + 1;
      int cab = cabOf(rr, lc);
      if (cab < 0) { hInt1[p][lc] = nullptr; hInt0[p][lc] = nullptr; continue; }

      TString hname1  = Form("h1_p%02d_ch%02d", rr, lc);
      TString htitle1 = Form("RAYRAW#%d ch%02d (cablenum%d);ADC Integral;Events", rr, lc, cab);
      hInt1[p][lc] = new TH1D(hname1, htitle1, HIST_NBIN, HIST_XMIN, HIST_XMAX);
      hInt1[p][lc]->SetLineWidth(1);

      TString hname0  = Form("h0_p%02d_ch%02d", rr, lc);
      TString htitle0 = Form("RAYRAW#%d ch%02d FIX[0..%d] (cablenum%d);ADC Integral;Events", rr, lc, FIX_RANGE, cab);
      hInt0[p][lc] = new TH1D(hname0, htitle0, HIST_NBIN, HIST_XMIN, HIST_XMAX);
      hInt0[p][lc]->SetLineWidth(1);
    }
  }

  // ---------- 積分ループ ----------
  const Long64_t nEntry = t->GetEntries();
  for (Long64_t iEntry = 0; iEntry < nEntry; ++iEntry) {
    t->GetEntry(iEntry);
    if (!waveform) continue;

    const Int_t nChAll = (Int_t) waveform->size();
    for (Int_t ch = 0; ch < nChAll; ++ch) {
      Int_t plane_id = (plane && (Int_t)plane->size() > ch) ? plane->at(ch) : (ch / N_CH_PER_PLANE);
      if (plane_id < 0 || plane_id >= N_RAYRAW) continue;

      const Int_t local_ch = ch % N_CH_PER_PLANE;
      TH1D* h1 = hInt1[plane_id][local_ch];
      TH1D* h0 = hInt0[plane_id][local_ch];
      if (!h1 || !h0) continue; // chmap 無い → スキップ

      const std::vector<double> &wf = (*waveform)[ch];
      const Int_t ns = (Int_t) wf.size();
      if (ns <= 0) continue;

      // ---- ベースライン（スライド平均の最小値）----
      Int_t startMaxEff = std::max(0, ns - BL_WIN);
      startMaxEff = std::min(startMaxEff, BL_START_MAX);

      double baseline = 0.0;
      {
        std::vector<double> bl_values;
        bl_values.reserve((startMaxEff / BL_STEP) + 2);
        for (Int_t start = 1; start <= startMaxEff; start += BL_STEP) {
          double sum = 0.0;
          const Int_t end = std::min(ns, start + BL_WIN);
          const Int_t width = end - start;
          if (width <= 0) continue;
          for (Int_t j = start; j < end; ++j) sum += wf[j];
          bl_values.push_back(sum / (double)width);
        }
        if (bl_values.empty()) continue;
        std::sort(bl_values.begin(), bl_values.end());
        baseline = bl_values.front();
      }

      // ---- 最大サンプル（可変窓の中心ピーク）----
      Int_t peak_idx = 0;
      double peak_val = -1e300;
      for (Int_t j = 0; j < ns; ++j) {
        if (wf[j] > peak_val) { peak_val = wf[j]; peak_idx = j; }
      }

      // ---- 可変窓の積分（ピーク中心 ±(INT_LEFT, INT_RIGHT)）----
      {
        Int_t jmin, jmax;
        if (peak_idx < INT_LEFT) {
          jmin = 0;
          jmax = std::min(ns - 1, INT_LEFT + INT_RIGHT - 1);
        } else if (peak_idx > ns - 1 - INT_RIGHT) {
          jmin = std::max(0, ns - (INT_LEFT + INT_RIGHT));
          jmax = ns - 1;
        } else {
          jmin = peak_idx - INT_LEFT;
          jmax = peak_idx + INT_RIGHT;
          if (jmin < 0) jmin = 0;
          if (jmax > ns - 1) jmax = ns - 1;
        }

        double integral1 = 0.0;
        for (Int_t j = jmin; j < jmax; ++j) integral1 += (wf[j] - baseline);
        if (peak_idx > INT_LEFT && peak_idx < ns - 1 - INT_RIGHT) {
          h1->Fill(integral1);
        }
      }

      // ---- 固定窓の積分（0 .. FIX_RANGE-1）----
      {
        const Int_t jmax_fix = std::min(ns, FIX_START + FIX_RANGE);
        if (jmax_fix > FIX_START) {
          double integral0 = 0.0;
          for (Int_t j = FIX_START; j < jmax_fix; ++j) integral0 += (wf[j] - baseline);
          h0->Fill(integral0);
        }
      }
    }
  }

  // ========= フィット結果の保管（後でCSVにまとめて出力） =========
  std::vector<std::vector<double>> mean0(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean0e(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean1(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean1e(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));

  std::vector<std::vector<int>> npeak(N_RAYRAW, std::vector<int>(N_CH_PER_PLANE, NAN)); //the number of peaks found by TSpectrum

  // ---------- プロット＆フィット（0pe：固定窓ヒストの最大ビン±10） ----------
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0);
  gStyle->SetStatW(0.22);
  gStyle->SetStatH(0.16);

  std::vector<TCanvas*> canv0(N_RAYRAW, (TCanvas*)nullptr);
  for (Int_t p = 0; p < N_RAYRAW; ++p) {
    canv0[p] = new TCanvas(Form("c0_%d", p+1), Form("RAYRAW#%d FIX", p+1), 2000, 1400);
    canv0[p]->Divide(8, 4);
    for (Int_t lc = 0; lc < N_CH_PER_PLANE; ++lc) {
      TH1D *h = hInt0[p][lc];
      canv0[p]->cd(lc + 1);
      if (!h || h->GetEntries() <= 0) continue;

      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.9);
      h->GetYaxis()->SetTitleOffset(0.9);
      h->SetFillColor(5);
      h->Draw();

      const double peakX0 = h->GetBinCenter(h->GetMaximumBin());
      const double fitmin0 = peakX0 - 10.0;
      const double fitmax0 = peakX0 + 10.0;

      TF1* f0 = new TF1(Form("f0_p%02d_ch%02d", p+1, lc), "gaus", fitmin0, fitmax0);
      f0->SetLineColor(kRed+1);
      TH1D* hcopy0 = (TH1D*)h->Clone(Form("h0fit_p%02d_ch%02d", p+1, lc));
      hcopy0->Fit(f0, "Q", "", fitmin0, fitmax0);

      h->Draw();
      f0->Draw("same");

      // 結果をパッドに描く & 保管
      double c  = f0->GetParameter(0), ce = f0->GetParError(0);
      double m  = f0->GetParameter(1), me = f0->GetParError(1);
      double s  = f0->GetParameter(2), se = f0->GetParError(2);

      TLatex *latex = new TLatex();
      latex->SetNDC(); latex->SetTextSize(0.05); latex->SetTextColor(kRed+1);
      latex->DrawLatex(0.45, 0.85, Form("Const = %.2f #pm %.2f", c, ce));
      latex->DrawLatex(0.45, 0.80, Form("Mean  = %.2f #pm %.2f", m, me));
      latex->DrawLatex(0.45, 0.75, Form("Sigma = %.2f #pm %.2f", s, se));

      mean0[p][lc]  = m;
      mean0e[p][lc] = me;
    }
    canv0[p]->Update();
    canv0[p]->SaveAs(Form("%s_RAYRAW#%02d_0pe.pdf", outputplot.Data(), p+1));
  }

  // ---------- プロット＆フィット（1pe：ピークサーチ→右側ピーク±10） ----------
  std::vector<TCanvas*> canv1(N_RAYRAW, (TCanvas*)nullptr);
  for (Int_t p = 0; p < N_RAYRAW; ++p) {
    canv1[p] = new TCanvas(Form("c1_%d", p+1), Form("RAYRAW#%d VAR", p+1), 2000, 1400);
    canv1[p]->Divide(8, 4);

    for (Int_t lc = 0; lc < N_CH_PER_PLANE; ++lc) {
      TH1D *h = hInt1[p][lc];
      canv1[p]->cd(lc + 1);
      if (!h || h->GetEntries() <= 0) continue;

      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.9);
      h->GetYaxis()->SetTitleOffset(0.9);
      h->SetFillColor(5);
      h->Draw();

      // TSpectrum で2ピーク探し、右側ピークを採用
      TSpectrum spec(PEAKS_MAX);
      Int_t nfound = spec.Search(h, SPEC_SIGMA, "", SPEC_THRESH);
      npeak[p][lc] = nfound;

      std::vector<double> xs;
      for (Int_t i = 0; i < nfound && i < 2; ++i) xs.push_back(spec.GetPositionX()[i]);
      std::sort(xs.begin(), xs.end()); // x昇順

      double peakX1 = xs.empty() ? h->GetBinCenter(h->GetMaximumBin())
                                 : (xs.size() >= 2 ? xs[1] : xs[0]);

      const double fitmin1 = peakX1 - 15.0;
      const double fitmax1 = peakX1 + 15.0;

      TF1* f1 = new TF1(Form("f1_p%02d_ch%02d", p+1, lc), "gaus", fitmin1, fitmax1);
      f1->SetLineColor(kRed+1);
      TH1D* hcopy1 = (TH1D*)h->Clone(Form("h1fit_p%02d_ch%02d", p+1, lc));
      hcopy1->Fit(f1, "Q", "", fitmin1, fitmax1);

      h->Draw();
      f1->Draw("same");

      // 描画 & 保管
      double c  = f1->GetParameter(0), ce = f1->GetParError(0);
      double m  = f1->GetParameter(1), me = f1->GetParError(1);
      double s  = f1->GetParameter(2), se = f1->GetParError(2);

      TLatex *latex = new TLatex();
      latex->SetNDC(); latex->SetTextSize(0.05); latex->SetTextColor(kRed+1);
      latex->DrawLatex(0.45, 0.85, Form("Const = %.2f #pm %.2f", c, ce));
      latex->DrawLatex(0.45, 0.80, Form("Mean  = %.2f #pm %.2f", m, me));
      latex->DrawLatex(0.45, 0.75, Form("Sigma = %.2f #pm %.2f", s, se));

      mean1[p][lc]  = m;
      mean1e[p][lc] = me;
    }
    canv1[p]->Update();
    canv1[p]->SaveAs(Form("%s_RAYRAW#%02d_1pe.pdf", outputplot.Data(), p+1));
  }

  // ---------- CSV 出力 ----------
  std::ofstream fout(outputcsv);
  if (!fout) {
    std::cerr << "[ERROR] CSV を作成できません: " << outputcsv.Data() << std::endl;
  } else {
    // 新フォーマット
    fout << "#cablenum,0pe,0pe_error,1pe,1pe_error,integralrange,badflag\n";
    for (Int_t p = 0; p < N_RAYRAW; ++p) {
      for (Int_t lc = 0; lc < N_CH_PER_PLANE; ++lc) {
        int cab = cabOf(p+1, lc);
        if (cab < 0) continue;

        double m0  = mean0[p][lc];
        double m0e = mean0e[p][lc];
        double m1v = mean1[p][lc];
        double m1e = mean1e[p][lc];

        int badflag = 0;
        if(m0>m1v || m1e>3.0 || npeak[p][lc] < 2 || m1v > 100) badflag=1; // badflag=0: good, 1: bad (due to missing ADC problem of RAYRAW)

        // NAN を空欄にしたい場合は適宜処理
        fout << cab << ","
             << m0  << "," << m0e << ","
             << m1v << "," << m1e << ","
             << FIX_RANGE << "," << badflag << "\n";
      }
    }
    fout.close();
    std::cout << "[INFO] CSV を保存しました: " << outputcsv.Data() << std::endl;
  }
}

// 呼び出し用ラッパ
void calibrayraw(const char* infile, const char* chmap_file)
{
  std::string path_cosmictest = "/group/nu/ninja/work/otani/FROST_cosmictest_dataandoutput/";
  TString filename_inroot = TString(path_cosmictest) + "data/rayraw/" + infile + ".root";
  TString filename_outcsv = TString(path_cosmictest) + "calibration/rayraw/calibresult/rayrawcalib_" + infile + ".csv";
  TString filename_outplot = TString(path_cosmictest) + "calibration/rayraw/plot/rayrawcalib_" + infile;
  TString filename_chmap = TString(path_cosmictest) + "chmap/" + chmap_file + ".txt";

  calibrayraw_(filename_inroot.Data(), filename_chmap.Data(), filename_outcsv.Data(), filename_outplot.Data());
}
