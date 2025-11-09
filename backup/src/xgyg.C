// xgyg.C
// 使い方：
//   root -l -b -q 'XYBarycenter.C+("run00044")'
// 入出力：
//   入力:  /.../lightyield/<run>_lightyield.root  の TTree "tree"（配列 272）
//   出力:  PDF: LightyieldHisto/xy_<run>.pdf
//          ROOT(optional): LightyieldHisto/xy_<run>.root に TTree "xy"（xg, yg）

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <limits>
#include <cmath>
#include <string>

static const int NOUT    = 272;   // 配列長（イベントあたり）
static inline bool isfinite_d(double x){ return std::isfinite(x); }

// cablenum -> 位置(mm) 対応
// 戻り値: true=有効なケーブル、false=対象外
inline bool CableToPosition(int cablenum, double& x, double& y) {
  x = std::numeric_limits<double>::quiet_NaN();
  y = std::numeric_limits<double>::quiet_NaN();

  if (1 <= cablenum && cablenum <= 140) {
    // y: 695, 685, ..., -695  (10mm刻み)
    y = 695.0 - 10.0 * (cablenum - 1);
    return true;
  }
  if (201 <= cablenum && cablenum <= 332) {
    // x: 655, 645, ..., -655 (10mm刻み)
    x = 655.0 - 10.0 * (cablenum - 201);
    return true;
  }
  return false; // その他のケーブルは重心計算の対象外
}

void xgyg_(const char* infile, const char* outprefix, double xgweight,
                   bool write_xy_root = true)
{
  TH1::AddDirectory(kFALSE);

  std::ostringstream oss_xgweight;
  oss_xgweight << std::fixed << std::setprecision(1) << xgweight;
  std::string s_xgweight = oss_xgweight.str();        // "2.0"

  // 入力を開く
  TFile f(infile, "READ");
  if (f.IsZombie()) { ::Error("XYBarycenter","Open failed: %s", infile); return; }
  TTree* t = (TTree*)f.Get("tree");
  if (!t) { ::Error("XYBarycenter","TTree 'tree' not found"); return; }

  Int_t    cablenum[NOUT];
  Int_t    rayraw[NOUT];
  Int_t    rayrawch[NOUT];
  Double_t lightyield[NOUT];
  Double_t unixtime[NOUT];

  std::vector<std::vector<double>>* leading  = nullptr;
  std::vector<std::vector<double>>* trailing = nullptr;

  // --- SetBranchAddress for input ---
  t->SetBranchAddress("cablenum",   cablenum);
  t->SetBranchAddress("rayraw",     rayraw);
  t->SetBranchAddress("rayrawch",   rayrawch);
  t->SetBranchAddress("lightyield", lightyield);
  t->SetBranchAddress("unixtime",   unixtime);
  t->SetBranchAddress("leading",   &leading);
  t->SetBranchAddress("trailing",  &trailing);

  // 2D ヒスト（ビン数はケーブル本数に合わせる）
  //   x: -655..655 (132本)  y: -695..695 (140本)
  TH2D* hXY = new TH2D("h", "; x_{g} [mm]; y_{g} [mm]",
                       132, -660.0, 660.0,
                       140, -700.0, 700.0);

  // 必要なら xg/yg を出力する TTree（オプション）
  TFile* fout = nullptr;
  TTree* tout = nullptr;
  double xg = std::numeric_limits<double>::quiet_NaN();
  double yg = std::numeric_limits<double>::quiet_NaN();
  double lightsum_x = std::numeric_limits<double>::quiet_NaN();
  double lightsum_y = std::numeric_limits<double>::quiet_NaN();
  if (write_xy_root) {
    TString outroot = TString::Format("%s_weight%s.root", outprefix,s_xgweight.c_str());
    fout = TFile::Open(outroot, "RECREATE");
    tout = new TTree("tree", "lightyield barycenter per event");
    tout->Branch("xg", &xg, "xg/D");
    tout->Branch("yg", &yg, "yg/D");
    tout->Branch("lightsum_x", &lightsum_x, "lightsum_x/D");
    tout->Branch("lightsum_y", &lightsum_y, "lightsum_y/D");

     // --- Define output branches with the same layout ---
    tout->Branch("cablenum",   cablenum,   Form("cablenum[%d]/I", NOUT));
    tout->Branch("rayraw",     rayraw,     Form("rayraw[%d]/I", NOUT));
    tout->Branch("rayrawch",   rayrawch,   Form("rayrawch[%d]/I", NOUT));
    tout->Branch("lightyield", lightyield, Form("lightyield[%d]/D", NOUT));
    tout->Branch("unixtime",   unixtime,   Form("unixtime[%d]/D", NOUT));

    // vector<vector<double>> はポインタのアドレスを渡す
    // （GetEntry後に先頭が有効なデータを指すので、その内容がそのまま書かれます）
    tout->Branch("leading",  &leading);
    tout->Branch("trailing", &trailing);
  }

  // イベントループ
  const Long64_t nent = t->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);

    // 分子・分母
    double x_num = 0.0, x_den = 0.0;
    double y_num = 0.0, y_den = 0.0;
    lightsum_x = 0.0, lightsum_y = 0.0;

    for (int i = 0; i < NOUT; ++i) {
      const double ly = lightyield[i];
      if (!(isfinite_d(ly) && ly > 0.0)) continue; // NaN/Inf/<=0 を除外

      double x_i, y_i;
      if (!CableToPosition(cablenum[i], x_i, y_i)) continue;

      if (isfinite_d(x_i)) { x_num += pow(ly,xgweight) * x_i; x_den += pow(ly,xgweight); lightsum_x += ly; }
      if (isfinite_d(y_i)) { y_num += pow(ly,xgweight) * y_i; y_den += pow(ly,xgweight); lightsum_y += ly; }
    }

    xg = (x_den > 0.0) ? (x_num / x_den) : std::numeric_limits<double>::quiet_NaN();
    yg = (y_den > 0.0) ? (y_num / y_den) : std::numeric_limits<double>::quiet_NaN();

    // 両方とも有限のときだけ 2D に Fill
    if (isfinite_d(xg) && isfinite_d(yg) && lightsum_x+lightsum_y>150) hXY->Fill(xg, yg);

    if (write_xy_root) tout->Fill();
  }

  // 描画
  // gStyle->SetOptStat(0);
  gStyle->SetStatX(0.99);       // 右上X位置（デフォルト0.99）
  gStyle->SetStatY(0.99);       // 右上Y位置
  gStyle->SetStatW(0.1);       // 幅（0〜1）
  gStyle->SetStatH(0.1);       // 高さ（0〜1）
  TCanvas* c = new TCanvas("c", "c", 1000, 900);
  hXY->SetContour(99);
  hXY->Draw("COLZ");
  hXY->GetYaxis()->SetTitleOffset(1.1);
  c->Modified(); c->Update();

  // PDF 保存
  TString outpdf = TString::Format("%s_weight%s.pdf", outprefix,s_xgweight.c_str());
  c->SaveAs(outpdf);

  // ルート出力
  if (write_xy_root) {
    fout->cd();
    hXY->Write(); // 2Dヒストも一緒に保存しておく
    tout->Write();
    fout->Close();
  }

  // 入力を閉じる（ヒストはファイル非依存のため消えない）
  f.Close();
}

// 呼び出しラッパ（あなたの既存パスに合わせています）
void xgyg(const char* infile, double weight)
{
  std::string base = "/group/nu/ninja/work/otani/FROST_cosmictest_dataandoutput/";
  TString filename_inroot = TString(base) + "dataaftercalib/rayraw/" + infile + "_lightyield.root";
  TString outprefix = TString(base) + "xgyg/xgyg_" + infile;

  xgyg_(filename_inroot.Data(), outprefix.Data(), weight, /*write_xy_root=*/true);
}
