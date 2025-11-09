// convert_lightyield.C
// root -l -q 'convert_lightyield.C("input.root", "calibration.csv", "output.root")'

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>

static const int NCH = 16;
static const int NSAMPLE = 150;   // 入力ブランチ waveform[*][150]/D を想定
static const int INTMIN = 25;     // 積分下限 (inclusive)
static const int INTMAX = 85;     // 積分上限 (exclusive) → 25..84 を積分
static const int BASE_WIN = 20;   // スライド平均の窓幅
static const int BASE_START_MAX = 70; // スライド開始 0,5,...,70
static const int BASE_STEP = 5;

// CSV は 1 行目ヘッダ、以降： name(or label), m0, m1, intrange
struct Calib {
  double m0[NCH] = {0.0};
  double m1[NCH] = {1.0};
  int intrange[NCH] = {1};
};

bool load_calib_csv(const std::string& csv, Calib& c){
  std::ifstream infile(csv);
  if(!infile.is_open()){
    std::cerr << "Error: cannot open calib csv: " << csv << std::endl;
    return false;
  }
  std::string line;
  if(!std::getline(infile, line)) { // header
    std::cerr << "Error: empty calib csv: " << csv << std::endl;
    return false;
  }
  int logical_idx = 0;
  while(std::getline(infile, line) && logical_idx < NCH){
    if(line.empty()) continue;
    std::stringstream ss(line);
    std::string item;
    int col = 0;
    int target = logical_idx; // 通常は行順
    std::string name;
    double m0=0, m1=1; int ir=1;
    while(std::getline(ss, item, ',')){
      if(col == 1){ m0 = std::stod(item);
      }else if(col == 2){ m1 = std::stod(item);
      }else if(col == 3){ ir = std::max(1, (int)std::stod(item));
      }
      ++col;
    }
    if(target < 0 || target >= NCH){
      std::cerr << "Warning: invalid target index parsed from csv row: " << line << std::endl;
    }else{
      c.m0[target] = m0;
      c.m1[target] = m1;
      c.intrange[target] = ir;
    }
    ++logical_idx;
  }
  return true;
}

// ベースライン推定
// スライド平均の最大
double estimate_baseline(const double* wf, int ch){
  std::vector<double> v;
  v.reserve((BASE_START_MAX/BASE_STEP)+1);
  for(int start=0; start<=BASE_START_MAX; start+=BASE_STEP){
    double sum=0.0;
    int end = start + BASE_WIN;
    if(end > NSAMPLE) end = NSAMPLE;
    for(int j=start; j<end; ++j) sum += wf[j];
    v.push_back(sum / std::max(1, end-start));
  }
  std::sort(v.begin(), v.end());
  return v.back();
}

// 1 イベント・1 チャンネルの積分（ベースライン差分）
double integrate_charge(const double* wf, int ch, double baseline){
  double value = 0.0;
  for(int j=INTMIN; j<INTMAX && j<NSAMPLE; ++j){
    value += (baseline - wf[j]);
  }
  return value;
}

// 電荷→p.e. 変換
double charge_to_pe(double value, const Calib& calib, int ch){
  const double m0 = calib.m0[ch];
  const double m1 = calib.m1[ch];
  const int ir = std::max(1, calib.intrange[ch]);
  const double width = (double)(INTMAX-INTMIN); // 積分サンプル数
  // 参照式: pe = ( value - m0 * width / intrange ) / (m1 - m0)
  const double denom = (m1 - m0);
  if(std::abs(denom) < 1e-12) return 0.0; // 安全策
  return ( value - m0 * (width / (double)ir) ) / denom;
}

void convert_lightyield(const char* in_root,
                        const char* calib_csv,
                        const char* out_root)
{
  // --- キャリブ読込 ---
  Calib calib;
  if(!load_calib_csv(calib_csv, calib)){
    std::cerr << "FATAL: calib CSV 読み込み失敗\n";
    return;
  }

  // --- 入力 ---
  TFile fin(in_root, "READ");
  if(fin.IsZombie()){ std::cerr << "FATAL: cannot open input root\n"; return; }
  TTree* tin = (TTree*)fin.Get("t");
  if(!tin){ std::cerr << "FATAL: tree 't' not found\n"; return; }

  // 必要ブランチのみ有効化（速度向上）
  tin->SetBranchStatus("*", 0);
  double wf[NCH][NSAMPLE] = {{0}};
  for(int ch=0; ch<NCH; ++ch){
    TString bname = Form("waveform%d", ch);
    if(tin->GetBranch(bname)){
      tin->SetBranchStatus(bname, 1);
      tin->SetBranchAddress(bname, wf[ch]);
    }else{
      std::cerr << "FATAL: missing branch " << bname << std::endl;
      return;
    }
  }
  double timestamp = 0.0;
  if(tin->GetBranch("timestamp")){
    tin->SetBranchStatus("timestamp", 1);
    tin->SetBranchAddress("timestamp", &timestamp);
  }else{
    std::cerr << "FATAL: missing branch timestamp\n";
    return;
  }

  // --- 出力 ---
  TFile fout(out_root, "RECREATE");
  TTree tout("t", "per-event light yield (p.e.) and timestamp");
  double lightyield[NCH] = {0.0};
  tout.Branch("lightyield", lightyield, Form("lightyield[%d]/D", NCH));
  tout.Branch("timestamp", &timestamp, "timestamp/D");

  // --- ループ ---
  const Long64_t nent = tin->GetEntries();
  for(Long64_t ie=0; ie<nent; ++ie){
    tin->GetEntry(ie);

    // 各 ch の p.e. を計算
    for(int ch=0; ch<NCH; ++ch){
      const double baseline = estimate_baseline(wf[ch], ch);
      const double q = integrate_charge(wf[ch], ch, baseline);
      lightyield[ch] = charge_to_pe(q, calib, ch);
    }

    tout.Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "Done. Wrote: " << out_root << std::endl;
}

// ROOT から関数形式で呼びやすいよう wrapper を用意
void convertlightyield_caen(const char* in_root,
                        const char* calib_csv)
{
  std::string path_cosmictest = "/group/nu/ninja/work/otani/FROST_cosmictest_dataandoutput/";
  TString filename_inroot = TString(path_cosmictest) + "data/caen/caen" + in_root + "wave.root";
  TString filename_scv = TString(path_cosmictest) + "calibration/caen/calibresult/caencalib" + calib_csv + ".csv";
  TString filename_outroot = TString(path_cosmictest) + "dataaftercalib/caen/caen" + in_root + "_aftercalib.root";
  convert_lightyield(filename_inroot, filename_scv, filename_outroot);
}
