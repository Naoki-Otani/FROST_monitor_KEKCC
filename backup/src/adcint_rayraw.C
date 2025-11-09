// adcint_rayraw.C
// 入力 tree: evnum/I, unixtime/D, waveform(vector<vector<double>>),
//            leading(vector<vector<double>>), trailing(vector<vector<double>>), plane(vector<int> 任意)
// 出力 tree: cablenum[272]/I, rayraw[272]/I, rayrawch[272]/I,
//            adcint[272]/D, unixtime[272]/D,
//            leading(vector<vector<double>>), trailing(vector<vector<double>>)

//root -l -q 'adcint_rayraw.C("input.root", "chmap.txt")'


#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <iostream>

static const int INT_MIN_SAMP = 15;
static const int INT_MAX_SAMP = 90;   // 両端含む (j <= 90)
static const int BL_WIN       = 15;
static const int BL_START_MAX = 70;
static const int BL_STEP      = 5;

static const int CH_PER_PLANE = 32;
static const int NOUT         = 272;

// cablenum の配列順: 1..140, 201..332（計272）
static std::vector<int> make_cab_order() {
  std::vector<int> v; v.reserve(NOUT);
  for (int c = 1; c <= 140; ++c) v.push_back(c);
  for (int c = 201; c <= 332; ++c) v.push_back(c);
  return v;
}

struct ChMaps { //rr:rayraw#, lc:local ch#
  std::map<long long, int> rrLc_to_cab;           // key(rr,lc)→cab
  std::map<int, std::pair<int,int>> cab_to_rrLc;  // cab→(rr,lc)
};

static inline long long key_rr_lc(int rr, int lc) {
  return ( (long long)rr << 32 ) | (unsigned)lc;
}

static ChMaps load_chmap(const char* fname) {
  ChMaps m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_chmap", "chmap.txt を開けませんでした: %s", fname);
    return m;
  }
  std::string line; bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; } // 先頭ヘッダ行を保険で無視
    std::istringstream iss(line);
    int cab, rr, lc;
    if (!(iss >> cab >> rr >> lc)) continue;
    m.rrLc_to_cab[key_rr_lc(rr, lc)] = cab;
    m.cab_to_rrLc[cab] = {rr, lc};
  }
  return m;
}

static inline int cabOf(const ChMaps& maps, int rr, int lc) {
  auto it = maps.rrLc_to_cab.find(key_rr_lc(rr, lc));
  return (it == maps.rrLc_to_cab.end()) ? -1 : it->second;
}

static int index_from_cab(const std::vector<int>& order, int cab) {
  static std::map<int,int> cache;
  auto it = cache.find(cab);
  if (it != cache.end()) return it->second;
  auto p = std::find(order.begin(), order.end(), cab);
  int idx = (p == order.end()) ? -1 : (int)std::distance(order.begin(), p);
  cache[cab] = idx;
  return idx;
}

static double estimate_baseline(const std::vector<double>& wf) {
  const int ns = (int)wf.size();
  if (ns <= 0) return 0.0;
  const int startMaxEff = std::min(std::max(0, ns - BL_WIN), BL_START_MAX);

  double best = std::numeric_limits<double>::infinity();
  for (int start = 0; start <= startMaxEff; start += BL_STEP) {
    const int end   = std::min(ns, start + BL_WIN);
    const int width = end - start;
    if (width <= 0) continue;
    double sum = 0.0;
    for (int j = start; j < end; ++j) sum += wf[j];
    const double ave = sum / (double)width;
    if (ave < best) best = ave;
  }
  return std::isfinite(best) ? best : 0.0;
}

static double integrate_adc_minus_baseline(const std::vector<double>& wf, double baseline) {
  const int ns = (int)wf.size();
  if (ns <= INT_MIN_SAMP) return 0.0;
  const int jmin = INT_MIN_SAMP;
  const int jmax = std::min(INT_MAX_SAMP, ns - 1);
  if (jmin > jmax) return 0.0;
  double integ = 0.0;
  for (int j = jmin; j <= jmax; ++j) {
    integ += (wf[j] - baseline);
  }
  return integ;
}

void adcint_rayraw_(const char* infile, const char* outfile, const char* chmap_file)
{
  const auto CAB_ORDER = make_cab_order();
  const auto maps      = load_chmap(chmap_file);

  // 入力
  TFile fin(infile, "READ");
  if (fin.IsZombie()) { ::Error("make_adcint_from_baseline", "Open failed: %s", infile); return; }
  TTree* tin = (TTree*)fin.Get("tree");
  if (!tin) { ::Error("make_adcint_from_baseline", "TTree 'tree' not found"); return; }

  tin->SetBranchStatus("*", 0);
  int evnum = 0; double unixtime = 0.0;
  std::vector<std::vector<double>> *waveform = nullptr;
  std::vector<std::vector<double>> *leading  = nullptr;
  std::vector<std::vector<double>> *trailing = nullptr;
  std::vector<int> *plane = nullptr;

  tin->SetBranchStatus("evnum", 1);      tin->SetBranchAddress("evnum", &evnum);
  tin->SetBranchStatus("unixtime", 1);   tin->SetBranchAddress("unixtime", &unixtime);
  tin->SetBranchStatus("waveform", 1);   tin->SetBranchAddress("waveform", &waveform);
  if (tin->GetBranch("leading"))  { tin->SetBranchStatus("leading", 1);  tin->SetBranchAddress("leading",  &leading);  }
  if (tin->GetBranch("trailing")) { tin->SetBranchStatus("trailing", 1); tin->SetBranchAddress("trailing", &trailing); }
  if (tin->GetBranch("plane"))    { tin->SetBranchStatus("plane", 1);    tin->SetBranchAddress("plane",    &plane);    }

  // 出力
  TFile fout(outfile, "RECREATE");
  TTree tout("tree", "ADC integral (15..90) per cablenum");

  int    out_cablenum[NOUT];
  int    out_rayraw[NOUT];
  int    out_rayrawch[NOUT];
  double out_adcint[NOUT];
  double out_unixtime[NOUT];
  auto   out_leading  = new std::vector<std::vector<double>>();
  auto   out_trailing = new std::vector<std::vector<double>>();

  tout.Branch("cablenum", out_cablenum,  Form("cablenum[%d]/I", NOUT));
  tout.Branch("rayraw",   out_rayraw,    Form("rayraw[%d]/I",   NOUT));
  tout.Branch("rayrawch", out_rayrawch,  Form("rayrawch[%d]/I", NOUT));
  tout.Branch("adcint",   out_adcint,    Form("adcint[%d]/D",   NOUT));
  tout.Branch("unixtime", out_unixtime,  Form("unixtime[%d]/D", NOUT));
  tout.Branch("leading",  &out_leading);
  tout.Branch("trailing", &out_trailing);

  // 固定情報（cablenum / rayraw / rayrawch）
  std::vector<int> fixed_rayraw(NOUT, -1), fixed_rayrawch(NOUT, -1);
  for (int i = 0; i < NOUT; ++i) {
    const int cab = CAB_ORDER[i];
    out_cablenum[i] = cab;
    auto it = maps.cab_to_rrLc.find(cab);
    if (it != maps.cab_to_rrLc.end()) {
      fixed_rayraw[i]   = it->second.first;   // RAYRAW (1始まり)
      fixed_rayrawch[i] = it->second.second;  // local ch (0..31)
    }
  }

  const Long64_t nent = tin->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tin->GetEntry(ie);

    // 初期化
    out_leading->assign(NOUT, std::vector<double>());   // 空ベクタで NOUT 分
    out_trailing->assign(NOUT, std::vector<double>());

    for (int i = 0; i < NOUT; ++i) {
      out_rayraw[i]   = fixed_rayraw[i];
      out_rayrawch[i] = fixed_rayrawch[i];
      out_adcint[i]   = std::numeric_limits<double>::quiet_NaN();
      out_unixtime[i] = unixtime;  // イベントのスカラを複製
    }

    if (waveform) {
      const int nChAll = (int)waveform->size();
      for (int ch = 0; ch < nChAll; ++ch) {
        // plane があれば使用、無ければ ch/32 で推定
        int plane_id = (plane && (int)plane->size() > ch) ? plane->at(ch) : (ch / CH_PER_PLANE);
        if (plane_id < 0) continue;
        const int rr = plane_id + 1;         // RAYRAW を 1 始まりに
        const int lc = ch % CH_PER_PLANE;

        const int cab = cabOf(maps, rr, lc);
        if (cab < 0) continue;
        const int idx = index_from_cab(CAB_ORDER, cab);
        if (idx < 0 || idx >= NOUT) continue;

        // baseline → 積分
        const auto& wf = waveform->at(ch);
        if ((int)wf.size() > INT_MIN_SAMP) {
          const double bl    = estimate_baseline(wf);
          const double integ = integrate_adc_minus_baseline(wf, bl);
          out_adcint[idx] = integ;
        }

        // leading/trailing: 入力の該当 ch ベクタをそのままコピー（無ければ空のまま）
        if (leading && (int)leading->size() > ch)  (*out_leading)[idx]  = leading->at(ch);
        if (trailing && (int)trailing->size() > ch) (*out_trailing)[idx] = trailing->at(ch);
      }
    }

    tout.Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "[OK] wrote: " << outfile << std::endl;
}

// ROOT から関数形式で呼びやすいよう wrapper を用意
void adcint_rayraw(const char* infile, const char* chmap_file)
{
  std::string path_cosmictest = "/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/";
  TString filename_inroot = TString(path_cosmictest) + "data/rayraw/" + infile + ".root";
  TString filename_outroot = TString(path_cosmictest) + "dataaftercalib/rayraw/adcintegral/" + infile + "_adcintegral.root";
  TString filename_chmap = TString(path_cosmictest) + "chmap/" + chmap_file + ".txt";
  adcint_rayraw_(filename_inroot.Data(), filename_outroot.Data(), filename_chmap.Data());
}
