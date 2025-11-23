// convertlightyield.cpp
// Standalone ROOT program to convert ADC integrals to light yield per channel/bunch.
// It scans calibration CSVs under calibration/calibresult/ and, for each run name X:
//   - requires rootfile/X.root to exist
//   - skips if rootfile_aftercalib/X_lightyield.root already exists
// Default: watch mode every 60 seconds (Ctrl-C to stop). Use --oneshot 1 for single scan.
//
// Build:
//   g++ -std=c++17 -O2 convertlightyield.cpp -o convertlightyield \
//       $(root-config --cflags --libs) -lSpectrum
//   # On older GCC add: -lstdc++fs
//
// Run (watch mode default 60s):
//   ./convertlightyield
//
// One-shot:
//   ./convertlightyield --oneshot 1
//
// Options:
//   -h, --help
//   --base <DIR>        (default /group/nu/ninja/work/otani/FROST_beamdata/test)
//   --chmap <FILE>      (default chmap_20251009.txt) under chmap/
//   --refgain <FILE>    (default refgain.csv) under calibration/ReferenceGain/
//   --watch <SEC>       (default 60; 0 disables watching)
//   --oneshot <0|1>     (default 0)

// ------- ROOT -------
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>

// ------- std -------
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <iostream>
#include <filesystem>
#include <regex>
#include <chrono>
#include <csignal>
#include <atomic>
#include <thread>
#include <cstdint>
#include <iomanip>

#include "/home/nu/notani/FROST_monitor/config/config.hpp"

namespace fs = std::filesystem;

// -----------------------------
// Parameters
// -----------------------------
static const Int_t INT_RANGE = FrostmonConfig::INT_RANGE;
static const Int_t BL_WIN       = FrostmonConfig::BL_WIN;
static const Int_t BL_START_MAX = FrostmonConfig::BL_START_MAX;
static const Int_t BL_STEP      = FrostmonConfig::BL_STEP;

static const Int_t SAMPLING_FIRST_BANCH = FrostmonConfig::SAMPLING_FIRST_BANCH;   // integration start index for bunch#0
static const Int_t ADC_THRESHOLD = FrostmonConfig::ADC_THRESHOLD; //threshold for leading/trailing_fromadc
static const Double_t BUNCH_INTERVAL = FrostmonConfig::BUNCH_INTERVAL;

static const Int_t CH_PER_PLANE = FrostmonConfig::N_CH_PER_PLANE;
static const Int_t NOUT         = FrostmonConfig::NOUT;
static const Int_t NBUNCH       = FrostmonConfig::NBUNCH;

static const Int_t BASEREF_CALIB_CABLENUM = FrostmonConfig::BASEREF_CALIB_CABLENUM;

static const Int_t SPILL_BIT_ADC_THRESHOLD = FrostmonConfig::SPILL_BIT_ADC_THRESHOLD;
static const Int_t SPILL_BIT_MIN_POINTS    = FrostmonConfig::SPILL_BIT_MIN_POINTS;

// -----------------------------
// Utilities & Models
// -----------------------------
static std::atomic<bool> g_stop{false};
static void handle_sigint(int){ g_stop = true; }

// file_time_type -> time_t (C++17-safe)
static std::time_t filetime_to_time_t(std::filesystem::file_time_type ftime) {
  using namespace std::chrono;
  auto sctp = time_point_cast<system_clock::duration>(
      ftime - std::filesystem::file_time_type::clock::now() + system_clock::now());
  return system_clock::to_time_t(sctp);
}

// cablenum order: 1..140, 201..332 (total 272)
static std::vector<int> make_cab_order() {
  std::vector<int> v; v.reserve(NOUT);
  for (int c = 1; c <= 140; ++c) v.push_back(c);
  for (int c = 201; c <= 332; ++c) v.push_back(c);
  return v;
}

struct ChMaps { // rr: rayraw#, lc: local ch#
  std::map<long long, int> rrLc_to_cab;           // key(rr,lc) → cab
  std::map<int, std::pair<int,int>> cab_to_rrLc;  // cab → (rr,lc)
};

static inline long long key_rr_lc(int rr, int lc) {
  return ( (long long)rr << 32 ) | (unsigned)lc;
}

struct LightYieldCorrMap {
  // cab → correction_factor
  std::map<int, double> cab_to_corr;
};

struct SpillBitMap {
  // key(rr,lc) → bit index (0-based; bit1→0, bit2→1, ...)
  std::map<long long, int> rrLc_to_bit0;
};

static SpillBitMap load_spillmap(const char* fname) {
  SpillBitMap m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_spillmap", "Cannot open spill chmap: %s", fname);
    return m;
  }
  std::string line;
  bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; } // ヘッダ行がある場合に備えて
    std::istringstream iss(line);
    int bit, rr, lc;
    if (!(iss >> bit >> rr >> lc)) continue;
    if (bit <= 0) continue;
    int bit0 = bit - 1; // bit1 → 0 (2^0), bit2 → 1 (2^1), ...
    m.rrLc_to_bit0[key_rr_lc(rr, lc)] = bit0;
  }
  return m;
}

static ChMaps load_chmap(const char* fname) {
  ChMaps m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_chmap", "Cannot open chmap: %s", fname);
    return m;
  }
  std::string line; bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; }
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

struct RefGainMap {
  // cab → (m0, m1)
  std::map<int, std::pair<double,double>> cab_to_m01_ref;
};

static RefGainMap load_calibref(const char* fname) {
  RefGainMap m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_calibref", "Cannot open reference gain: %s", fname);
    return m;
  }
  std::string line;
  bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; }
    for (char &c : line) if (c == ',') c = ' ';
    std::istringstream iss(line);
    int cab = -1;
    double m0 = 0.0, m1 = 0.0;
    if (!(iss >> cab >> m0 >> m1)) continue;
    m.cab_to_m01_ref[cab] = std::make_pair(m0, m1);
  }
  return m;
}

struct RayrawCalibMap {
  // cab → (m0, m1)
  std::map<int, std::pair<double,double>> cab_to_m01;
  std::map<int, int> cab_to_intrange;
};

static RayrawCalibMap load_calib(const char* fname) {
  RayrawCalibMap m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_calib", "Cannot open calibration CSV: %s", fname);
    return m;
  }
  std::string line;
  bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; }
    for (char &c : line) if (c == ',') c = ' ';
    std::istringstream iss(line);
    int cab = -1, integral = 0, bad = 0;
    double m0 = 0.0, m0e = 0.0, m1 = 0.0, m1e = 0.0;
    if (!(iss >> cab >> m0 >> m0e >> m1 >> m1e >> integral)) continue;
    if (iss >> bad) { /* optional bad flag */ }
    m.cab_to_m01[cab] = std::make_pair(m0, m1);
    m.cab_to_intrange[cab] = integral;
  }
  return m;
}

static LightYieldCorrMap load_lightyield_corr(const char* fname) {
  LightYieldCorrMap m;
  std::ifstream fin(fname);
  if (!fin) {
    ::Warning("load_lightyield_corr", "Cannot open lightyield correction CSV: %s", fname);
    return m;
  }
  std::string line;
  bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { first = false; continue; }
    if (first) { first = false; continue; }
    for (char &c : line) if (c == ',') c = ' ';
    std::istringstream iss(line);
    int cab = -1;
    double corr = 1.0;
    if (!(iss >> cab >> corr)) continue;
    m.cab_to_corr[cab] = corr;
  }
  return m;
}

static double estimate_baseline(const std::vector<double>& wf) {
  const int ns = (int)wf.size();
  if (ns <= 0) return 0.0;
  const int startMaxEff = std::min(std::max(0, ns - BL_WIN), BL_START_MAX);

  double best = std::numeric_limits<double>::infinity();
  for (int start = 1; start <= startMaxEff; start += BL_STEP) {
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

static double integrate_adc_minus_baseline(const std::vector<double>& wf, double baseline, int int_min) {
  const int ns = (int)wf.size();
  if (ns <= INT_RANGE) return 0.0;
  const int jmin = int_min;
  const int jmax = std::min(int_min + INT_RANGE, ns - 1);
  if (jmin > jmax) return 0.0;
  double integ = 0.0;
  for (int j = jmin; j < jmax; ++j) {
    integ += (wf[j] - baseline);
  }
  return integ;
}

static bool spillbit_from_waveform(const std::vector<double>& wf) {
  int count = 0;
  for (double v : wf) {
    if (v >= SPILL_BIT_ADC_THRESHOLD) {
      ++count;
      if (count >= SPILL_BIT_MIN_POINTS) return true;
    }
  }
  return false;
}

static std::vector<double> calculate_leading_fromadc(const std::vector<double>& wf) {
  std::vector<double> leading_adc;
  const int ns = (int)wf.size();
  for(int j=0; j<ns-1; ++j){
    if(wf[j] < ADC_THRESHOLD && wf[j+1] > ADC_THRESHOLD){
      leading_adc.push_back((double)(j+1));
    }
  }
  return leading_adc;
}

static std::vector<double> calculate_trailing_fromadc(const std::vector<double>& wf) {
  std::vector<double> trailing_adc;
  const int ns = (int)wf.size();
  for(int j=0; j<ns-1; ++j){
    if(wf[j] > ADC_THRESHOLD && wf[j+1] < ADC_THRESHOLD){
      trailing_adc.push_back((double)(j+1));
    }
  }
  return trailing_adc;
}

// -----------------------------
// Core converter (single run)
// -----------------------------
static void convertlightyield_rayraw_(const char* infile, const char* outfile,
                                      const char* chmap_file,
                                      const char* referencegain_csv,
                                      const char* rayrawcalib_csv,
                                      const char* lycorr_csv,
                                      const char* spillchmap_file)
{
  const auto CAB_ORDER = make_cab_order();
  const auto maps      = load_chmap(chmap_file);
  const auto refgain   = load_calibref(referencegain_csv);
  const auto calib     = load_calib(rayrawcalib_csv);
  const auto spillmap  = load_spillmap(spillchmap_file);
  const auto lycorr    = load_lightyield_corr(lycorr_csv);

  // Input ROOT
  TFile fin(infile, "READ");
  if (fin.IsZombie()) { ::Error("convertlightyield", "Open failed: %s", infile); return; }
  TTree* tin = (TTree*)fin.Get("tree");
  if (!tin) { ::Error("convertlightyield", "TTree 'tree' not found"); return; }

  tin->SetBranchStatus("*", 0);
  int evnum = 0; double unixtime = 0.0;
  std::vector<std::vector<double>> *waveform = nullptr;
  std::vector<std::vector<double>> *leading  = nullptr;
  std::vector<std::vector<double>> *trailing = nullptr;
  std::vector<int> *plane = nullptr;

  tin->SetBranchStatus("evnum", 1);      tin->SetBranchAddress("evnum", &evnum);
  tin->SetBranchStatus("unixtime", 1);   tin->SetBranchAddress("unixtime", &unixtime);
  tin->SetBranchStatus("waveform", 1);   tin->SetBranchAddress("waveform", &waveform);
  tin->SetBranchStatus("leading", 1);    tin->SetBranchAddress("leading",  &leading);
  tin->SetBranchStatus("trailing", 1);   tin->SetBranchAddress("trailing", &trailing);
  tin->SetBranchStatus("plane", 1);      tin->SetBranchAddress("plane",    &plane);

  // Output ROOT
  TFile fout(outfile, "RECREATE");
  TTree tout("tree", "tree after calibration");

  int    out_evnum = 0;
  int    out_cablenum[NOUT];
  int    out_rayraw[NOUT];
  int    out_rayrawch[NOUT];
  double out_lightyield[NOUT][NBUNCH];
  double out_unixtime[NOUT];
  int out_spillnum;
  auto   out_leading  = new std::vector<std::vector<double>>();
  auto   out_trailing = new std::vector<std::vector<double>>();
  auto   out_leading_fromadc  = new std::vector<std::vector<double>>();
  auto   out_trailing_fromadc = new std::vector<std::vector<double>>();

  tout.Branch("evnum", &out_evnum,  "evnum/I");
  tout.Branch("cablenum", out_cablenum,  Form("cablenum[%d]/I", NOUT));
  tout.Branch("rayraw",   out_rayraw,    Form("rayraw[%d]/I",   NOUT));
  tout.Branch("rayrawch", out_rayrawch,  Form("rayrawch[%d]/I", NOUT));
  tout.Branch("lightyield", out_lightyield, Form("lightyield[%d][%d]/D",   NOUT, NBUNCH));
  tout.Branch("unixtime", out_unixtime,  Form("unixtime[%d]/D", NOUT));
  tout.Branch("spillnum", &out_spillnum,  "spillnum/I");
  tout.Branch("leading",  &out_leading);
  tout.Branch("trailing", &out_trailing);
  tout.Branch("leading_fromadc",  &out_leading_fromadc);
  tout.Branch("trailing_fromadc", &out_trailing_fromadc);

  // Fixed channel info
  std::vector<int> fixed_rayraw(NOUT, -1), fixed_rayrawch(NOUT, -1);
  for (int i = 0; i < NOUT; ++i) {
    const int cab = CAB_ORDER[i];
    out_cablenum[i] = cab;
    auto it = maps.cab_to_rrLc.find(cab);
    if (it != maps.cab_to_rrLc.end()) {
      fixed_rayraw[i]   = it->second.first;   // RAYRAW (1-based)
      fixed_rayrawch[i] = it->second.second;  // local ch (0..31)
    }
  }

  const Long64_t nent = tin->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tin->GetEntry(ie);

    out_leading->assign(NOUT, std::vector<double>());
    out_trailing->assign(NOUT, std::vector<double>());
    out_leading_fromadc->assign(NOUT, std::vector<double>());
    out_trailing_fromadc->assign(NOUT, std::vector<double>());

    for (int i = 0; i < NOUT; ++i) {
      out_rayraw[i]   = fixed_rayraw[i];
      out_rayrawch[i] = fixed_rayrawch[i];
      out_unixtime[i] = unixtime;
      for(int bunch=0; bunch<NBUNCH; ++bunch){
        out_lightyield[i][bunch] = std::numeric_limits<double>::quiet_NaN();
      }
    }

    out_evnum = evnum;

    int spill_bits[15] = {0};

    if (waveform) {
      const int nChAll = (int)waveform->size();
      for (int ch = 0; ch < nChAll; ++ch) {
        int plane_id = (plane && (int)plane->size() > ch) ? plane->at(ch) : (ch / CH_PER_PLANE);
        if (plane_id < 0) continue;
        const int rr = plane_id + 1;
        const int lc = ch % CH_PER_PLANE;

        if (!spillmap.rrLc_to_bit0.empty()) {
          auto itb = spillmap.rrLc_to_bit0.find(key_rr_lc(rr, lc));
          if (itb != spillmap.rrLc_to_bit0.end()) {
            int bit0 = itb->second;
            if (bit0 >= 0 && bit0 < 15) {
              const auto& wf_spill = waveform->at(ch);
              bool is_one = spillbit_from_waveform(wf_spill);
              spill_bits[bit0] = is_one ? 1 : 0;
            }
          }
        }

        const int cab = cabOf(maps, rr, lc);
        if (cab < 0) continue;
        const int idx = index_from_cab(CAB_ORDER, cab);
        if (idx < 0 || idx >= NOUT) continue;

        double adcint[NBUNCH] = {0.0};
        const auto& wf = waveform->at(ch);
        if ((int)wf.size() > SAMPLING_FIRST_BANCH) {
          const double bl    = estimate_baseline(wf);
          for(int bunch=0; bunch<NBUNCH; ++bunch) {
            const int int_min = SAMPLING_FIRST_BANCH + (int)(bunch * BUNCH_INTERVAL);
            adcint[bunch] = integrate_adc_minus_baseline(wf, bl, int_min);
          }
        }

        // Fetch calibration values (with basic guard/defaults)
        double zero_pe = 50.0, one_pe_base = 0.0, zero_pe_base = 0.0;
        int intrange_calib = INT_RANGE;

        if (auto it = calib.cab_to_m01.find(cab); it != calib.cab_to_m01.end()) {
          zero_pe = it->second.first;
          if(zero_pe < 5.0 || zero_pe > 100.0) zero_pe = 50.0;
        }
        if (auto it = calib.cab_to_intrange.find(cab); it != calib.cab_to_intrange.end()) {
          intrange_calib = it->second;
          if (intrange_calib <= 0) intrange_calib = INT_RANGE;
        }
        if (auto it = calib.cab_to_m01.find(BASEREF_CALIB_CABLENUM); it != calib.cab_to_m01.end()) {
          zero_pe_base = it->second.first;
          one_pe_base  = it->second.second;
        }

        double m0_ref = 0.0, m1_ref = 1.0, m0_ref_base = 0.0, m1_ref_base = 1.0;
        if (auto it = refgain.cab_to_m01_ref.find(cab); it != refgain.cab_to_m01_ref.end()) {
          m0_ref = it->second.first; m1_ref = it->second.second;
        }
        if (auto it = refgain.cab_to_m01_ref.find(BASEREF_CALIB_CABLENUM); it != refgain.cab_to_m01_ref.end()) {
          m0_ref_base = it->second.first; m1_ref_base = it->second.second;
        }

        double corr_factor = 1.0;
        if (auto it = lycorr.cab_to_corr.find(cab); it != lycorr.cab_to_corr.end()) {
          corr_factor = it->second;
        }

        for(int bunch=0; bunch<NBUNCH; ++bunch){
          const double ly_raw =
            (adcint[bunch] - zero_pe * (double)(INT_RANGE)/(double)intrange_calib)
            * (m1_ref_base - m0_ref_base)
            / ((one_pe_base - zero_pe_base) * (m1_ref - m0_ref));
          out_lightyield[idx][bunch] = ly_raw * corr_factor;
        }

        // leading/trailing_fromadc
        (*out_leading_fromadc)[idx]  = calculate_leading_fromadc(wf);
        (*out_trailing_fromadc)[idx] = calculate_trailing_fromadc(wf);

        // copy existing leading/trailing if available
        if (leading && (int)leading->size() > ch)   (*out_leading)[idx]  = leading->at(ch);
        if (trailing && (int)trailing->size() > ch) (*out_trailing)[idx] = trailing->at(ch);
      }
    }
    out_spillnum = 0;
    for (int ib = 0; ib < 15; ++ib) {
      if (spill_bits[ib]) {
        out_spillnum |= (1 << ib);  // bit1(2^0)〜bit15(2^14)
      }
    }
    tout.Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "[OK] wrote: " << outfile << std::endl;
}

// -----------------------------
// Batch scan & process
// -----------------------------
static bool fileExists(const fs::path& p) {
  std::error_code ec; return fs::exists(p, ec);
}

struct CsvItem {
  fs::path path;         // .../calibration/calibresult/calib_<RUN>.csv
  std::string run;       // <RUN>
  std::time_t mtime = 0;
  std::uintmax_t size = 0; // file size in bytes
};

static std::vector<CsvItem> listCalibCsv(const fs::path& calibDir) {
  std::vector<CsvItem> v;
  if (!fileExists(calibDir)) return v;
  std::regex re(R"(calib_(.+)\.csv$)");
  for (const auto& de : fs::directory_iterator(calibDir)) {
    if (!de.is_regular_file()) continue;
    const auto& p = de.path();
    if (p.extension() != ".csv") continue;
    std::smatch m;
    auto name = p.filename().string();
    if (!std::regex_match(name, m, re)) continue;
    CsvItem it;
    it.path = p;
    it.run = m[1].str();
    std::error_code ec;
    auto ft = fs::last_write_time(p, ec);
    it.mtime = (!ec) ? filetime_to_time_t(ft) : 0;

    std::error_code ec2;
    auto fsize = fs::file_size(p, ec2);
    if (!ec2) {
      it.size = fsize;
    } else it.size = 0;
    v.push_back(it);
  }
  std::sort(v.begin(), v.end(), [](const CsvItem& a, const CsvItem& b){
    if (a.mtime != b.mtime) return a.mtime < b.mtime;
    return a.run < b.run;
  });
  return v;
}

static int scanAndConvert(const std::string& base,
                          const std::string& chmapFile,
                          const std::string& refgainCsv,
                          const std::string& spillchmapFile,
                          const std::string& lycorrCsv)
{
  const fs::path rootDir   = fs::path(base) / "rootfile";
  const fs::path outDir    = fs::path(base) / "rootfile_aftercalib";
  const fs::path calibDir  = fs::path(base) / "calibration" / "calibresult";
  const fs::path chmapPath = fs::path(base) / "chmap" / chmapFile;
  const fs::path refgainPath = fs::path(base) / "calibration" / "ReferenceGain" / refgainCsv;
  const fs::path lycorrPath  = fs::path(base) / "calibration" / "lightyield_correctionfactor" /lycorrCsv;
  const fs::path spillmapPath = fs::path(base) / "chmap" / spillchmapFile;

  // CSV must be >= 1KB and "stable" (no update) for at least 60 seconds
  constexpr std::uintmax_t MIN_CSV_SIZE_BYTES = 1024ull;
  constexpr std::time_t    CSV_STABLE_SEC     = 60;

  std::error_code ec;
  fs::create_directories(outDir, ec);

  auto csvs = listCalibCsv(calibDir);
  if (csvs.empty()) {
    std::cout << "[INFO] No calibration CSVs found in " << calibDir << std::endl;
    return 0;
  }

  // current time for stability judgement
  std::time_t now_t = std::time(nullptr);

  int converted = 0;
  for (const auto& ci : csvs) {
    // --- CSV size & stability check ---
    bool stable = false;
    if (ci.size >= MIN_CSV_SIZE_BYTES && ci.mtime != 0) {
      if (std::difftime(now_t, ci.mtime) >= CSV_STABLE_SEC) {
        stable = true;
      }
    }
    if (!stable) {
      std::cout << "[SKIP] CSV not ready (too small or recently updated): "
                << ci.path << " (size=" << ci.size << " bytes)\n";
      continue;
    }

    const std::string run = ci.run; // e.g., run00097_0_9999
    const fs::path inRoot  = rootDir   / (run + ".root");
    const fs::path outRoot = outDir    / (run + "_lightyield.root");
    const fs::path calibCsv= ci.path;

    if (!fileExists(inRoot)) {
      std::cout << "[SKIP] Input ROOT missing: " << inRoot << std::endl;
      continue;
    }
    if (fileExists(outRoot)) {
      // Already converted
      continue;
    }

    std::cout << "[RUN] " << run << std::endl;
    convertlightyield_rayraw_(inRoot.string().c_str(),
                              outRoot.string().c_str(),
                              chmapPath.string().c_str(),
                              refgainPath.string().c_str(),
                              calibCsv.string().c_str(),
                              lycorrPath.string().c_str(),
                              spillmapPath.string().c_str());
    ++converted;
    if (g_stop) break;
  }
  std::cout << "[DONE] Converted " << converted << " run(s)." << std::endl;
  return converted;
}

// -----------------------------
// CLI & main
// -----------------------------
static void printUsage(const char* prog) {
  std::cout <<
    "Usage: " << prog << " [options]\n"
    "Options:\n"
    "  -h, --help           Show this help and exit\n"
    "  --base <DIR>         Base directory (contains rootfile/, rootfile_aftercalib/, calibration/, chmap/)\n"
    "                       Default: /group/nu/ninja/work/otani/FROST_beamdata/test\n"
    "  --chmap <FILENAME>   Chmap filename under chmap/ (default: chmap_20251009.txt)\n"
    "  --refgain <FILE>     Reference gain CSV under calibration/ReferenceGain/ (default: refgain.csv)\n"
    "  --spillchmap <FILE>  Spillnum chmap filename under chmap/ (default: chmap_spillnum20251111.txt)\n"
    "  --lycorr <FILE>      Lightyield correction CSV under calibration/lightyield_correctionfactor (default: lightyield_correctionfactor.csv)\n"
    "  --watch <SEC>        Watch interval seconds (0 disables; default 60)\n"
    "  --oneshot <0|1>      If 1, run one scan and exit (default 0)\n";
}

int main(int argc, char** argv) {
  std::string base = FrostmonConfig::OUTPUT_DIR;
  std::string chmapFile = FrostmonConfig::CHMAP_FILE;
  std::string refgainCsv = FrostmonConfig::REFGAIN_CSV;
  std::string spillchmapFile = FrostmonConfig::SPILL_CHMAP_FILE;
  std::string lycorrCsv = FrostmonConfig::LIGHTYIELD_CORR_CSV;
  int watchSec = 60;    // default watch mode ON every 60s
  bool oneshot = false;

  std::signal(SIGINT,  handle_sigint);
  std::signal(SIGTERM, handle_sigint);

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") { printUsage(argv[0]); return 0; }
    auto next = [&](int& i)->std::string {
      if (i + 1 < argc) return std::string(argv[++i]);
      std::cerr << "Missing value after " << a << "\n";
      printUsage(argv[0]); std::exit(1);
    };
    if (a == "--base")        base       = next(i);
    else if (a == "--chmap")  chmapFile  = next(i);
    else if (a == "--refgain")refgainCsv = next(i);
    else if (a == "--spillchmap") spillchmapFile = next(i);
    else if (a == "--lycorr") lycorrCsv  = next(i);
    else if (a == "--watch")  watchSec   = std::stoi(next(i));
    else if (a == "--oneshot"){ std::string v = next(i); oneshot = (v=="1"||v=="true"||v=="yes"); }
    else { std::cerr << "Unknown option: " << a << "\n"; printUsage(argv[0]); return 1; }
  }

  std::cout << "[CFG] base=" << base
            << " chmap=" << chmapFile
            << " refgain=" << refgainCsv
            << " spillchmap=" << spillchmapFile
            << " lycorr=" << lycorrCsv
            << " watch=" << watchSec << "s"
            << " oneshot=" << (oneshot?1:0) << "\n";

  if (oneshot || watchSec <= 0) {
    scanAndConvert(base, chmapFile, refgainCsv, spillchmapFile, lycorrCsv);
    return 0;
  }

  std::cout << "[WATCH] Start watching every " << watchSec << " seconds. Press Ctrl-C to stop.\n";
  while (!g_stop) {
    std::time_t now = std::time(nullptr);
    std::cout << "\n[WATCH] Scan at " << std::put_time(std::localtime(&now), "%F %T") << std::endl;
    scanAndConvert(base, chmapFile, refgainCsv, spillchmapFile, lycorrCsv);

    for (int s = 0; s < watchSec && !g_stop; ++s) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
  }
  std::cout << "\n[WATCH] Stopped.\n";
  return 0;
}
