#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TDatime.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TError.h>

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <ctime>
#include <cctype>
#include <regex>
#include <memory>
#include "/home/nu/notani/FROST_monitor/config/config.hpp"

// ------------------------
// Global configuration
// ------------------------

static const std::string PATH_BSD_DIR =
  FrostmonConfig::OUTPUT_DIR + "/bsd/";
static const std::string PATH_LY_DIR  =
  FrostmonConfig::OUTPUT_DIR + "/rootfile_aftercalib/";

static const std::string PATH_OUT_BASE =
  FrostmonConfig::OUTPUT_DIR + "/dataquality_withBSD/";

static const std::string PATH_POT_PLOT =
  PATH_OUT_BASE + "pot_plot/";

static const std::string PATH_POT_INFO =
  PATH_OUT_BASE + "pot_info/";

static const std::string PATH_POT_10MIN_CSV =
  PATH_POT_INFO + "pot_spill_10min.csv";
static const std::string PATH_MISSED_SPILLS_CSV =
  PATH_POT_INFO + "missed_spills.csv";

static const std::string PATH_POT_POINTS =
  PATH_POT_PLOT+ "pot_points.tsv";
static const std::string PATH_POT_PROC   =
  PATH_POT_PLOT + "processed_lightyield.tsv";
static const std::string PATH_OUT_PDF    =
  PATH_POT_PLOT + "accumulated_pot_withBSD.pdf";

static const std::string PATH_EVENTRATE_DIR =
  PATH_OUT_BASE + "eventrate_plot/";
static const std::string PATH_EVENT  =
  PATH_EVENTRATE_DIR + "eventrate.tsv";
static const std::string PATH_EVENTRATE_PDF =
  PATH_EVENTRATE_DIR + "eventrate_plot.pdf";

// For "max 1 event per spill" counting
static const std::string PATH_EVENT_SPILL =
  PATH_EVENTRATE_DIR + "eventrate_spill.tsv";
static const std::string PATH_EVENTRATE_PDF_SPILL =
  PATH_EVENTRATE_DIR + "eventrate_plot_spill.pdf";

static const std::string PATH_ACQUIRED_BUNCH_DIR =
  PATH_OUT_BASE + "acquired_bunch/";

static const char* BSD_TREE_NAME  = "bsd";
static const char* LY_TREE_NAME   = "tree";

// file readiness
static const Long_t MIN_LY_SIZE_BYTES = 10L * 1024; // 10 KB
static const Long_t STABLE_SEC        = 60;               // consider stable if mtime older than 60 s

// matching
static const int SPILL_MOD = FrostmonConfig::SPILL_MOD; // 15-bit modulo for spillnum matching

static const int MAX_TIME_DIFF = FrostmonConfig::MAX_TIME_DIFF; // seconds

static const Double_t LIGHTMAX_MIN = FrostmonConfig::LIGHTMAX_MIN; //threshold for event

// ------------------------
// BSD time cut (hardcoded Unix time)
// ------------------------

// Keep only spills with trg_sec >= 2026/01/10 00:00:00 JST.
// 2026/01/10 00:00:00 JST = 2026/01/09 15:00:00 UTC = 1767970800
static const long long BSD_CUT_START_SEC = 1767970800LL;

// ------------------------
// Event-rate plot time cut (hardcoded Unix time)
// ------------------------
// Keep only points with time >= 2026/01/17 00:00:00 JST.
// 2026/01/17 00:00:00 JST = 2026/01/16 15:00:00 UTC = 1768575600
static const long long EVENT_CUT_START_SEC = 1768575600LL;


// ------------------------
// Small helpers
// ------------------------

struct FileInfoLite {
  std::string path;
  Long_t size  = 0;
  Long_t mtime = 0;
};

static bool GetPathInfo(const std::string& p,
                        Long_t& id, Long_t& size,
                        Long_t& flags, Long_t& mtime)
{
  return gSystem->GetPathInfo(p.c_str(), &id, &size, &flags, &mtime) == 0;
}

static bool IsLikelyStableFile(const std::string& path, int wait_ms = 500)
{
  Long_t id1=0, sz1=0, fl1=0, mt1=0;
  if (!GetPathInfo(path, id1, sz1, fl1, mt1)) return false;
  gSystem->Sleep(wait_ms);
  Long_t id2=0, sz2=0, fl2=0, mt2=0;
  if (!GetPathInfo(path, id2, sz2, fl2, mt2)) return false;
  return (sz1 == sz2);
}

static bool IsReadyLightyieldFile(const FileInfoLite& fi)
{
  if (fi.size < MIN_LY_SIZE_BYTES) return false;
  if (!IsLikelyStableFile(fi.path, 500)) return false;

  const Long_t nowSec = static_cast<Long_t>(time(nullptr));
  if (nowSec - fi.mtime < STABLE_SEC) return false;

  return true;
}

static std::vector<FileInfoLite> ListLightyieldFiles(const std::string& dir)
{
  std::vector<FileInfoLite> out;
  TSystemDirectory d("lydir", dir.c_str());
  TList* files = d.GetListOfFiles();
  if (!files) return out;

  TIter it(files);
  while (TSystemFile* f = (TSystemFile*)it()) {
    TString name = f->GetName();
    if (f->IsDirectory()) continue;
    if (!name.EndsWith("_lightyield.root")) continue;

    std::string path = dir + name.Data();
    Long_t id=0, sz=0, fl=0, mt=0;
    if (GetPathInfo(path, id, sz, fl, mt)) {
      out.push_back({path, sz, mt});
    }
  }

  std::sort(out.begin(), out.end(),
            [](const FileInfoLite& a, const FileInfoLite& b){
              return a.mtime < b.mtime; // older first
            });
  return out;
}

static std::vector<FileInfoLite> ListBsdFiles(const std::string& dir)
{
  std::vector<FileInfoLite> out;
  TSystemDirectory d("bsddir", dir.c_str());
  TList* files = d.GetListOfFiles();
  if (!files) return out;

  TIter it(files);
  while (TSystemFile* f = (TSystemFile*)it()) {
    TString name = f->GetName();
    if (f->IsDirectory()) continue;
    if (!name.EndsWith(".root")) continue;
    if (!name.BeginsWith("bsd_run")) continue;

    std::string path = dir + name.Data();
    Long_t id=0, sz=0, fl=0, mt=0;
    if (GetPathInfo(path, id, sz, fl, mt)) {
      out.push_back({path, sz, mt});
    }
  }

  std::sort(out.begin(), out.end(),
            [](const FileInfoLite& a, const FileInfoLite& b){
              return a.mtime < b.mtime;
            });
  return out;
}

// ------------------------
// Utilities for date formatting, POT formatting, and JS output
// ------------------------

// Convert Unix time → "YYYY/MM/DD"
static std::string FormatDateYMD(std::time_t t)
{
  std::tm lt{};
#if defined(_WIN32)
  localtime_s(&lt, &t);
#else
  lt = *std::localtime(&t);
#endif

  std::ostringstream oss;
  oss << (lt.tm_year + 1900) << "/";
  oss << std::setw(2) << std::setfill('0') << (lt.tm_mon + 1) << "/";
  oss << std::setw(2) << std::setfill('0') << lt.tm_mday;
  return oss.str();
}

// Convert Unix time → "YYYY/MM/DD/H:MM" (hour is not zero-padded)
static std::string FormatDateTimeYMDHM(std::time_t t)
{
  std::tm lt{};
#if defined(_WIN32)
  localtime_s(&lt, &t);
#else
  lt = *std::localtime(&t);
#endif

  std::ostringstream oss;
  oss << (lt.tm_year + 1900) << "/";
  oss << std::setw(2) << std::setfill('0') << (lt.tm_mon + 1) << "/";
  oss << std::setw(2) << std::setfill('0') << lt.tm_mday << "/";
  oss << (lt.tm_hour) << ":";
  oss << std::setw(2) << std::setfill('0') << lt.tm_min;
  return oss.str();
}

// Format POT as scientific notation for HTML:
// Example: 4.41 × 10^19 becomes "4.41 &times; 10<sup>19</sup>"
static std::string FormatSciForHtml(double val)
{
  if (!std::isfinite(val) || val <= 0.0) return "0";

  int exp = static_cast<int>(std::floor(std::log10(val)));
  double mant = val / std::pow(10.0, exp);

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(4) << mant
      << " &times; 10<sup>" << exp << "</sup>";
  return oss.str();
}

// Write JS file containing: document.write(" ... ");
static void WriteJsFile(const std::string& path, const std::string& payload)
{
  std::ofstream ofs(path, std::ios::trunc);
  if (!ofs) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open JS output file: %s", path.c_str());
    return;
  }

  ofs << "document.write(\"" << payload << "\");\n";
}

// ------------------------
// BSD spills
// ------------------------

// Rule for acquired bunches per DAQ run, read from:
//   PATH_ACQUIRED_BUNCH_DIR + "acquired_bunch_rules.txt"
// Format example:
//   # run_max  acquired_bunch
//   5    1,2
//   9999 0
// Here, runs 1-5 have only bunch 1 and 2 acquired, and runs 6-9999 have all 8 bunches.
struct AcquiredBunchRule {
  int run_max = 0;
  bool all_bunches = false;          // true if all 8 bunches are acquired
  std::vector<int> bunch_list;       // 1-based bunch indices (1..8) when not all_bunches
};

static void LoadAcquiredBunchRules(std::vector<AcquiredBunchRule>& rules)
{
  rules.clear();

  const std::string path = PATH_ACQUIRED_BUNCH_DIR + FrostmonConfig::ACQUIRED_BUNCH_RULES_FILE;
  std::ifstream fin(path);
  if (!fin) {
    ::Info("dataqualityplot_withBSD",
           "No acquired_bunch_rules.txt found at %s. Assuming all 8 bunches for all runs.",
           path.c_str());
    return;
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    int run_max = 0;
    std::string bunch_str;
    if (!(iss >> run_max >> bunch_str)) continue;

    AcquiredBunchRule rule;
    rule.run_max = run_max;
    if (bunch_str == "0") {
      // 0 means "all 8 bunches"
      rule.all_bunches = true;
    } else {
      rule.all_bunches = false;
      std::stringstream ss(bunch_str);
      std::string token;
      while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        int b = std::stoi(token);
        if (b >= 1 && b <= 8) {
          rule.bunch_list.push_back(b);
        }
      }
    }
    rules.push_back(std::move(rule));
  }

  std::sort(rules.begin(), rules.end(),
            [](const AcquiredBunchRule& a, const AcquiredBunchRule& b){
              return a.run_max < b.run_max;
            });
}

// Return bit mask for acquired bunches for a given run number.
// bit 0 corresponds to bunch 1, ..., bit 7 to bunch 8.
// If no rule matches, or run <= 0, return 0xFF (all 8 bunches).
static unsigned char GetAcquiredBunchMaskForRun(int run,
                                                const std::vector<AcquiredBunchRule>& rules)
{
  if (run <= 0 || rules.empty()) return 0xFF;

  for (const auto& r : rules) {
    if (run <= r.run_max) {
      if (r.all_bunches || r.bunch_list.empty()) return 0xFF;
      unsigned char mask = 0u;
      for (int b : r.bunch_list) {
        if (b >= 1 && b <= 8) {
          mask |= static_cast<unsigned char>(1u << (b - 1));
        }
      }
      if (mask == 0u) mask = 0xFF;
      return mask;
    }
  }
  return 0xFF;
}


struct BsdSpill {
  int spillnum_full = 0;   // full spillnumber
  int spillnum_mod  = 0;   // low 15 bits
  int good_spill_flag = 0;
  int spill_flag = 0;
  int trg_sec = 0;         // trg_sec[2]
  // Total POT over all 8 bunches (ct_np[4][0])
  double pot = 0.0;
  // POT per bunch (ct_np[4][1]..ct_np[4][8]):
  // pot_bunch[0] -> bunch 1, ..., pot_bunch[7] -> bunch 8
  double pot_bunch[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
};

static void LoadAllBsdSpills(std::vector<BsdSpill>& bsd_spills,
                             std::unordered_map<int, std::vector<size_t>>& index_by_mod)
{
  bsd_spills.clear();
  index_by_mod.clear();

  auto files = ListBsdFiles(PATH_BSD_DIR);
  if (files.empty()) {
    ::Info("dataqualityplot_withBSD", "No BSD root files in %s", PATH_BSD_DIR.c_str());
    return;
  }

  for (const auto& fi : files) {
    TFile f(fi.path.c_str(), "READ");
    if (f.IsZombie()) {
      ::Warning("dataqualityplot_withBSD", "Failed to open BSD file: %s", fi.path.c_str());
      continue;
    }

    TTree* t = (TTree*)f.Get(BSD_TREE_NAME);
    if (!t) {
      ::Warning("dataqualityplot_withBSD", "No tree '%s' in BSD file: %s",
                BSD_TREE_NAME, fi.path.c_str());
      f.Close();
      continue;
    }

    t->SetBranchStatus("*", 0);

    Int_t spillnum = 0;
    Int_t trg_sec_arr[3] = {0,0,0};
    Int_t good_spill_flag = 0;
    Int_t spill_flag = 0;
    Double_t ct_np[5][9];

    t->SetBranchStatus("spillnum", 1);
    t->SetBranchAddress("spillnum", &spillnum);

    t->SetBranchStatus("trg_sec", 1);
    t->SetBranchAddress("trg_sec", trg_sec_arr);

    t->SetBranchStatus("good_spill_flag", 1);
    t->SetBranchAddress("good_spill_flag", &good_spill_flag);

    t->SetBranchStatus("spill_flag", 1);
    t->SetBranchAddress("spill_flag", &spill_flag);

    t->SetBranchStatus("ct_np", 1);
    t->SetBranchAddress("ct_np", ct_np);

    const Long64_t nent = t->GetEntries();
    for (Long64_t ie = 0; ie < nent; ++ie) {
      t->GetEntry(ie);

      // If you want to include bad spills also, comment out next 2 lines
      // if (good_spill_flag == 0) continue;
      if (spill_flag == 0) continue;

      BsdSpill s;
      s.spillnum_full   = spillnum;
      s.spillnum_mod    = spillnum & (SPILL_MOD - 1); // 15 bits
      s.good_spill_flag = good_spill_flag;
      s.spill_flag      = spill_flag;
      s.trg_sec         = trg_sec_arr[2];             // trg_sec[2]

      // Skip spills earlier than the configured start time.
      if (static_cast<long long>(s.trg_sec) < BSD_CUT_START_SEC) {
        continue;
      }

      s.pot             = ct_np[4][0];                // total POT (all bunches)

      // Store POT per bunch: ct_np[4][1]..ct_np[4][8]
      for (int b = 0; b < 8; ++b) {
        s.pot_bunch[b] = ct_np[4][b + 1];
      }

      size_t idx = bsd_spills.size();
      bsd_spills.push_back(s);
      index_by_mod[s.spillnum_mod].push_back(idx);
    }

    f.Close();
  }

  ::Info("dataqualityplot_withBSD",
         "Loaded %zu BSD spills from %zu files.",
         bsd_spills.size(), files.size());
}

// ------------------------
// Processed lightyield files (cache)
// ------------------------

static void LoadProcessedLyFiles(std::unordered_map<std::string, Long_t>& seen)
{
  seen.clear();
  std::ifstream fin(PATH_POT_PROC);
  if (!fin) return;

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string path;
    Long_t mt = 0;
    if (!(iss >> std::quoted(path) >> mt)) continue;
    seen[path] = mt;
  }
}

static void SaveProcessedLyFiles(const std::unordered_map<std::string, Long_t>& seen)
{
  std::ofstream fout(PATH_POT_PROC, std::ios::trunc);
  fout << "#path\tmtime\n";
  for (const auto& kv : seen) {
    fout << std::quoted(kv.first) << "\t" << kv.second << "\n";
  }
}

// ------------------------
// POT points file
// ------------------------

static void AppendPotPoints(const std::vector<std::pair<double,double>>& new_points)
{
  // This function rewrites the entire POT points file
  // using the provided list. This ensures that all available lightyield files
  // are always fully reflected in the cache.

  gSystem->mkdir(PATH_POT_PLOT.c_str(), true);
  std::ofstream fout(PATH_POT_POINTS, std::ios::trunc);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open pot_points file for append: %s",
              PATH_POT_POINTS.c_str());
    return;
  }

  // Always write header
  fout << "#time pot_increment\n";

  fout << std::setprecision(16);
  for (const auto& p : new_points) {
    fout << p.first << " " << p.second << "\n";
  }
}

// Append event timestamps (Unix time [s]) to PATH_EVENT.
// Each line:  time
static void AppendNuEvents(const std::vector<double>& nu_times)
{
  // This function rewrites the entire events
  // file using the provided list. This ensures that all available lightyield
  // files are always fully reflected in the cache.
  gSystem->mkdir(PATH_EVENTRATE_DIR.c_str(), true);
  std::ofstream fout(PATH_EVENT, std::ios::trunc);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open events file for append: %s",
              PATH_EVENT.c_str());
    return;
  }

  // Always write header
  fout << "#time\n";

  fout << std::setprecision(16);
  for (double t : nu_times) {
    fout << t << "\n";
  }
}

// Append event timestamps (Unix time [s]) to PATH_EVENT_SPILL.
// Each line: time
// Here, at most one entry is stored per spill (if any bunch had a hit).
static void AppendNuEventsSpill(const std::vector<double>& nu_times_spill)
{
  gSystem->mkdir(PATH_EVENTRATE_DIR.c_str(), true);
  std::ofstream fout(PATH_EVENT_SPILL, std::ios::trunc);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open eventsrate (per spill) file for write: %s",
              PATH_EVENT_SPILL.c_str());
    return;
  }

  // Always write header
  fout << "#time\n";

  fout << std::setprecision(16);
  for (double t : nu_times_spill) {
    fout << t << "\n";
  }
}

// ------------------------
// Process lightyield file -> find POT increments
// ------------------------

// Extract DAQ run number from a lightyield file path.
// Example: ".../run00003_140000_149999_lightyield.root" -> 3
static int ExtractRunNumberFromLyPath(const std::string& path)
{
  std::regex re("run(\\d{5})_");
  std::smatch m;
  if (std::regex_search(path, m, re)) {
    try {
      return std::stoi(m[1].str());
    } catch (...) {
      return -1;
    }
  }
  return -1;
}


static void ProcessLightyieldFileForPOT(
    const std::string& path,
    const std::vector<BsdSpill>& bsd_spills,
    const std::unordered_map<int, std::vector<size_t>>& index_by_mod,
    std::vector<std::pair<double,double>>& out_new_points,
    std::vector<double>& out_new_nu_times,
    std::vector<double>& out_new_nu_times_spill,
    const std::vector<AcquiredBunchRule>& acq_rules,
    std::unordered_set<size_t>* out_matched_bsd_indices)
{
  TFile f(path.c_str(), "READ");
  if (f.IsZombie()) {
    ::Warning("dataqualityplot_withBSD", "Failed to open lightyield file: %s", path.c_str());
    return;
  }

  TTree* t = (TTree*)f.Get(LY_TREE_NAME);
  if (!t) {
    ::Warning("dataqualityplot_withBSD", "No tree '%s' in lightyield file: %s",
              LY_TREE_NAME, path.c_str());
    f.Close();
    return;
  }

  t->SetBranchStatus("*", 0);

  Int_t spillnum = 0;
  Double_t unixtime_arr[272];
  Double_t lightyield_arr[272][8];
  Int_t cablenum_arr[272];  // channel mapping (Y: 1-140, X: 201-332)

  t->SetBranchStatus("spillnum",  1);
  t->SetBranchAddress("spillnum", &spillnum);

  t->SetBranchStatus("unixtime",  1);
  t->SetBranchAddress("unixtime", unixtime_arr);

  t->SetBranchStatus("lightyield", 1);
  t->SetBranchAddress("lightyield", lightyield_arr);

  // cablenum gives the mapping of each channel:
  //  - 1..140  : Y-plane
  //  - 201..332: X-plane
  t->SetBranchStatus("cablenum", 1);
  t->SetBranchAddress("cablenum", cablenum_arr);

  const Long64_t nent = t->GetEntries();

  // Determine DAQ run number and acquired bunch mask for this lightyield file.
  // The run number is taken from the lightyield file name
  // (e.g. run00003_... -> run = 3).
  const int run_number = ExtractRunNumberFromLyPath(path);
  const unsigned char acquired_bunch_mask = GetAcquiredBunchMaskForRun(run_number, acq_rules);

  std::unordered_set<size_t> used_bsd_indices;
  used_bsd_indices.reserve(1024);

  Long64_t match_count   = 0;
  Long64_t event_checked = 0;

  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    ++event_checked;

    const double ut_d = unixtime_arr[0];
    if (!std::isfinite(ut_d)) continue;

    int ut = static_cast<int>(std::llround(ut_d)); // sec
    const int key = spillnum & (SPILL_MOD - 1);

    auto it = index_by_mod.find(key);
    if (it == index_by_mod.end()) continue;

    int best_dt = MAX_TIME_DIFF + 1;
    size_t best_idx = std::numeric_limits<size_t>::max();

    for (size_t idx : it->second) {
      const int dt = std::abs(bsd_spills[idx].trg_sec - ut);
      if (dt <= MAX_TIME_DIFF && dt < best_dt) {
        best_dt = dt;
        best_idx = idx;
      }
    }
    if (best_idx == std::numeric_limits<size_t>::max()) continue;

    const BsdSpill& s = bsd_spills[best_idx];

    // Count "events" per bunch:
    // For each bunch, compute the maximum lightyield in the Y-plane and X-plane
    // separately, using cablenum mapping:
    //   Y-plane: cablenum in [1, 140]
    //   X-plane: cablenum in [201, 332]
    // A bunch is counted as one event if both maxima are >= LIGHTMAX_MIN.
    // When only a subset of bunches is acquired in this run, we consider
    // only those bunches.
    // In addition, we also keep track of whether this spill had at least one
    // hit bunch, to build a "max 1 event per spill" event-rate plot.
    int  n_nu_events_this_spill   = 0;     // per-bunch counting
    bool has_nu_event_this_spill  = false; // per-spill flag

    if (acquired_bunch_mask == 0xFF) {
      // All 8 bunches acquired: use all bunch indices (0..7)
      for (int b = 0; b < 8; ++b) {
        double ly_max_y = -1.0; // max on Y-plane
        double ly_max_x = -1.0; // max on X-plane

        for (int ch = 0; ch < 272; ++ch) {
          const int cable = cablenum_arr[ch];
          const double v  = lightyield_arr[ch][b];
          if (!std::isfinite(v)) continue;

          // Y-plane: cablenum 1..140
          if (cable >= 1 && cable <= 140) {
            if (v > ly_max_y) ly_max_y = v;
          }
          // X-plane: cablenum 201..332
          else if (cable >= 201 && cable <= 332) {
            if (v > ly_max_x) ly_max_x = v;
          }
        }

        // Require both X and Y maxima above threshold
        if (ly_max_x >= LIGHTMAX_MIN && ly_max_y >= LIGHTMAX_MIN) {
          ++n_nu_events_this_spill;
          has_nu_event_this_spill = true;
        }
      }
    } else {
      // Only specific bunches acquired: consider only those bunch indices.
      for (int b = 0; b < 8; ++b) {
        if (!(acquired_bunch_mask & (1u << b))) continue;
        double ly_max_y = -1.0; // max on Y-plane
        double ly_max_x = -1.0; // max on X-plane

        for (int ch = 0; ch < 272; ++ch) {
          const int cable = cablenum_arr[ch];
          const double v  = lightyield_arr[ch][b];
          if (!std::isfinite(v)) continue;

          if (cable >= 1 && cable <= 140) {
            if (v > ly_max_y) ly_max_y = v;
          } else if (cable >= 201 && cable <= 332) {
            if (v > ly_max_x) ly_max_x = v;
          }
        }

        if (ly_max_x >= LIGHTMAX_MIN && ly_max_y >= LIGHTMAX_MIN) {
          ++n_nu_events_this_spill;
          has_nu_event_this_spill = true;
        }
      }
    }

    // Each event (per bunch) is stored as one entry
    // at this spill time.
    for (int i = 0; i < n_nu_events_this_spill; ++i) {
      out_new_nu_times.push_back(static_cast<double>(s.trg_sec));
    }

    // For the "per-spill" counting, we store at most one entry per spill
    // (if any bunch had a hit).
    if (has_nu_event_this_spill) {
      out_new_nu_times_spill.push_back(static_cast<double>(s.trg_sec));
    }

    // For POT accumulation, avoid double counting the same BSD spill within this file.
    // The POT used here depends on the acquired bunch configuration of this run.
    double pot_increment = 0.0;
    if (acquired_bunch_mask == 0xFF) {
      // All bunches: use total POT
      pot_increment = s.pot;
    } else {
      // Only specific bunches: sum POT over those bunches
      for (int b = 0; b < 8; ++b) {
        if (acquired_bunch_mask & (1u << b)) {
          pot_increment += s.pot_bunch[b];
        }
      }
    }

    if (used_bsd_indices.insert(best_idx).second) {
      out_new_points.emplace_back(static_cast<double>(s.trg_sec), pot_increment);
      ++match_count;

      if (out_matched_bsd_indices) {
        out_matched_bsd_indices->insert(best_idx);
      }
    }
  }

  f.Close();

  ::Info("dataqualityplot_withBSD",
         "Processed LY file: %s (events=%lld, matched_spills=%lld)",
         path.c_str(), event_checked, match_count);
}

// ------------------------
// Build and save Accumulated POT plot (Delivered + Recorded)
// ------------------------
static void BuildAndSaveAccumulatedPotPlot(
    const std::vector<BsdSpill>& bsd_spills)
{
  // 1) Recorded by FROST: read pot_points.tsv
  std::vector<std::pair<double,double>> points; // (time, pot_increment)

  {
    std::ifstream fin(PATH_POT_POINTS);
    if (fin) {
      std::string line;
      while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double t = 0.0, p = 0.0;
        if (!(iss >> t >> p)) continue;
        if (!std::isfinite(t) || !std::isfinite(p)) continue;
        points.emplace_back(t, p);
      }
      fin.close();
    }
  }

  // Recorded-by-FROST curve (red line)
  std::vector<double> vx_rec, vy_rec;
  if (!points.empty()) {
    std::sort(points.begin(), points.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    const int n = (int)points.size();
    vx_rec.resize(n);
    vy_rec.resize(n);

    double accum = 0.0;
    for (int i = 0; i < n; ++i) {
      accum += points[i].second;
      vx_rec[i] = points[i].first;
      vy_rec[i] = accum;
    }
  }

  // 2) Delivered POT from all BSD spills (black line)
  std::vector<BsdSpill> bsd_sorted = bsd_spills;
  std::sort(bsd_sorted.begin(), bsd_sorted.end(),
            [](const BsdSpill& a, const BsdSpill& b) {
              return a.trg_sec < b.trg_sec;
            });

  std::vector<double> vx_del, vy_del;
  if (!bsd_sorted.empty()) {
    const int n = (int)bsd_sorted.size();
    vx_del.resize(n);
    vy_del.resize(n);

    double accum = 0.0;
    for (int i = 0; i < n; ++i) {
      accum += bsd_sorted[i].pot;
      vx_del[i] = (double)bsd_sorted[i].trg_sec;
      vy_del[i] = accum;
    }
  }

  // Scalar summary for JS outputs
  const int delivered_spills = static_cast<int>(bsd_sorted.size());
  const int recorded_spills  = static_cast<int>(vx_rec.size());

  const double delivered_pot_total =
      vx_del.empty() ? 0.0 : vy_del.back();
  const double recorded_pot_total  =
      vx_rec.empty() ? 0.0 : vy_rec.back();

  // Create output directories
  gSystem->mkdir(PATH_OUT_BASE.c_str(), true);
  gSystem->mkdir(PATH_POT_INFO.c_str(), true);

  // period.js: "YYYY/MM/DD~YYYY/MM/DD"
  if (!bsd_sorted.empty()) {
    std::time_t t0 = static_cast<std::time_t>(bsd_sorted.front().trg_sec);
    std::time_t t1 = static_cast<std::time_t>(bsd_sorted.back().trg_sec);
    std::string period_str = FormatDateYMD(t0) + "~" + FormatDateYMD(t1);
    WriteJsFile(PATH_POT_INFO + "period.js", period_str);
  }

  // delivered_spills.js / recorded_spills.js
  WriteJsFile(PATH_POT_INFO + "delivered_spills.js",
              std::to_string(delivered_spills));
  WriteJsFile(PATH_POT_INFO + "recorded_spills.js",
              std::to_string(recorded_spills));

  // delivered_pot.js / recorded_pot.js / efficiency.js
  if (delivered_pot_total > 0.0) {
    WriteJsFile(PATH_POT_INFO + "delivered_pot.js",
                FormatSciForHtml(delivered_pot_total));
    WriteJsFile(PATH_POT_INFO + "recorded_pot.js",
                FormatSciForHtml(recorded_pot_total));

    double eff = (recorded_pot_total / delivered_pot_total) * 100.0;
    std::ostringstream oe;
    oe << std::fixed << std::setprecision(4) << eff << " %";
    WriteJsFile(PATH_POT_INFO + "efficiency.js", oe.str());
  } else {
    WriteJsFile(PATH_POT_INFO + "delivered_pot.js", "0");
    WriteJsFile(PATH_POT_INFO + "recorded_pot.js", "0");
    WriteJsFile(PATH_POT_INFO + "efficiency.js",   "0 %");
  }

  if (vx_del.empty() && vx_rec.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No Delivered nor Recorded POT points to plot yet.");
    return;
  }

  // 3) Draw Accumulated POT plot
  gStyle->SetOptStat(0);

  TCanvas c("c_pot", "Accumulated POT (Delivered vs FROST)", 1600, 800);
  c.SetGrid();
  c.SetLeftMargin(0.12);

  bool hasDelivered = !vx_del.empty();
  bool hasRecorded  = !vx_rec.empty();

  TGraph* gr_del = nullptr;
  TGraph* gr_rec = nullptr;

  if (hasDelivered) {
    gr_del = new TGraph((int)vx_del.size(), vx_del.data(), vy_del.data());
    gr_del->SetTitle("Accumulated POT;time;Accumulated POT");
    gr_del->SetLineColor(kBlack);
    gr_del->SetLineWidth(2);
    gr_del->SetMarkerStyle(0);
    gr_del->SetMinimum(0);

    TAxis* ax = gr_del->GetXaxis();
    ax->SetTimeDisplay(1);
    ax->SetTimeOffset(0, "local");
    ax->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    ax->SetLabelSize(0.03);
    ax->SetLabelOffset(0.02);
    ax->SetTitleOffset(1.5);
    ax->SetNdivisions(520);

    gr_del->Draw("AL");
  }

  if (hasRecorded) {
    gr_rec = new TGraph((int)vx_rec.size(), vx_rec.data(), vy_rec.data());
    if (!hasDelivered) {
      gr_rec->SetTitle("Accumulated POT;time;Accumulated POT");
      TAxis* ax = gr_rec->GetXaxis();
      ax->SetTimeDisplay(1);
      ax->SetTimeOffset(0, "local");
      ax->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      ax->SetLabelSize(0.03);
      ax->SetLabelOffset(0.02);
      ax->SetTitleOffset(1.5);
      ax->SetNdivisions(520);
      gr_rec->Draw("AL");
    } else {
      gr_rec->Draw("L SAME");
    }
    gr_rec->SetLineColor(kRed);
    gr_rec->SetLineWidth(2);
    gr_rec->SetMarkerStyle(0);
    gr_rec->SetMinimum(0);
  }

  // Legend
  TLegend leg(0.15, 0.75, 0.45, 0.85);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  if (hasDelivered) leg.AddEntry(gr_del, "Delivered", "l");
  if (hasRecorded)  leg.AddEntry(gr_rec, "Recorded by FROST", "l");
  leg.Draw();

  c.SaveAs(PATH_OUT_PDF.c_str());

  ::Info("dataqualityplot_withBSD",
         "Saved Accumulated POT plot (Delivered + Recorded) to %s",
         PATH_OUT_PDF.c_str());

  delete gr_del;
  delete gr_rec;
}

static void BuildAndSavePotSpillCsv10Min(const std::vector<BsdSpill>& bsd_spills)
{
  // Recorded points (time, pot_increment) from pot_points.tsv
  std::vector<std::pair<int,double>> rec_points; // (trg_sec, pot_inc)
  {
    std::ifstream fin(PATH_POT_POINTS);
    if (fin) {
      std::string line;
      while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double t=0.0, p=0.0;
        if (!(iss >> t >> p)) continue;
        if (!std::isfinite(t) || !std::isfinite(p)) continue;
        rec_points.emplace_back((int)std::llround(t), p);
      }
    }
  }
  std::sort(rec_points.begin(), rec_points.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });

  // Delivered spills from BSD (time, pot_inc=pot)
  std::vector<BsdSpill> bsd_sorted = bsd_spills;
  std::sort(bsd_sorted.begin(), bsd_sorted.end(),
            [](const BsdSpill& a, const BsdSpill& b){ return a.trg_sec < b.trg_sec; });

  if (bsd_sorted.empty() && rec_points.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No delivered nor recorded data: skip 10-min CSV.");
    return;
  }

  // Determine time range
  const int BIN_SEC = 600; // 10 min
  int t_min = std::numeric_limits<int>::max();
  int t_max = 0;

  if (!bsd_sorted.empty()) {
    t_min = std::min(t_min, bsd_sorted.front().trg_sec);
    t_max = std::max(t_max, bsd_sorted.back().trg_sec);
  }
  if (!rec_points.empty()) {
    t_min = std::min(t_min, rec_points.front().first);
    t_max = std::max(t_max, rec_points.back().first);
  }

  // Align to 10-min boundaries
  int t0 = (t_min / BIN_SEC) * BIN_SEC;                 // floor
  int t1 = ((t_max + BIN_SEC - 1) / BIN_SEC) * BIN_SEC; // ceil

  gSystem->mkdir(PATH_POT_INFO.c_str(), true);

  std::ofstream fout(PATH_POT_10MIN_CSV, std::ios::trunc);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open 10-min CSV: %s", PATH_POT_10MIN_CSV.c_str());
    return;
  }

  fout << "#time, delivered spill, recorded spill, delivered POT, recorded POT, efficiency[%]\n";
  fout << std::setprecision(16);

  // Sweep with two pointers (linear, fast)
  size_t i_del = 0, i_rec = 0;
  long long delivered_spill = 0;
  long long recorded_spill  = 0;
  double delivered_pot = 0.0;
  double recorded_pot  = 0.0;

  for (int t = t0; t <= t1; t += BIN_SEC) {
    while (i_del < bsd_sorted.size() && bsd_sorted[i_del].trg_sec <= t) {
      delivered_spill++;
      delivered_pot += bsd_sorted[i_del].pot;
      ++i_del;
    }
    while (i_rec < rec_points.size() && rec_points[i_rec].first <= t) {
      recorded_spill++;
      recorded_pot += rec_points[i_rec].second;
      ++i_rec;
    }

    fout << FormatDateTimeYMDHM((std::time_t)t) << ","
         << delivered_spill << ","
         << recorded_spill << ","
         << delivered_pot << ","
         << recorded_pot << ","
         << 100.*(recorded_pot/delivered_pot) << "\n";
  }

  ::Info("dataqualityplot_withBSD",
         "Saved 10-min CSV to %s", PATH_POT_10MIN_CSV.c_str());
}

static void BuildAndSaveMissedSpillsCsv(
    const std::vector<BsdSpill>& bsd_spills,
    const std::unordered_set<size_t>& matched_bsd_indices)
{
  gSystem->mkdir(PATH_POT_INFO.c_str(), true);

  std::ofstream fout(PATH_MISSED_SPILLS_CSV, std::ios::trunc);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open missed spills CSV: %s", PATH_MISSED_SPILLS_CSV.c_str());
    return;
  }

  fout << "#time, unixtime, spillnum, POT\n";
  fout << std::setprecision(16);

  // 1) Collect indices of missed spills (i.e., not matched by any lightyield event)
  std::vector<size_t> missed;
  missed.reserve(bsd_spills.size());
  for (size_t i = 0; i < bsd_spills.size(); ++i) {
    if (matched_bsd_indices.find(i) != matched_bsd_indices.end()) continue;
    missed.push_back(i);
  }

  // 2) Sort by unixtime (trg_sec). If the time is identical, sort by spill number for stability.
  std::sort(missed.begin(), missed.end(),
            [&](size_t a, size_t b){
              const auto& sa = bsd_spills[a];
              const auto& sb = bsd_spills[b];
              if (sa.trg_sec != sb.trg_sec) return sa.trg_sec < sb.trg_sec;
              return sa.spillnum_full < sb.spillnum_full;
            });

  // 3) Write out missed spills in chronological order
  for (size_t idx : missed) {
    const auto& s = bsd_spills[idx];
    fout << FormatDateTimeYMDHM((std::time_t)s.trg_sec) << ", "
         << s.trg_sec << ", "
         << s.spillnum_full << ", "
         << s.pot << "\n";
  }

  ::Info("dataqualityplot_withBSD",
         "Saved missed spills CSV to %s", PATH_MISSED_SPILLS_CSV.c_str());
}

// ------------------------
// Build and save daily event-rate plot as a histogram:
//   X: date (one bin per day, labeled as "MM/DD")
//   Y: Number of events / 10^15 POT (for that day)
//   Drawn with error bars only (no filled bins). Also performs a constant fit.
// ------------------------
static void BuildAndSaveEventRatePlot()
{
  // 1) Load daily POT from pot_points.tsv
  std::ifstream fin_pot(PATH_POT_POINTS);
  if (!fin_pot) {
    ::Warning("dataqualityplot_withBSD",
              "No pot_points file found: %s (cannot build event-rate plot).",
              PATH_POT_POINTS.c_str());
    return;
  }

  std::unordered_map<long long, double> pot_per_day;   // key: YYYYMMDD, value: POT sum
  std::string line;
  while (std::getline(fin_pot, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double t = 0.0, p = 0.0;
    if (!(iss >> t >> p)) continue;
    if (!std::isfinite(t) || !std::isfinite(p)) continue;
    if (static_cast<long long>(t) < EVENT_CUT_START_SEC) continue; // Apply time cut for event-rate plot

    std::time_t tt = static_cast<std::time_t>(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    lt = *std::localtime(&tt);
#endif
    long long day_key =
      (static_cast<long long>(lt.tm_year + 1900) * 10000LL) +
      (static_cast<long long>(lt.tm_mon + 1) * 100LL) +
      static_cast<long long>(lt.tm_mday);

    pot_per_day[day_key] += p;
  }
  fin_pot.close();

  // 2) Load daily event counts from eventrate.tsv
  std::ifstream fin_nu(PATH_EVENT);
  if (!fin_nu) {
    ::Warning("dataqualityplot_withBSD",
              "No eventrate file found: %s (cannot build event-rate plot).",
              PATH_EVENT.c_str());
    return;
  }

  std::unordered_map<long long, long long> nu_per_day; // key: YYYYMMDD, value: event count
  while (std::getline(fin_nu, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double t = 0.0;
    if (!(iss >> t)) continue;
    if (!std::isfinite(t)) continue;
    if (static_cast<long long>(t) < EVENT_CUT_START_SEC) continue; // Apply time cut for event-rate plot

    std::time_t tt = static_cast<std::time_t>(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    lt = *std::localtime(&tt);
#endif
    long long day_key =
      (static_cast<long long>(lt.tm_year + 1900) * 10000LL) +
      (static_cast<long long>(lt.tm_mon + 1) * 100LL) +
      static_cast<long long>(lt.tm_mday);

    nu_per_day[day_key] += 1;
  }
  fin_nu.close();

  if (pot_per_day.empty() || nu_per_day.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No valid POT or neutrino data for event-rate plot.");
    return;
  }

  // 3) Collect valid days (POT > 0 and events > 0)
  std::vector<long long> day_keys_all;
  day_keys_all.reserve(pot_per_day.size());
  for (const auto& kv : pot_per_day) day_keys_all.push_back(kv.first);
  std::sort(day_keys_all.begin(), day_keys_all.end());

  std::vector<long long> day_keys;
  day_keys.reserve(day_keys_all.size());
  for (long long dk : day_keys_all) {
    auto it_pot = pot_per_day.find(dk);
    auto it_nu  = nu_per_day.find(dk);
    if (it_pot == pot_per_day.end()) continue;
    if (it_nu  == nu_per_day.end())  continue;
    const double pot = it_pot->second;
    const long long n_events = it_nu->second;
    if (!(pot > 0.0) || n_events <= 0) continue;
    day_keys.push_back(dk);
  }

  if (day_keys.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No valid daily points for event-rate plot (after cuts).");
    return;
  }

  const int nbins = static_cast<int>(day_keys.size());

  gSystem->mkdir(PATH_EVENTRATE_DIR.c_str(), true);
  gStyle->SetOptStat(0);

  TH1D h("h_eventrate",
         "Event rate of FROST;Date;Number of events / 10^{15} POT",
         nbins, 0.5, nbins + 0.5);
  h.SetStats(0);
  h.SetFillStyle(0);        // no fill
  h.SetLineColor(kBlack);
  h.SetLineWidth(2);
  h.SetMarkerStyle(0);      // no marker
  h.SetMarkerSize(0.0);

  // Fill histogram bins and set bin labels (MM/DD)
  for (int i = 0; i < nbins; ++i) {
    long long dk = day_keys[i];
    auto it_pot = pot_per_day.find(dk);
    auto it_nu  = nu_per_day.find(dk);
    if (it_pot == pot_per_day.end() || it_nu == nu_per_day.end()) continue;

    const double pot        = it_pot->second;
    const long long n_events = it_nu->second;
    if (!(pot > 0.0) || n_events <= 0) continue;

    const double rate     = static_cast<double>(n_events) * 1.0e15 / pot;
    const double err_rate = std::sqrt(static_cast<double>(n_events)) * 1.0e15 / pot;

    // std::cout << "rate" << rate << std::endl;
    // std::cout << "err_rate" << err_rate << std::endl;
    // std::cout << "n_events" << n_events << std::endl;
    // std::cout << "pot" << pot << std::endl;

    h.SetBinContent(i + 1, rate);
    h.SetBinError(i + 1, err_rate);
    // h.SetBinError(i + 1, 10);

    int y = static_cast<int>(dk / 10000LL);
    int m = static_cast<int>((dk / 100LL) % 100LL);
    int d = static_cast<int>(dk % 100LL);

    std::ostringstream lab;
    lab << std::setw(2) << std::setfill('0') << m
        << "/"
        << std::setw(2) << std::setfill('0') << d;
    h.GetXaxis()->SetBinLabel(i + 1, lab.str().c_str());
  }

  TCanvas c("c_eventrate", "Event rate of FROST", 1200, 800);
  c.SetGrid();
  c.SetLeftMargin(0.12);

  h.GetXaxis()->SetLabelSize(0.04);
  h.GetXaxis()->SetTitleOffset(1.5);
  h.GetXaxis()->LabelsOption("v");  // vertical labels if many days

  double ymax = h.GetMaximum();
  double ymin = h.GetMinimum();
  h.SetMinimum(ymin * 0.9);
  h.SetMaximum(ymax * 1.1); // 20% margin at the top

  h.Draw("E1");  // error bars only (no filled boxes)

  // Constant fit: y = p0 over all bins with content.
  // Note: we filled the histogram with SetBinContent(), so GetEntries()
  // is not reliable here. Instead, we explicitly check for non-empty bins.
  int n_fit_bins = 0;
  for (int ibin = 1; ibin <= h.GetNbinsX(); ++ibin) {
    double y  = h.GetBinContent(ibin);
    double ey = h.GetBinError(ibin);
    if (y > 0.0 && ey > 0.0 && std::isfinite(y) && std::isfinite(ey)) {
      ++n_fit_bins;
    }
  }


  TF1* f_const = nullptr;
  TPaveText* pt = nullptr;

  if (n_fit_bins > 0) {
    // constant function
    f_const = new TF1("f_const", "[0]",
                      h.GetXaxis()->GetXmin(),
                      h.GetXaxis()->GetXmax());
    f_const->SetLineWidth(2);
    f_const->SetLineColor(kRed);

    // Quiet fit (no printout)
    h.Fit(f_const, "Q0");  // "0" = do not draw automatically

    // Draw fit function on top of histogram
    f_const->Draw("SAME");

    // Stats box
    double p0      = f_const->GetParameter(0);
    double p0_err  = f_const->GetParError(0);
    double chi2    = f_const->GetChisquare();
    int    ndf     = f_const->GetNDF();

    pt = new TPaveText(0.80, 0.90, 0.99, 0.99, "NDC");
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetLineColor(1);
    pt->SetLineWidth(1);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->AddText(Form("#chi^{2} / ndf = %.1f / %d", chi2, ndf));
    pt->AddText(Form("p_{0} = %.3f #pm %.3f", p0, p0_err));
    pt->Draw();
  }

  c.SaveAs(PATH_EVENTRATE_PDF.c_str());

  ::Info("dataqualityplot_withBSD",
         "Saved event-rate histogram (per bunch) to %s",
         PATH_EVENTRATE_PDF.c_str());
}

// Build and save daily event-rate plot as a histogram (per spill):
//   X: date (one bin per day, labeled as "MM/DD")
//   Y: Number of events / 10^15 POT (for that day)
//   Here, each spill contributes at most one event (if any bunch fired).
//   Drawn with error bars only (no filled bins). Also performs a constant fit.
static void BuildAndSaveEventRatePlotSpill()
{
  // 1) Load daily POT from pot_points.tsv
  std::ifstream fin_pot(PATH_POT_POINTS);
  if (!fin_pot) {
    ::Warning("dataqualityplot_withBSD",
              "No pot_points file found: %s (cannot build per-spill event-rate plot).",
              PATH_POT_POINTS.c_str());
    return;
  }

  std::unordered_map<long long, double> pot_per_day;   // key: YYYYMMDD, value: POT sum
  std::string line;
  while (std::getline(fin_pot, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double t = 0.0, p = 0.0;
    if (!(iss >> t >> p)) continue;
    if (!std::isfinite(t) || !std::isfinite(p)) continue;
    if (static_cast<long long>(t) < EVENT_CUT_START_SEC) continue; // Apply time cut for event-rate plot

    std::time_t tt = static_cast<std::time_t>(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    lt = *std::localtime(&tt);
#endif
    long long day_key =
      (static_cast<long long>(lt.tm_year + 1900) * 10000LL) +
      (static_cast<long long>(lt.tm_mon + 1) * 100LL) +
      static_cast<long long>(lt.tm_mday);

    pot_per_day[day_key] += p;
  }
  fin_pot.close();

  // 2) Load daily event counts from eventrate_spill.tsv
  std::ifstream fin_nu(PATH_EVENT_SPILL);
  if (!fin_nu) {
    ::Warning("dataqualityplot_withBSD",
              "No per-spill events file found: %s (cannot build per-spill event-rate plot).",
              PATH_EVENT_SPILL.c_str());
    return;
  }

  std::unordered_map<long long, long long> nu_per_day; // key: YYYYMMDD, value: event count
  while (std::getline(fin_nu, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double t = 0.0;
    if (!(iss >> t)) continue;
    if (!std::isfinite(t)) continue;
    if (static_cast<long long>(t) < EVENT_CUT_START_SEC) continue; // Apply time cut for event-rate plot

    std::time_t tt = static_cast<std::time_t>(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    lt = *std::localtime(&tt);
#endif
    long long day_key =
      (static_cast<long long>(lt.tm_year + 1900) * 10000LL) +
      (static_cast<long long>(lt.tm_mon + 1) * 100LL) +
      static_cast<long long>(lt.tm_mday);

    nu_per_day[day_key] += 1;
  }
  fin_nu.close();

  if (pot_per_day.empty() || nu_per_day.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No valid POT or per-spill neutrino data for event-rate plot.");
    return;
  }

  // 3) Collect valid days (POT > 0 and events > 0)
  std::vector<long long> day_keys_all;
  day_keys_all.reserve(pot_per_day.size());
  for (const auto& kv : pot_per_day) day_keys_all.push_back(kv.first);
  std::sort(day_keys_all.begin(), day_keys_all.end());

  std::vector<long long> day_keys;
  day_keys.reserve(day_keys_all.size());
  for (long long dk : day_keys_all) {
    auto it_pot = pot_per_day.find(dk);
    auto it_nu  = nu_per_day.find(dk);
    if (it_pot == pot_per_day.end()) continue;
    if (it_nu  == nu_per_day.end())  continue;
    const double pot = it_pot->second;
    const long long n_events = it_nu->second;
    if (!(pot > 0.0) || n_events <= 0) continue;
    day_keys.push_back(dk);
  }

  if (day_keys.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No valid daily points for per-spill event-rate plot (after cuts).");
    return;
  }

  const int nbins = static_cast<int>(day_keys.size());

  gSystem->mkdir(PATH_EVENTRATE_DIR.c_str(), true);
  gStyle->SetOptStat(0);

  TH1D h("h_eventrate_spill",
         "Event rate of FROST (max 1 event per spill);Date;Number of events / 10^{15} POT",
         nbins, 0.5, nbins + 0.5);
  h.SetStats(0);
  h.SetFillStyle(0);        // no fill
  h.SetLineColor(kBlack);
  h.SetLineWidth(2);
  h.SetMarkerStyle(0);      // no marker
  h.SetMarkerSize(0.0);

  // Fill histogram bins and set bin labels (MM/DD)
  for (int i = 0; i < nbins; ++i) {
    long long dk = day_keys[i];
    auto it_pot = pot_per_day.find(dk);
    auto it_nu  = nu_per_day.find(dk);
    if (it_pot == pot_per_day.end() || it_nu == nu_per_day.end()) continue;

    const double pot        = it_pot->second;
    const long long n_events = it_nu->second;
    if (!(pot > 0.0) || n_events <= 0) continue;

    const double rate     = static_cast<double>(n_events) * 1.0e15 / pot;
    const double err_rate = std::sqrt(static_cast<double>(n_events)) * 1.0e15 / pot;

    h.SetBinContent(i + 1, rate);
    h.SetBinError(i + 1, err_rate);

    int y = static_cast<int>(dk / 10000LL);
    int m = static_cast<int>((dk / 100LL) % 100LL);
    int d = static_cast<int>(dk % 100LL);

    std::ostringstream lab;
    lab << std::setw(2) << std::setfill('0') << m
        << "/"
        << std::setw(2) << std::setfill('0') << d;
    h.GetXaxis()->SetBinLabel(i + 1, lab.str().c_str());
  }

  TCanvas c("c_eventrate_spill", "Event rate of FROST (max 1 event per spill)", 1200, 800);
  c.SetGrid();
  c.SetLeftMargin(0.12);

  h.GetXaxis()->SetLabelSize(0.04);
  h.GetXaxis()->LabelsOption("v");  // vertical labels if many days

  double ymax = h.GetMaximum();
  double ymin = h.GetMinimum();
  h.SetMinimum(ymin * 0.9);
  h.SetMaximum(ymax * 1.1); // 20% margin at the top

  h.Draw("E1");  // error bars only (no filled boxes)

  // Constant fit: y = p0 over all bins with content.
  int n_fit_bins = 0;
  for (int ibin = 1; ibin <= h.GetNbinsX(); ++ibin) {
    double y  = h.GetBinContent(ibin);
    double ey = h.GetBinError(ibin);
    if (y > 0.0 && ey > 0.0 && std::isfinite(y) && std::isfinite(ey)) {
      ++n_fit_bins;
    }
  }

  TF1* f_const = nullptr;
  TPaveText* pt = nullptr;

  if (n_fit_bins > 0) {
    f_const = new TF1("f_const_spill", "[0]",
                      h.GetXaxis()->GetXmin(),
                      h.GetXaxis()->GetXmax());
    f_const->SetLineWidth(2);
    f_const->SetLineColor(kRed);

    // Quiet fit (no printout)
    h.Fit(f_const, "Q0");  // "0" = do not draw automatically

   // Draw fit function on top of histogram
   f_const->Draw("SAME");
   // Stats box
   double p0      = f_const->GetParameter(0);
   double p0_err  = f_const->GetParError(0);
   double chi2    = f_const->GetChisquare();
   int    ndf     = f_const->GetNDF();
   pt = new TPaveText(0.80, 0.90, 0.99, 0.99, "NDC");
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetLineColor(1);
   pt->SetLineWidth(1);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.03);
   pt->AddText(Form("#chi^{2} / ndf = %.1f / %d", chi2, ndf));
   pt->AddText(Form("p_{0} = %.3f #pm %.3f", p0, p0_err));
   pt->Draw();
 }
 c.SaveAs(PATH_EVENTRATE_PDF_SPILL.c_str());
 ::Info("dataqualityplot_withBSD",
        "Saved per-spill event-rate histogram to %s",
        PATH_EVENTRATE_PDF_SPILL.c_str());
 }

// ------------------------
// One iteration
// ------------------------

static void OneIteration(int refresh_ms)
{
  // Ensure output base dir exists
  gSystem->mkdir(PATH_OUT_BASE.c_str(), true);

  // 1) Load all BSD spills into memory
  std::vector<BsdSpill> bsd_spills;
  std::unordered_map<int, std::vector<size_t>> index_by_mod;
  LoadAllBsdSpills(bsd_spills, index_by_mod);
  if (bsd_spills.empty()) {
    ::Warning("dataqualityplot_withBSD",
              "No BSD spills loaded. Sleeping...");
    gSystem->Sleep(refresh_ms);
    return;
  }

  // 2) List lightyield files
  auto ly_files = ListLightyieldFiles(PATH_LY_DIR);
  if (ly_files.empty()) {
    ::Info("dataqualityplot_withBSD",
           "No *_lightyield.root files in %s", PATH_LY_DIR.c_str());
    // Even in this case we can still draw the Delivered-only plot
    BuildAndSaveAccumulatedPotPlot(bsd_spills);
    gSystem->Sleep(refresh_ms);
    return;
  }

  // 3) Load acquired bunch rules
  std::vector<AcquiredBunchRule> acq_rules;
  LoadAcquiredBunchRules(acq_rules);

  // In the updated version we always rebuild the POT and eventrate caches
  // from all available (ready) lightyield files. This guarantees that the
  // accumulated POT and event-rate plots are always based on the full dataset,
  // not only on newly added files.

  std::vector<std::pair<double,double>> all_points;
  std::vector<double>                   all_nu_times;        // per-bunch counting
  std::vector<double>                   all_nu_times_spill;  // max 1 event per spill

  // 4) Process all "ready" lightyield files
  std::unordered_set<size_t> matched_bsd_indices;
  matched_bsd_indices.reserve(bsd_spills.size());

  for (const auto& fi : ly_files) {
    if (!IsReadyLightyieldFile(fi)) continue;

    ProcessLightyieldFileForPOT(fi.path,
                                bsd_spills,
                                index_by_mod,
                                all_points,
                                all_nu_times,
                                all_nu_times_spill,
                                acq_rules,
                                &matched_bsd_indices);
  }

  // 5) Rewrite cache files so they represent all current data
  AppendPotPoints(all_points);
  // Per-bunch events (original definition)
  AppendNuEvents(all_nu_times);
  // Per-spill events (max 1 event per spill)
  AppendNuEventsSpill(all_nu_times_spill);
  // 10-minute interval CSV output (Delivered/Recorded spills & POT)
  BuildAndSavePotSpillCsv10Min(bsd_spills);
  // Missed spills CSV output
  BuildAndSaveMissedSpillsCsv(bsd_spills, matched_bsd_indices);

  // 6) Build & save Accumulated POT plot (Delivered + Recorded)
  BuildAndSaveAccumulatedPotPlot(bsd_spills);

  // 7) Build & save daily event-rate plots
  //    (1) Per-bunch counting (original definition)
  //    (2) Per-spill counting (max 1 event per spill)
  BuildAndSaveEventRatePlot();
  BuildAndSaveEventRatePlotSpill();

  gSystem->Sleep(refresh_ms);
}

// ------------------------
// main
// ------------------------

int main(int argc, char** argv)
{
  // Default: infinite loop; refresh 10 min
  int max_loops  = -1;
  int refresh_ms = 600000; // 10 min

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if ((a == "-n" || a == "--loops") && i+1 < argc) {
      max_loops = std::stoi(argv[++i]);
    } else if ((a == "-r" || a == "--refresh") && i+1 < argc) {
      refresh_ms = std::stoi(argv[++i]);
    } else if (a == "-h" || a == "--help") {
      std::cout
        << "Usage: " << argv[0] << " [-n loops] [-r refresh_ms]\n"
        << "  -n, --loops     Number of iterations; negative for infinite (default -1)\n"
        << "  -r, --refresh   Refresh interval in milliseconds (default 600000 = 10 min)\n";
      return 0;
    }
  }

  int iter = 0;
  while (max_loops < 0 || iter < max_loops) {
    ++iter;
    OneIteration(refresh_ms);
  }

  ::Info("dataqualityplot_withBSD", "Stopped.");
  return 0;
}
