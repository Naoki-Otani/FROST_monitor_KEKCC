// dataqualityplot_withBSD.cpp
// Purpose:
//   Standalone C++17 program using ROOT to generate:
//     (1) Accumulated POT vs time plot
//         - Delivered POT from BSD (black curve, "Delivered")
//         - POT recorded by FROST (red curve, "Recorded by FROST")
//     (2) Event-rate plot of neutrino-like events
//         - Neutrino-like event defined as: for each bunch (0–7),
//           max(lightyield[ch][b]) >= 10 → counts as 1 event
//         - Event rate = (#events per day) / (Recorded POT of the day)
//         - Units: events / 1×10^14 POT
//         - Output: eventrate_plot/eventrate_frost.pdf
//     (3) JS summary files for HTML dashboard
//
// BSD Input:
//   Directory: /group/nu/ninja/work/otani/FROST_beamdata/test/bsd/
//   Tree name: "bsd"
//   Used branches:
//     * spillnum          (full 32-bit spill number)
//     * trg_sec[3]        → use trg_sec[2] (Unix time)
//     * good_spill_flag   (0=bad; 1,–1,2,3,–3,100 are good beams)
//     * ct_np[5][9]       → use ct_np[4][0] (POT)
//
// Lightyield Input:
//   Directory: /group/nu/ninja/work/otani/FROST_beamdata/test/rootfile_aftercalib/
//   Tree name: "tree"
//   Used branches:
//     * spillnum          (15-bit; 0–32767 wraps)
//     * unixtime[272]     → use unixtime[0]
//     * lightyield[272][8]  // used for neutrino-like event counting
//
// Matching rule (per reconstructed event):
//   1) key = (spillnum & 0x7FFF)                 // 15-bit modulo
//   2) Find BSD spills sharing this key
//   3) Among them, select the one satisfying:
//         |bsd.trg_sec − round(unixtime[0])| ≤ 5 seconds
//      choosing the minimal |Δt|
//   4) Each BSD spill contributes its POT at most once per lightyield file
//   5) Neutrino-like events are counted per bunch (0–7)
//
// Incremental update strategy:
//   - BSD files: always fully re-read (small data volume)
//   - Lightyield files: incremental processing using:
//       * processed_lightyield.tsv  → records (path, mtime)
//       * pot_points.tsv            → stores accumulated POT increments
//       * nu_events.tsv             → stores times of neutrino-like events
//
// Outputs (all under /group/nu/ninja/work/otani/FROST_beamdata/test/dataquality_withBSD):
//
//   1. Cache files:
//        pot_points.tsv
//        processed_lightyield.tsv
//        nu_events.tsv
//
//   2. PDF plots:
//        accumulated_pot_withBSD.pdf
//        eventrate_plot/eventrate_frost.pdf
//
//   3. JS summary files (for the web dashboard):
//        pot_info/period.js
//        pot_info/delivered_spills.js
//        pot_info/recorded_spills.js
//        pot_info/delivered_pot.js
//        pot_info/recorded_pot.js
//        pot_info/efficiency.js
//
// Build example:
//   g++ -O2 -std=c++17 dataqualityplot_withBSD.cpp -o dataqualityplot_withBSD \
//       $(root-config --cflags --libs)
//
// Run examples:
//   ./dataqualityplot_withBSD                   // infinite loop, refresh every 10 min
//   ./dataqualityplot_withBSD -n 1              // run once
//   ./dataqualityplot_withBSD -n -1 -r 5000     // infinite loop, refresh every 5 sec



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

static const std::string PATH_POT_POINTS =
  PATH_POT_PLOT+ "pot_points.tsv";
static const std::string PATH_POT_PROC   =
  PATH_POT_PLOT + "processed_lightyield.tsv";
static const std::string PATH_OUT_PDF    =
  PATH_POT_PLOT + "accumulated_pot_withBSD.pdf";

static const std::string PATH_EVENTRATE_DIR =
  PATH_OUT_BASE + "eventrate_plot/";
static const std::string PATH_NU_EVENTS  =
  PATH_EVENTRATE_DIR + "neutrino_events.tsv";
static const std::string PATH_EVENTRATE_PDF =
  PATH_EVENTRATE_DIR + "eventrate_plot.pdf";

static const char* BSD_TREE_NAME  = "bsd";
static const char* LY_TREE_NAME   = "tree";

// file readiness
static const Long_t MIN_LY_SIZE_BYTES = 10L * 1024; // 10 KB
static const Long_t STABLE_SEC        = 60;               // consider stable if mtime older than 60 s

// matching
static const int SPILL_MOD = FrostmonConfig::SPILL_MOD; // 15-bit modulo for spillnum matching

static const int MAX_TIME_DIFF = FrostmonConfig::MAX_TIME_DIFF; // seconds

static const Double_t LIGHTMAX_MIN = FrostmonConfig::LIGHTMAX_MIN; //threshold for neutrino event

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

// Format POT as scientific notation for HTML:
// Example: 4.41 × 10^19 becomes "4.41 &times; 10<sup>19</sup>"
static std::string FormatSciForHtml(double val)
{
  if (!std::isfinite(val) || val <= 0.0) return "0";

  int exp = static_cast<int>(std::floor(std::log10(val)));
  double mant = val / std::pow(10.0, exp);

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << mant
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

struct BsdSpill {
  int spillnum_full = 0;   // full spillnumber
  int spillnum_mod  = 0;   // low 15 bits
  int good_spill_flag = 0;
  int trg_sec = 0;         // trg_sec[2]
  double pot = 0.0;        // ct_np[4][0]
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
    Double_t ct_np[5][9];

    t->SetBranchStatus("spillnum", 1);
    t->SetBranchAddress("spillnum", &spillnum);

    t->SetBranchStatus("trg_sec", 1);
    t->SetBranchAddress("trg_sec", trg_sec_arr);

    t->SetBranchStatus("good_spill_flag", 1);
    t->SetBranchAddress("good_spill_flag", &good_spill_flag);

    t->SetBranchStatus("ct_np", 1);
    t->SetBranchAddress("ct_np", ct_np);

    const Long64_t nent = t->GetEntries();
    for (Long64_t ie = 0; ie < nent; ++ie) {
      t->GetEntry(ie);

      // If you want to include bad spills also, comment out next 2 lines
      if (good_spill_flag == 0) continue;

      BsdSpill s;
      s.spillnum_full   = spillnum;
      s.spillnum_mod    = spillnum & (SPILL_MOD - 1); // 15 bits
      s.good_spill_flag = good_spill_flag;
      s.trg_sec         = trg_sec_arr[2];             // trg_sec[2]
      s.pot             = ct_np[4][0];                // ct_np[4][0]

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
  if (new_points.empty()) return;

  std::ofstream fout(PATH_POT_POINTS, std::ios::app);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open pot_points file for append: %s",
              PATH_POT_POINTS.c_str());
    return;
  }

  // check header
  static bool header_written = false;
  static bool header_checked = false;

  if (!header_checked) {
    std::ifstream fin(PATH_POT_POINTS);
    header_written = false;
    if (fin) {
      std::string firstLine;
      if (std::getline(fin, firstLine)) {
        if (!firstLine.empty() && firstLine[0] == '#') header_written = true;
      }
    }
    header_checked = true;
  }

  if (!header_written) {
    fout << "#time pot_increment\n";
    header_written = true;
  }

  fout << std::setprecision(16);
  for (const auto& p : new_points) {
    fout << p.first << " " << p.second << "\n";
  }
}

// Append neutrino event timestamps (Unix time [s]) to PATH_NU_EVENTS.
// Each line:  time
static void AppendNuEvents(const std::vector<double>& nu_times)
{
  if (nu_times.empty()) return;

  // Ensure the parent directory exists before opening the file
  gSystem->mkdir(PATH_EVENTRATE_DIR.c_str(), true);
  std::ofstream fout(PATH_NU_EVENTS, std::ios::app);
  if (!fout) {
    ::Warning("dataqualityplot_withBSD",
              "Failed to open neutrino events file for append: %s",
              PATH_NU_EVENTS.c_str());
    return;
  }

  static bool header_written = false;
  static bool header_checked = false;

  if (!header_checked) {
    std::ifstream fin(PATH_NU_EVENTS);
    header_written = false;
    if (fin) {
      std::string firstLine;
      if (std::getline(fin, firstLine)) {
        if (!firstLine.empty() && firstLine[0] == '#') header_written = true;
      }
    }
    header_checked = true;
  }

  if (!header_written) {
    fout << "#time\n";
    header_written = true;
  }

  fout << std::setprecision(16);
  for (double t : nu_times) {
    fout << t << "\n";
  }
}

// ------------------------
// Process lightyield file -> find POT increments
// ------------------------

static void ProcessLightyieldFileForPOT(
    const std::string& path,
    const std::vector<BsdSpill>& bsd_spills,
    const std::unordered_map<int, std::vector<size_t>>& index_by_mod,
    std::vector<std::pair<double,double>>& out_new_points,
    std::vector<double>& out_new_nu_times)
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

  t->SetBranchStatus("spillnum",  1);
  t->SetBranchAddress("spillnum", &spillnum);

  t->SetBranchStatus("unixtime",  1);
  t->SetBranchAddress("unixtime", unixtime_arr);

  t->SetBranchStatus("lightyield", 1);
  t->SetBranchAddress("lightyield", lightyield_arr);

  const Long64_t nent = t->GetEntries();

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

    // Count "neutrino events" per bunch:
    // For each bunch (0..7), take the maximum lightyield over all channels.
    // If max(lightyield[ch][b]) >= 10, that bunch contributes 1 event.
    int n_nu_events_this_spill = 0;
    for (int b = 0; b < 8; ++b) {
      double ly_max_b = -1.0;
      for (int ch = 0; ch < 272; ++ch) {
        double v = lightyield_arr[ch][b];
        if (std::isfinite(v) && v > ly_max_b) {
          ly_max_b = v;
        }
        // if(v > LIGHTMAX_MIN){
        //   std::cout << "ch" << ch << std::endl;
        // }
      }
      if (ly_max_b >= LIGHTMAX_MIN) {
        ++n_nu_events_this_spill;
      }
    }

    // Each neutrino event is stored as one entry at this spill time
    for (int i = 0; i < n_nu_events_this_spill; ++i) {
      out_new_nu_times.push_back(static_cast<double>(s.trg_sec));
    }

    // For POT accumulation, avoid double counting the same BSD spill within this file
    if (used_bsd_indices.insert(best_idx).second) {
      out_new_points.emplace_back(static_cast<double>(s.trg_sec), s.pot);
      ++match_count;
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

  TCanvas c("c_pot", "Accumulated POT (Delivered vs FROST)", 1200, 800);
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

    TAxis* ax = gr_del->GetXaxis();
    ax->SetTimeDisplay(1);
    ax->SetTimeOffset(0, "local");
    ax->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    ax->SetLabelSize(0.03);
    ax->SetLabelOffset(0.02);
    ax->SetTitleOffset(1.3);

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
      ax->SetTitleOffset(1.3);
      gr_rec->Draw("AL");
    } else {
      gr_rec->Draw("L SAME");
    }
    gr_rec->SetLineColor(kRed);
    gr_rec->SetLineWidth(2);
    gr_rec->SetMarkerStyle(0);
  }

  // Legend
  TLegend leg(0.15, 0.75, 0.45, 0.85);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  if (hasDelivered) leg.AddEntry(gr_del, "Delivered", "l");
  if (hasRecorded)  leg.AddEntry(gr_rec, "Recorded by FROST", "l");
  leg.Draw();

  // Timestamp
  TDatime now;
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.03);
  lat.DrawLatex(0.7, 0.93,
                Form("updated: %04d-%02d-%02d %02d:%02d:%02d",
                     now.GetYear(), now.GetMonth(), now.GetDay(),
                     now.GetHour(), now.GetMinute(), now.GetSecond()));

  c.SaveAs(PATH_OUT_PDF.c_str());

  ::Info("dataqualityplot_withBSD",
         "Saved Accumulated POT plot (Delivered + Recorded) to %s",
         PATH_OUT_PDF.c_str());

  delete gr_del;
  delete gr_rec;
}

// ------------------------
// Build and save daily event-rate plot as a histogram:
//   X: date (one bin per day, labeled as "MM/DD")
//   Y: Number of events / 10^14 POT (for that day)
//   Drawn with error bars only (no filled bins).
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

  // 2) Load daily neutrino event counts from neutrino_events.tsv
  std::ifstream fin_nu(PATH_NU_EVENTS);
  if (!fin_nu) {
    ::Warning("dataqualityplot_withBSD",
              "No neutrino events file found: %s (cannot build event-rate plot).",
              PATH_NU_EVENTS.c_str());
    return;
  }

  std::unordered_map<long long, long long> nu_per_day; // key: YYYYMMDD, value: event count
  while (std::getline(fin_nu, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double t = 0.0;
    if (!(iss >> t)) continue;
    if (!std::isfinite(t)) continue;

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
         "Event rate of FROST;Date;Number of events / 10^{14} POT",
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

    const double rate     = static_cast<double>(n_events) * 1.0e14 / pot;
    const double err_rate = std::sqrt(static_cast<double>(n_events)) * 1.0e14 / pot;

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
  h.GetXaxis()->LabelsOption("v");  // vertical labels if many days

  h.Draw("E1");  // error bars only (no filled boxes)

  c.SaveAs(PATH_EVENTRATE_PDF.c_str());

  ::Info("dataqualityplot_withBSD",
         "Saved event-rate histogram to %s",
         PATH_EVENTRATE_PDF.c_str());
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
    // それでも Delivered だけは書けるので plot だけ実行
    BuildAndSaveAccumulatedPotPlot(bsd_spills);
    gSystem->Sleep(refresh_ms);
    return;
  }

  // 3) Load processed LY file list
  std::unordered_map<std::string, Long_t> seen;
  LoadProcessedLyFiles(seen);

  bool changed_seen = false;
  std::vector<std::pair<double,double>> new_points;
  std::vector<double>                   new_nu_times;

  // 4) Process only new/updated and "ready" LY files
  for (const auto& fi : ly_files) {
    if (!IsReadyLightyieldFile(fi)) continue;

    auto it = seen.find(fi.path);
    if (it != seen.end() && it->second == fi.mtime) {
      // already processed with same mtime
      continue;
    }

    ProcessLightyieldFileForPOT(fi.path,
                                bsd_spills,
                                index_by_mod,
                                new_points,
                                new_nu_times);
    seen[fi.path] = fi.mtime;
    changed_seen = true;
  }

  // 5) Update caches
  if (!new_points.empty()) {
    AppendPotPoints(new_points);
  }
  if (!new_nu_times.empty()) {
    AppendNuEvents(new_nu_times);
  }
  if (changed_seen) {
    SaveProcessedLyFiles(seen);
  }

  // 6) Build & save Accumulated POT plot (Delivered + Recorded)
  BuildAndSaveAccumulatedPotPlot(bsd_spills);

  // 7) Build & save daily event-rate plot
  BuildAndSaveEventRatePlot();

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
