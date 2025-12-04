//calibration.cpp
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
#include <regex>
#include <filesystem>
#include <chrono>
#include <limits>
#include <iomanip>
#include <cctype>
#include <thread>
#include <atomic>
#include <csignal>
#include <ctime>

#include "/home/nu/notani/FROST_monitor/config/config.hpp"

// graceful shutdown flag
static std::atomic<bool> g_stop{false};
static void handle_sigint(int) { g_stop = true; }

namespace fs = std::filesystem;

// ---------------------------------------------
// Core calibration function (converted from macro)
// ---------------------------------------------
void calibrayraw_(TString filename, TString chmapPath, TString outputcsv, TString outputplot) {
  // Batch mode (no GUI)
  gROOT->SetBatch(kTRUE);

  // ---------- Input constants ----------
  const Int_t N_RAYRAW = FrostmonConfig::N_RAYRAW;          // plane count (1..11)
  const Int_t N_CH_PER_PLANE = FrostmonConfig::N_CH_PER_PLANE;    // channels per plane (0..31)
  const Int_t HIST_NBIN = FrostmonConfig::HIST_NBIN;
  const Double_t HIST_XMIN =  FrostmonConfig::HIST_XMIN;
  const Double_t HIST_XMAX =  FrostmonConfig::HIST_XMAX;

  // Baseline estimation
  const Int_t BL_WIN = FrostmonConfig::BL_WIN;
  const Int_t BL_START_MAX = FrostmonConfig::BL_START_MAX;
  const Int_t BL_STEP = FrostmonConfig::BL_STEP;

  // Integration windows
  const Int_t INT_LEFT  = FrostmonConfig::INT_LEFT;
  const Int_t INT_RIGHT = FrostmonConfig::INT_RIGHT;

  // Fixed window for 0pe: [FIX_START, FIX_START + FIX_RANGE)
  const Int_t FIX_START = FrostmonConfig::FIX_START;
  const Int_t FIX_RANGE = FrostmonConfig::FIX_RANGE;

  // TSpectrum peak search (for 1pe)
  const Int_t    PEAKS_MAX    = FrostmonConfig::PEAKS_MAX;
  const Double_t SPEC_SIGMA   = FrostmonConfig::SPEC_SIGMA;
  const Double_t SPEC_THRESH  = FrostmonConfig::SPEC_THRESH;

  // ADC thresholds (exclude unphysical low values)
  const double ADC_MIN = FrostmonConfig::ADC_MIN;

  // ---------- chmap load ----------
  std::map<long long, int> cabMap; // key = ((long long)rr << 32) | (unsigned)lc
  {
    std::ifstream fin(chmapPath.Data());
    if (fin) {
      std::string line;
      bool first = true;
      while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') { first = false; continue; }
        if (first) { first = false; continue; } // skip header
        std::istringstream iss(line);
        int cab, rr, lc;
        if (!(iss >> cab >> rr >> lc)) continue;
        long long key = ( (long long)rr << 32 ) | (unsigned)lc;
        if (cabMap.find(key) == cabMap.end()) cabMap[key] = cab;
      }
    } else {
      std::cerr << "[WARN] Cannot open chmap: " << chmapPath << "\n"
                << "      Missing/undefined channels will be skipped.\n";
    }
  }
  auto cabOf = [&](int rr, int lc)->int {
    long long key = ( (long long)rr << 32 ) | (unsigned)lc;
    auto it = cabMap.find(key);
    if (it == cabMap.end()) return -1;
    return it->second;
  };

  // ---------- Input ROOT ----------
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

  // ---------- Branches ----------
  std::vector<std::vector<double>> *waveform = nullptr;  // [ch][sample]
  std::vector<int> *plane = nullptr;                      // [ch] -> 0..10 (RAYRAW#-1)
  t->SetBranchAddress("waveform", &waveform);
  if (t->GetBranch("plane")) t->SetBranchAddress("plane", &plane);
  t->GetEntry(0);

  // ---------- Histograms (plane x channel) ----------
  // hInt1: variable window (peak±(20,25)) → 1pe fit
  // hInt0: fixed window ([1..FIX_RANGE]) → 0pe fit
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

  // ---------- Integration loop ----------
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
      if (!h1 || !h0) continue; // chmap missing → skip

      const std::vector<double> &wf = (*waveform)[ch];
      const Int_t ns = (Int_t) wf.size();
      if (ns <= 0) continue;

      // ---- Baseline (minimum of sliding-mean windows, ignoring values below ADC_MIN_BASELINE) ----
      Int_t startMaxEff = std::max(0, ns - BL_WIN);
      startMaxEff = std::min(startMaxEff, BL_START_MAX);

      double baseline = 0.0;
      {
        std::vector<double> bl_values;
        bl_values.reserve((startMaxEff / BL_STEP) + 2);

        for (Int_t start = 1; start <= startMaxEff; start += BL_STEP) {
          double sum   = 0.0;
          int    count = 0;

          const Int_t end   = std::min(ns, start + BL_WIN);
          const Int_t width = end - start;
          if (width <= 0) continue;

          for (Int_t j = start; j < end; ++j) {
            const double v = wf[j];
            if (v <= ADC_MIN) continue;  // skip samples below ADC threshold
            sum += v;
            ++count;
          }

          if (count == 0) continue;  // this window has no valid samples
          bl_values.push_back(sum / static_cast<double>(count));
        }

        // if no valid window exists, skip this channel/event
        if (bl_values.empty()) continue;
        std::sort(bl_values.begin(), bl_values.end());
        baseline = bl_values.front();
      }


      // ---- Global maximum sample (variable window center) ----
      Int_t peak_idx = 0;
      double peak_val = -1e300;
      for (Int_t j = 0; j < ns; ++j) {
        if (wf[j] > peak_val) { peak_val = wf[j]; peak_idx = j; }
      }

      // ---- Variable-window integration around the peak (ignore samples below ADC_MIN_INTEGRAL) ----
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
          if (jmin < 0)      jmin = 0;
          if (jmax > ns - 1) jmax = ns - 1;
        }

        double integral1 = 0.0;
        for (Int_t j = jmin; j < jmax; ++j) {
          const double v = wf[j];
          if (v <= ADC_MIN) continue;  // skip samples below ADC threshold
          integral1 += (v - baseline);
        }

        if (peak_idx > INT_LEFT && peak_idx < ns - 1 - INT_RIGHT) {
          h1->Fill(integral1);
        }
      }


      // ---- Fixed-window integration [FIX_START .. FIX_START + FIX_RANGE) (ignore samples below ADC_MIN_INTEGRAL) ----
      {
        const Int_t jmax_fix = std::min(ns, FIX_START + FIX_RANGE);
        if (jmax_fix > FIX_START) {
          double integral0 = 0.0;
          for (Int_t j = FIX_START; j < jmax_fix; ++j) {
            const double v = wf[j];
            if (v <= ADC_MIN) continue;  // skip samples below ADC threshold
            integral0 += (v - baseline);
          }
          h0->Fill(integral0);
        }
      }

    }
  }

  // ========= Store fit results for CSV =========
  std::vector<std::vector<double>> mean0(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean0e(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean1(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<double>> mean1e(N_RAYRAW, std::vector<double>(N_CH_PER_PLANE, NAN));
  std::vector<std::vector<int>> npeak(N_RAYRAW, std::vector<int>(N_CH_PER_PLANE, -1));

  // ---------- Plot & fit (0pe: fixed window, ±10 around max bin) ----------
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

  // ---------- Plot & fit (1pe: peak search → right-side peak ±15) ----------
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

      TSpectrum spec(PEAKS_MAX);
      Int_t nfound = spec.Search(h, SPEC_SIGMA, "", SPEC_THRESH);
      npeak[p][lc] = nfound;

      std::vector<double> xs;
      for (Int_t i = 0; i < nfound && i < 2; ++i) xs.push_back(spec.GetPositionX()[i]);
      std::sort(xs.begin(), xs.end());

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

  // ---------- CSV output ----------
  std::ofstream fout(outputcsv);
  if (!fout) {
    std::cerr << "[ERROR] Cannot create CSV: " << outputcsv.Data() << std::endl;
  } else {
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
        if (m0 > m1v || m1e > 3.0 || npeak[p][lc] < 2 || m1v > 100) badflag = 1;

        fout << cab << ","
             << m0  << "," << m0e << ","
             << m1v << "," << m1e << ","
             << FIX_RANGE << "," << badflag << "\n";
      }
    }
    fout.close();
    std::cout << "[INFO] Saved CSV: " << outputcsv.Data() << std::endl;
  }
}

// ---------------------------------------------
// Helpers for batch driver
// ---------------------------------------------
struct FileInfo {
  fs::path path;
  std::string base;   // runNNNNN[_start_end]
  std::uintmax_t size = 0; // current file size in bytes
  int run = -1;
  long long segStart = -1; // parsed from _start_end
  std::time_t mtime = 0;
  bool parsed = false;
};

static bool fileExists(const fs::path& p) {
  std::error_code ec;
  return fs::exists(p, ec);
}

// C++17-safe conversion: file_time_type -> time_t
static std::time_t filetime_to_time_t(std::filesystem::file_time_type ftime) {
  using namespace std::chrono;
  auto sctp = time_point_cast<system_clock::duration>(
      ftime - std::filesystem::file_time_type::clock::now() + system_clock::now());
  return system_clock::to_time_t(sctp);
}

// Try to parse run number and segment start from filename
// Supported:
//   run00097_0_9999.root  -> run=97, segStart=0
//   run00097.root         -> run=97, segStart=-1
static FileInfo parseInfo(const fs::directory_entry& de) {
  FileInfo info;
  info.path = de.path();
  info.base = info.path.stem().string(); // without .root

  // mtime
  std::error_code ec;
  auto ftime = fs::last_write_time(info.path, ec);
  if (!ec) {
    info.mtime = filetime_to_time_t(ftime);
  } else {
    info.mtime = 0;
  }

  // size
  std::error_code ec2;
  auto fsize = fs::file_size(info.path, ec2);
  if (!ec2) {
    info.size = fsize;
  } else {
    info.size = 0;
  }

  const std::string name = info.path.filename().string();
  std::regex reFull(R"(run(\d{5})_(\d+)_\d+\.root)");
  std::regex reSimple(R"(run(\d{5})\.root)");
  std::smatch m;

  if (std::regex_match(name, m, reFull)) {
    info.run = std::stoi(m[1].str());
    info.segStart = std::stoll(m[2].str());
    info.parsed = true;
  } else if (std::regex_match(name, m, reSimple)) {
    info.run = std::stoi(m[1].str());
    info.segStart = -1;
    info.parsed = true;
  } else {
    info.parsed = false;
  }
  return info;
}

// Sort by (run asc, segStart asc), fallback to mtime asc, then name
static bool infoLess(const FileInfo& a, const FileInfo& b) {
  if (a.parsed && b.parsed) {
    if (a.run != b.run) return a.run < b.run;
    if (a.segStart != b.segStart) return a.segStart < b.segStart;
    if (a.mtime != b.mtime) return a.mtime < b.mtime;
    return a.base < b.base;
  }
  // If only one parsed, put parsed first
  if (a.parsed != b.parsed) return a.parsed;
  // Neither parsed → mtime then name
  if (a.mtime != b.mtime) return a.mtime < b.mtime;
  return a.base < b.base;
}

// Trim leading/trailing spaces
static std::string trim(std::string s) {
  auto issp = [](unsigned char c){ return std::isspace(c); };
  while (!s.empty() && issp(s.front())) s.erase(s.begin());
  while (!s.empty() && issp(s.back()))  s.pop_back();
  return s;
}

struct RunFileRule {
  int run_max;           // This rule applies to runs with run <= run_max
  std::string filename;  // File name for this rule
};

// Load run→file rules from a text file.
// Expected format per line (whitespace separated):
//   <run_max> <filename>
// Example:
//   5      chmap_20251009.txt
//   20     chmap_20251122.txt
static std::vector<RunFileRule> loadRunFileRules(const std::filesystem::path& rulePath) {
  std::vector<RunFileRule> rules;
  std::ifstream fin(rulePath);
  if (!fin) return rules;  // If file cannot be opened, return empty

  std::string line;
  while (std::getline(fin, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    int runmax;
    std::string fname;
    if (!(iss >> runmax >> fname)) continue;

    rules.push_back(RunFileRule{runmax, fname});
  }

  // Ensure rules are sorted by run_max in ascending order
  std::sort(rules.begin(), rules.end(),
            [](const RunFileRule& a, const RunFileRule& b) {
              return a.run_max < b.run_max;
            });

  return rules;
}

// Choose a file name for a given run using the rule list.
// If no rule matches (or run is invalid), the defaultFile is returned.
static std::string chooseFileForRun(const std::vector<RunFileRule>& rules,
                                    int run,
                                    const std::string& defaultFile)
{
  if (run < 0 || rules.empty()) return defaultFile;
  for (const auto& r : rules) {
    if (run <= r.run_max) return r.filename;
  }
  return defaultFile;
}


static void printUsage(const char* prog) {
  std::cout <<
    "Usage: " << prog << " [options]\n"
    "Options:\n"
    "  -h, --help           Show this help and exit\n"
    "  --base <DIR>         Base directory (contains rootfile/, calibration/, chmap/)\n"
    "                       Default: /group/nu/ninja/work/otani/FROST_beamdata/test\n"
    "  --chmap <FILENAME>   Chmap filename under chmap/ (e.g., chmap_20251122.txt)\n"
    "                       Default: chmap_20251122.txt\n"
    "  --limit <N>          Process at most N files (0 = no limit; default 0)\n"
    "  --dry-run <0|1>      List files to process without running (default 0)\n"
    "  --watch <SEC>        Watch mode: rescan every SEC seconds until Ctrl-C (default 60)\n";
}

// One-shot scan and process. Returns number of processed files.
static int scanAndProcess(const std::string& base,
                          const std::string& chmapFileCli,
                          int limit,
                          bool dryRun)
{
  namespace fs = std::filesystem;

  fs::path rootDir = fs::path(base) / "rootfile";
  fs::path csvDir  = fs::path(base) / "calibration" / "calibresult";
  fs::path plotDir = fs::path(base) / "calibration" / "plot";

  // Determine the default chmap:
  //  - If the CLI option --chmap is given, always use that file.
  //  - Otherwise, use the config default FrostmonConfig::CHMAP_FILE.
  std::string defaultChmap = chmapFileCli.empty()
                           ? FrostmonConfig::CHMAP_FILE
                           : chmapFileCli;

  // Load chmap rules only when the chmap was not fixed by CLI.
  std::vector<RunFileRule> chmapRules;
  if (chmapFileCli.empty()) {
    fs::path rulePath = fs::path(base) / "chmap" / FrostmonConfig::CHMAP_RULE_FILE;
    chmapRules = loadRunFileRules(rulePath);
    if (!chmapRules.empty()) {
      std::cout << "[CFG] Loaded " << chmapRules.size()
                << " chmap rule(s) from " << rulePath << "\n";
    } else {
      std::cout << "[CFG] No chmap rule file (" << rulePath
                << "). Using default chmap=" << defaultChmap << "\n";
    }
  } else {
    std::cout << "[CFG] --chmap specified. Chmap rules are disabled. chmap="
              << defaultChmap << "\n";
  }

  // Threshold: only files >= 10 KB and "stable" for >= 60 s are considered
  constexpr std::uintmax_t MIN_SIZE_BYTES = 10ull * 1024ull;
  constexpr std::time_t    STABLE_SEC     = 60;

  // Ensure output directories exist
  std::error_code ec;
  fs::create_directories(csvDir, ec);
  fs::create_directories(plotDir, ec);

  // Collect .root files
  std::vector<FileInfo> files;
  if (!fs::exists(rootDir)) {
    std::cerr << "[ERROR] rootDir does not exist: " << rootDir << "\n";
    return 0;
  }
  for (const auto& de : fs::directory_iterator(rootDir)) {
    if (!de.is_regular_file()) continue;
    if (de.path().extension() != ".root") continue;
    files.push_back(parseInfo(de));
  }
  if (files.empty()) {
    std::cout << "[INFO] No .root files in " << rootDir << "\n";
    return 0;
  }

  // Sort by run / segment / mtime
  std::sort(files.begin(), files.end(), infoLess);

  // Current time for "stability" judgement
  std::time_t now_t = std::time(nullptr);

  // Filter "unprocessed" (CSV missing) and "stable" files
  std::vector<FileInfo> todo;
  for (const auto& fi : files) {
    // Size & stability check
    bool stable = false;
    if (fi.size >= MIN_SIZE_BYTES && fi.mtime != 0) {
      if (std::difftime(now_t, fi.mtime) >= STABLE_SEC) {
        stable = true;
      }
    }
    if (!stable) {
      // Probably still being written (or too small) → skip in this round
      continue;
    }

    const std::string infile = fi.base; // e.g., run00097_0_9999
    fs::path outCsv = csvDir / ("calib_" + infile + ".csv");
    if (!fileExists(outCsv)) todo.push_back(fi);
  }

  if (todo.empty()) {
    std::cout << "[INFO] Nothing to do.\n";
    return 0;
  }

  std::cout << "[INFO] To process (" << todo.size() << ") "
            << "(only files >= 10KB and stable for >= 60s are considered):\n";
  for (const auto& fi : todo) {
    std::cout << "  - " << fi.path.filename().string()
              << "  [run="      << (fi.parsed ? std::to_string(fi.run)      : "NA")
              << ", segStart="  << (fi.parsed ? std::to_string(fi.segStart) : "NA")
              << ", size="      << (fi.size / (1024.0 * 1024.0)) << "MB"
              << "]\n";
  }

  if (dryRun) {
    std::cout << "[INFO] Dry-run: listing only.\n";
    return 0;
  }

  int processed = 0;
  for (const auto& fi : todo) {
    if (g_stop) break;                 // allow graceful stop
    if (limit > 0 && processed >= limit) break;

    const std::string infile = fi.base;

    // Select chmap file for this run
    std::string chmapFileThisRun = defaultChmap;
    if (!chmapRules.empty() && fi.parsed) {
      chmapFileThisRun = chooseFileForRun(chmapRules, fi.run, defaultChmap);
    }

    const TString filename_inroot  =
        TString((fs::path(base) / "rootfile" / (infile + ".root")).string());
    const TString filename_outcsv  =
        TString((fs::path(base) / "calibration" / "calibresult"
                 / ("calib_" + infile + ".csv")).string());
    const TString filename_outplot =
        TString((fs::path(base) / "calibration" / "plot"
                 / ("calib_" + infile)).string());
    const TString filename_chmap   =
        TString((fs::path(base) / "chmap" / chmapFileThisRun).string());

    std::cout << "[RUN] " << infile
              << " (run=" << (fi.parsed ? std::to_string(fi.run) : std::string("NA"))
              << ", chmap=" << chmapFileThisRun << ")\n";

    calibrayraw_(filename_inroot, filename_chmap, filename_outcsv, filename_outplot);
    ++processed;
  }

  std::cout << "[DONE] Processed " << processed << " file(s).\n";
  return processed;
}

// ---------------------------------------------
// Batch driver (main)
// ---------------------------------------------
int main(int argc, char** argv) {
  std::string base = FrostmonConfig::OUTPUT_DIR;
  // Empty chmapFileCli means "--chmap" was not given; config default + rules are used.
  std::string chmapFileCli;
  int limit = 0;
  bool dryRun = false;
  int watchSec = 60;  // 0 => one-shot

  // handle Ctrl-C
  std::signal(SIGINT, handle_sigint);
  std::signal(SIGTERM, handle_sigint);

  // ---- parse CLI ----
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") { printUsage(argv[0]); return 0; }

    auto next = [&](int& i)->std::string {
      if (i + 1 < argc) return std::string(argv[++i]);
      std::cerr << "Missing value after " << a << "\n";
      printUsage(argv[0]);
      std::exit(1);
    };

    if (a == "--base")        base       = next(i);
    else if (a == "--chmap")  chmapFileCli = next(i);
    else if (a == "--limit")  limit     = std::stoi(next(i));
    else if (a == "--dry-run"){ std::string v = next(i); dryRun = (v=="1"||v=="true"||v=="yes"); }
    else if (a == "--watch")  watchSec  = std::stoi(next(i));
    else {
      std::cerr << "Unknown option: " << a << "\n";
      printUsage(argv[0]);
      return 1;
    }
  }

// This is the effective default chmap used when no per-run rule overrides it.
  const std::string effectiveChmap =
      chmapFileCli.empty() ? FrostmonConfig::CHMAP_FILE : chmapFileCli;

  std::cout << "[CFG] base=" << base
            << " chmap=" << effectiveChmap
            << " limit=" << limit
            << " dryRun=" << (dryRun?1:0)
            << " watch=" << watchSec << "s\n";

  if (watchSec <= 0) {
    // one-shot
    (void)scanAndProcess(base, chmapFileCli, limit, dryRun);
    return 0;
  }

  // watch mode
  std::cout << "[WATCH] Start watching every " << watchSec << " seconds. Press Ctrl-C to stop.\n";

  while (!g_stop) {
    // show timestamp
    std::time_t now = std::time(nullptr);
    std::cout << "\n[WATCH] Scan at " << std::put_time(std::localtime(&now), "%F %T") << std::endl;

    (void)scanAndProcess(base, chmapFileCli, /*limit=*/0, /*dryRun=*/dryRun);
    // In watch mode, limit is ignored (always process all pending)

    // sleep in small chunks so Ctrl-C is responsive
    for (int s = 0; s < watchSec && !g_stop; ++s) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
  }

  std::cout << "\n[WATCH] Stopped.\n";
  return 0;
}
