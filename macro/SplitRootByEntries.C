// SplitRootByEntries.C
// ROOT macro to split a TTree ROOT file into chunks of fixed number of events.
//
// Example:
//   root -l -q 'SplitRootByEntries.C("run00032_0_86305_aftershift.root","tree","/path/to/outdir",10000)'
//
// Output files:
//   <outdir>/run00032_0_9999.root
//   <outdir>/run00032_10000_19999.root
//   ...
//   <outdir>/run00032_80000_86303.root   (last chunk may be smaller)

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <algorithm>

static std::string BasenameNoExt(const std::string& path) {
  // Remove directory part
  std::string base = path;
  const size_t slash = base.find_last_of("/\\");
  if (slash != std::string::npos) base = base.substr(slash + 1);

  // Remove .root extension if present
  const std::string ext = ".root";
  if (base.size() >= ext.size() &&
      base.compare(base.size() - ext.size(), ext.size(), ext) == 0) {
    base = base.substr(0, base.size() - ext.size());
  }
  return base;
}

static std::string JoinPath(const std::string& dir, const std::string& file) {
  if (dir.empty()) return file;
  if (dir.back() == '/' || dir.back() == '\\') return dir + file;
  return dir + "/" + file;
}

void SplitRootByEntries(const char* inFileName = "run00032_0_86305_aftershift.root",
                        const char* outDir     = ".",
                        const char* treeName   = "tree",
                        Long64_t chunkSize     = 10000)
{
  if (chunkSize <= 0) {
    std::cerr << "[ERROR] chunkSize must be > 0" << std::endl;
    return;
  }

  // Ensure output directory exists (mkdir -p behavior)
  if (gSystem->AccessPathName(outDir)) {
    const int ok = gSystem->mkdir(outDir, /*recursive=*/kTRUE);
    if (ok != 0) {
      std::cerr << "[ERROR] Cannot create output directory: " << outDir << std::endl;
      return;
    }
  }

  // Open input
  TFile* fin = TFile::Open(inFileName, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inFileName << std::endl;
    return;
  }

  TTree* tin = dynamic_cast<TTree*>(fin->Get(treeName));
  if (!tin) {
    std::cerr << "[ERROR] Cannot find TTree '" << treeName << "' in " << inFileName << std::endl;
    fin->Close();
    return;
  }

  const Long64_t nEntries = tin->GetEntries();
  if (nEntries <= 0) {
    std::cerr << "[WARN] No entries found in tree." << std::endl;
    fin->Close();
    return;
  }

  // Enable all branches (default). You can disable/enable selectively for speed if needed.
  // tin->SetBranchStatus("*", 1);

  // Split loop
  const std::string inBase = BasenameNoExt(inFileName);

  // If you strictly want names based on run00032_* pattern, we ignore inBase and hardcode prefix.
  // But to match your examples, we derive prefix up to the first "_0_" if present.
  std::string prefix = inBase;
  {
    // Try to reduce "run00032_0_86305_aftershift" -> "run00032"
    // by taking substring before first "_0_" pattern.
    const std::string key = "_0_";
    const size_t pos = inBase.find(key);
    if (pos != std::string::npos) {
      prefix = inBase.substr(0, pos);
    } else {
      // Fallback: take substring before first '_' if exists
      const size_t p2 = inBase.find('_');
      if (p2 != std::string::npos) prefix = inBase.substr(0, p2);
    }
  }

  Long64_t start = 0;
  while (start < nEntries) {
    const Long64_t endInclusive = std::min(start + chunkSize - 1, nEntries - 1);

    // Create output file name: prefix_start_end.root
    const std::string outName =
      prefix + "_" + std::to_string(start) + "_" + std::to_string(endInclusive) + ".root";
    const std::string outPath = JoinPath(outDir, outName);

    TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
      std::cerr << "[ERROR] Cannot create output file: " << outPath << std::endl;
      if (fout) fout->Close();
      fin->Close();
      return;
    }

    // Copy tree structure and only selected entries
    // TTree::CopyTree supports selection strings, not entry ranges, so we use CopyEntries with an event list.
    // The simplest robust way is CloneTree(0) + loop Fill() (safe for big vector branches).
    TTree* tout = tin->CloneTree(0); // clone structure only
    if (!tout) {
      std::cerr << "[ERROR] Failed to clone tree structure." << std::endl;
      fout->Close();
      fin->Close();
      return;
    }

    for (Long64_t i = start; i <= endInclusive; ++i) {
      tin->GetEntry(i);
      tout->Fill();
    }

    fout->cd();
    tout->Write();
    fout->Close();

    std::cout << "[INFO] Wrote: " << outPath
              << " (entries " << start << " .. " << endInclusive << ")" << std::endl;

    start = endInclusive + 1;
  }

  fin->Close();
  std::cout << "[INFO] Done. Total entries: " << nEntries
            << ", chunkSize: " << chunkSize
            << ", outputDir: " << outDir << std::endl;
}
