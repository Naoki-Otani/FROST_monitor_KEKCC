// ShiftPlanes.C
// ROOT macro to shift planes (0,1,2,3,6,7,9,10) by +1 event (forward) and drop edge events.
// Compile & run:
//   root -l -q 'ShiftPlanes.C("input.root","output.root","tree")'

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <vector>
#include <string>
#include <iostream>

static inline bool IsShiftPlane(int p) {
  return (p == 0 || p == 1 || p == 2 || p == 3 || p == 6 || p == 7 || p == 9 || p == 10);
}

void ShiftPlanes(const char* inFileName  = "input.root",
                 const char* outFileName = "output.root",
                 const char* treeName    = "tree")
{
  // --- Open input ---
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
  if (nEntries < 3) {
    std::cerr << "[WARN] Need at least 3 entries to shift and drop edges. nEntries=" << nEntries << std::endl;
    fin->Close();
    return;
  }

  // --- Input branch addresses ---
  Int_t  evnum_in   = 0;
  Double_t unixtime_in = 0.0;

  std::vector<std::vector<double>>* waveform_in = nullptr;
  std::vector<double>* max_adc_in = nullptr;
  std::vector<double>* integral_in = nullptr;
  std::vector<double>* integral_ped_in = nullptr;
  std::vector<std::vector<double>>* leading_in = nullptr;
  std::vector<std::vector<double>>* trailing_in = nullptr;
  std::vector<double>* tdc_first_in = nullptr;
  std::vector<double>* de_maxadc_in = nullptr;     // "dE(MaxADC)"
  std::vector<double>* de_integral_in = nullptr;   // "dE(Integral)"
  std::vector<int>* plane_in = nullptr;

  tin->SetBranchAddress("evnum", &evnum_in);
  tin->SetBranchAddress("unixtime", &unixtime_in);

  tin->SetBranchAddress("waveform", &waveform_in);
  tin->SetBranchAddress("max_adc", &max_adc_in);
  tin->SetBranchAddress("integral", &integral_in);
  tin->SetBranchAddress("integral_ped", &integral_ped_in);
  tin->SetBranchAddress("leading", &leading_in);
  tin->SetBranchAddress("trailing", &trailing_in);
  tin->SetBranchAddress("tdc_first", &tdc_first_in);

  // Branch names contain parentheses; that's okay for ROOT, but not for C++ identifiers.
  tin->SetBranchAddress("dE(MaxADC)", &de_maxadc_in);
  tin->SetBranchAddress("dE(Integral)", &de_integral_in);

  tin->SetBranchAddress("plane", &plane_in);

  // --- Create output ---
  TFile* fout = TFile::Open(outFileName, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFileName << std::endl;
    fin->Close();
    return;
  }

  TTree* tout = new TTree(treeName, "shifted planes 4,5,8 by +1 event");

  // Output variables
  Int_t  evnum_out = 0;
  Double_t unixtime_out = 0.0;

  std::vector<std::vector<double>> waveform_out;
  std::vector<double> max_adc_out;
  std::vector<double> integral_out;
  std::vector<double> integral_ped_out;
  std::vector<std::vector<double>> leading_out;
  std::vector<std::vector<double>> trailing_out;
  std::vector<double> tdc_first_out;
  std::vector<double> de_maxadc_out;
  std::vector<double> de_integral_out;
  std::vector<int> plane_out;

  // Output branches (same names/types as input)
  tout->Branch("evnum", &evnum_out, "evnum/I");
  tout->Branch("unixtime", &unixtime_out, "unixtime/D");

  tout->Branch("waveform", &waveform_out);
  tout->Branch("max_adc", &max_adc_out);
  tout->Branch("integral", &integral_out);
  tout->Branch("integral_ped", &integral_ped_out);
  tout->Branch("leading", &leading_out);
  tout->Branch("trailing", &trailing_out);
  tout->Branch("tdc_first", &tdc_first_out);
  tout->Branch("dE(MaxADC)", &de_maxadc_out);
  tout->Branch("dE(Integral)", &de_integral_out);
  tout->Branch("plane", &plane_out);

  // --- Cache (stores only channels belonging to shift planes, from the previous input event) ---
  // We keep one-event cache for shifted channels, so we can build output with a single GetEntry per event.
  std::vector<std::vector<double>> cache_waveform;
  std::vector<double> cache_max_adc;
  std::vector<double> cache_integral;
  std::vector<double> cache_integral_ped;
  std::vector<std::vector<double>> cache_leading;
  std::vector<std::vector<double>> cache_trailing;
  std::vector<double> cache_tdc_first;
  std::vector<double> cache_de_maxadc;
  std::vector<double> cache_de_integral;

  auto resize_all = [&](size_t nch) {
    waveform_out.resize(nch);
    max_adc_out.resize(nch);
    integral_out.resize(nch);
    integral_ped_out.resize(nch);
    leading_out.resize(nch);
    trailing_out.resize(nch);
    tdc_first_out.resize(nch);
    de_maxadc_out.resize(nch);
    de_integral_out.resize(nch);

    cache_waveform.resize(nch);
    cache_max_adc.resize(nch);
    cache_integral.resize(nch);
    cache_integral_ped.resize(nch);
    cache_leading.resize(nch);
    cache_trailing.resize(nch);
    cache_tdc_first.resize(nch);
    cache_de_maxadc.resize(nch);
    cache_de_integral.resize(nch);
  };

  auto clear_cache_nonshift = [&](const std::vector<int>& pl) {
    // Keep memory usage low: clear (empty) vectors for non-shift channels.
    // Shift channels will be overwritten with real data.
    const size_t nch = pl.size();
    for (size_t ch = 0; ch < nch; ++ch) {
      if (!IsShiftPlane(pl[ch])) {
        cache_waveform[ch].clear();
        cache_leading[ch].clear();
        cache_trailing[ch].clear();
        // For scalar vectors, leaving old values is okay if we never read them for non-shift,
        // but we still set them defensively to 0.
        cache_max_adc[ch] = 0.0;
        cache_integral[ch] = 0.0;
        cache_integral_ped[ch] = 0.0;
        cache_tdc_first[ch] = 0.0;
        cache_de_maxadc[ch] = 0.0;
        cache_de_integral[ch] = 0.0;
      }
    }
  };

  auto update_cache_from_current = [&]() {
    // Copy only shift-plane channels from the *current* input entry into cache.
    if (!plane_in) return;
    const size_t nch = plane_in->size();
    if (nch == 0) return;

    // Defensive checks on input vectors
    auto safeSizeVV = [&](const std::vector<std::vector<double>>* v) -> size_t { return v ? v->size() : 0; };
    auto safeSizeV  = [&](const std::vector<double>* v) -> size_t { return v ? v->size() : 0; };

    const size_t sz_wave = safeSizeVV(waveform_in);
    const size_t sz_lead = safeSizeVV(leading_in);
    const size_t sz_trai = safeSizeVV(trailing_in);

    const size_t sz_max  = safeSizeV(max_adc_in);
    const size_t sz_int  = safeSizeV(integral_in);
    const size_t sz_ped  = safeSizeV(integral_ped_in);
    const size_t sz_tdc  = safeSizeV(tdc_first_in);
    const size_t sz_de_m = safeSizeV(de_maxadc_in);
    const size_t sz_de_i = safeSizeV(de_integral_in);

    clear_cache_nonshift(*plane_in);

    for (size_t ch = 0; ch < nch; ++ch) {
      if (!IsShiftPlane((*plane_in)[ch])) continue;

      if (waveform_in && ch < sz_wave)  cache_waveform[ch] = (*waveform_in)[ch];
      if (leading_in  && ch < sz_lead)  cache_leading[ch]  = (*leading_in)[ch];
      if (trailing_in && ch < sz_trai)  cache_trailing[ch] = (*trailing_in)[ch];

      if (max_adc_in       && ch < sz_max)  cache_max_adc[ch]       = (*max_adc_in)[ch];
      if (integral_in      && ch < sz_int)  cache_integral[ch]      = (*integral_in)[ch];
      if (integral_ped_in  && ch < sz_ped)  cache_integral_ped[ch]  = (*integral_ped_in)[ch];
      if (tdc_first_in     && ch < sz_tdc)  cache_tdc_first[ch]     = (*tdc_first_in)[ch];
      if (de_maxadc_in     && ch < sz_de_m) cache_de_maxadc[ch]     = (*de_maxadc_in)[ch];
      if (de_integral_in   && ch < sz_de_i) cache_de_integral[ch]   = (*de_integral_in)[ch];
    }
  };

  // --- Initialize cache with entry 0 (we will NOT output entry 0) ---
  tin->GetEntry(0);
  if (!plane_in) {
    std::cerr << "[ERROR] plane branch is null at entry 0." << std::endl;
    fout->Close(); fin->Close();
    return;
  }
  resize_all(plane_in->size());
  update_cache_from_current();

  // --- Main loop: output entries correspond to input 1..N-2 (drop first and last) ---
  Long64_t outCount = 0;

  for (Long64_t i = 1; i <= nEntries - 2; ++i) {
    tin->GetEntry(i);

    if (!plane_in) {
      std::cerr << "[WARN] plane is null at entry " << i << " (skipped)" << std::endl;
      continue;
    }

    const size_t nch = plane_in->size();
    if (nch == 0) {
      std::cerr << "[WARN] nch=0 at entry " << i << " (skipped)" << std::endl;
      continue;
    }

    // Ensure output containers match current event channel count
    resize_all(nch);

    // Copy "as-is" branches
    evnum_out = evnum_in;
    unixtime_out = unixtime_in;
    plane_out = *plane_in;

    // Defensive sizes
    const size_t sz_wave = (waveform_in ? waveform_in->size() : 0);
    const size_t sz_lead = (leading_in  ? leading_in->size()  : 0);
    const size_t sz_trai = (trailing_in ? trailing_in->size() : 0);

    const size_t sz_max  = (max_adc_in      ? max_adc_in->size()      : 0);
    const size_t sz_int  = (integral_in     ? integral_in->size()     : 0);
    const size_t sz_ped  = (integral_ped_in ? integral_ped_in->size() : 0);
    const size_t sz_tdc  = (tdc_first_in    ? tdc_first_in->size()    : 0);
    const size_t sz_de_m = (de_maxadc_in    ? de_maxadc_in->size()    : 0);
    const size_t sz_de_i = (de_integral_in  ? de_integral_in->size()  : 0);

    for (size_t ch = 0; ch < nch; ++ch) {
      const bool doShift = IsShiftPlane(plane_out[ch]);

      // waveform / leading / trailing (vector<vector<double>>)
      if (doShift) {
        // Use cached values from previous input event (i-1)
        waveform_out[ch] = cache_waveform[ch];
        leading_out[ch]  = cache_leading[ch];
        trailing_out[ch] = cache_trailing[ch];
      } else {
        // Use current input event (i)
        waveform_out[ch] = (waveform_in && ch < sz_wave) ? (*waveform_in)[ch] : std::vector<double>{};
        leading_out[ch]  = (leading_in  && ch < sz_lead) ? (*leading_in)[ch]  : std::vector<double>{};
        trailing_out[ch] = (trailing_in && ch < sz_trai) ? (*trailing_in)[ch] : std::vector<double>{};
      }

      // vector<double> branches
      if (doShift) {
        max_adc_out[ch]       = cache_max_adc[ch];
        integral_out[ch]      = cache_integral[ch];
        integral_ped_out[ch]  = cache_integral_ped[ch];
        tdc_first_out[ch]     = cache_tdc_first[ch];
        de_maxadc_out[ch]     = cache_de_maxadc[ch];
        de_integral_out[ch]   = cache_de_integral[ch];
      } else {
        max_adc_out[ch]       = (max_adc_in      && ch < sz_max)  ? (*max_adc_in)[ch]      : 0.0;
        integral_out[ch]      = (integral_in     && ch < sz_int)  ? (*integral_in)[ch]     : 0.0;
        integral_ped_out[ch]  = (integral_ped_in && ch < sz_ped)  ? (*integral_ped_in)[ch] : 0.0;
        tdc_first_out[ch]     = (tdc_first_in    && ch < sz_tdc)  ? (*tdc_first_in)[ch]    : 0.0;
        de_maxadc_out[ch]     = (de_maxadc_in    && ch < sz_de_m) ? (*de_maxadc_in)[ch]    : 0.0;
        de_integral_out[ch]   = (de_integral_in  && ch < sz_de_i) ? (*de_integral_in)[ch]  : 0.0;
      }
    }

    tout->Fill();
    ++outCount;

    // Update cache using the current entry i (for the next output event i+1)
    update_cache_from_current();
  }

  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();

  std::cout << "[INFO] Done. Input entries: " << nEntries
            << ", Output entries: " << outCount
            << " (kept input 1..N-2, shifted planes 4/5/8 forward by 1)" << std::endl;
}
