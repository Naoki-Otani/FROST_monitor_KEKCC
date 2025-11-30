// eventdisplay_FROST.C
// Usage:
//   root -l
//   .L eventdisplay_FROST.C+
//   eventdisplay_FROST(5, 0);  // run number and starting "sub-event" index
//
// Display order:
//   Event 0 Bunch 0 → Event 0 Bunch 1 → ... → Event 0 Bunch 7
//   → Event 1 Bunch 7 → ...
//
// Key operation:
//   Enter = next sub-event
//   q     = quit

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>    // convert unixtime to date/time string

#include "/home/nu/notani/FROST_monitor/config/config.hpp"

// ------- Constants -------
static const int NCH_TOTAL = 272;
static const int N_BUNCH   = 8;   // lightyield[272][8]

// ------- Histogram settings for RAYRAW -------
static const int    N_BINS_X   = 132;
static const double X_MIN_X    = -660.0;
static const double X_MAX_X    =  660.0;

static const int    N_BINS_Y   = 140;
static const double X_MIN_Y    = -700.0;
static const double X_MAX_Y    =  700.0;

// Cable number → histogram bin mapping
// X: cab 332 → bin1 (left = -660), 201 → bin132 (right = +660)
static inline int cabToBinX(int cab) { return 333 - cab; }

// Y: cab 140 → bin1 (left = -700), 1 → bin140 (right = +700)
static inline int cabToBinY(int cab) { return 141 - cab; }

// ================================================================
// Main function (file name specified)
// ================================================================
void eventdisplay_FROST_(const char* rayraw_file,
                   Long64_t event_start = 0, Long64_t bunch_start = 0)
{
  // ----- Load file -----
  TFile fr(rayraw_file, "READ");
  if (fr.IsZombie()) {
    std::cerr << "Failed to open rayraw file: " << rayraw_file << "\n";
    return;
  }
  TTree* tr = (TTree*)fr.Get("tree");
  if (!tr) {
    std::cerr << "TTree 'tree' not found in rayraw file.\n";
    return;
  }

  // ----- Branch settings -----
  tr->SetBranchStatus("*", 0);

  int    cablenum[NCH_TOTAL];
  double lightyield[NCH_TOTAL][N_BUNCH];
  double unixtime[NCH_TOTAL];
  int    spillnum;
  int evnum;

  tr->SetBranchStatus("cablenum",  1);
  tr->SetBranchStatus("lightyield",1);
  tr->SetBranchStatus("unixtime",  1);
  tr->SetBranchStatus("spillnum",  1);
  tr->SetBranchStatus("evnum",  1);

  tr->SetBranchAddress("cablenum",  cablenum);
  tr->SetBranchAddress("lightyield",lightyield);
  tr->SetBranchAddress("unixtime",  unixtime);
  tr->SetBranchAddress("spillnum",  &spillnum);
  tr->SetBranchAddress("evnum",  &evnum);

  // ----- Number of events -----
  const Long64_t nEvent = tr->GetEntries();
  const Long64_t nSub   = nEvent * N_BUNCH;  // total sub-events: event × bunch

  if (event_start < 0) event_start = 0;
  if (bunch_start < 0) bunch_start = 0;
  if (event_start >= nEvent) {
    std::cerr << "event_start >= nEvent (" << event_start << " >= " << nEvent << ")\n";
    return;
  }
  if (bunch_start >= N_BUNCH) {
    std::cerr << "bunch_start >= N_BUNCH (" << bunch_start << " >= " << N_BUNCH << ")\n";
    return;
  }

  // ----- Canvas -----
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_rayraw",
                           "Event display: RAYRAW X/Y with bunch, time, spill",
                           1200, 800);
  // Upper half: X histogram, lower half: Y
  c->Divide(1, 2);

  std::string line;

  Long64_t sub_start = event_start * N_BUNCH + bunch_start;
  // Iterate over sub-events (event × bunch)
  for (Long64_t isub = sub_start; isub < nSub; ++isub) {

    Long64_t iev    = isub / N_BUNCH;    // event index
    int      ibunch = isub % N_BUNCH;    // bunch index 0..7 → displayed as 1..8

    tr->GetEntry(iev);

    // ----- Create histograms -----
    TH1D* hX = new TH1D(Form("hX_%lld_%d", iev, ibunch),
                        Form("x;"
                             "Fiber position (x) [mm];Light yield [p.e.]"),
                        N_BINS_X, X_MIN_X, X_MAX_X);

    TH1D* hY = new TH1D(Form("hY_%lld_%d", iev, ibunch),
                        Form("y;"
                             "Fiber position (y) [mm];Light yield [p.e.]"),
                        N_BINS_Y, X_MIN_Y, X_MAX_Y);

    hX->SetDirectory(nullptr);
    hY->SetDirectory(nullptr);

    // ----- Fill histograms -----
    for (int k = 0; k < NCH_TOTAL; ++k) {
      const int cab = cablenum[k];
      const double ly = lightyield[k][ibunch];
      if (!std::isfinite(ly)) continue;

      if (cab >= 201 && cab <= 332) {
        int bin = cabToBinX(cab);
        if (bin >= 1 && bin <= N_BINS_X) hX->SetBinContent(bin, ly);
      }
      else if (cab >= 1 && cab <= 140) {
        int bin = cabToBinY(cab);
        if (bin >= 1 && bin <= N_BINS_Y) hY->SetBinContent(bin, ly);
      }
    }

    // Style
    hX->SetLineColor(kBlue+1); hX->SetLineWidth(2);
    hY->SetLineColor(kRed+1);  hY->SetLineWidth(2);

    // ----- Draw histograms -----
    c->cd(1); hX->Draw("HIST");
    c->cd(2); hY->Draw("HIST");

    // ----- Draw text (bunch, unixtime, spill) -----
    c->cd(1);

    // Convert unixtime → Y/M/D H:M:S (use ch0 as representative)
    char time_str[64] = "";
    {
      time_t tt = (time_t)unixtime[0];
      struct tm* ptm = localtime(&tt);
      if (ptm) {
        strftime(time_str, sizeof(time_str), "%Y/%m/%d %H:%M:%S", ptm);
      } else {
        snprintf(time_str, sizeof(time_str), "N/A");
      }
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.12, 0.94,
                    Form("Event %d, Bunch %d, Spill %d, Time %s",
                         evnum, ibunch, spillnum, time_str));

    c->Update();
    gSystem->ProcessEvents();

    // ----- Console navigation -----
    std::cout << "Event " << evnum
              << "  Bunch " << (ibunch)
              << "  —  <Enter>=next, 'q'=quit : " << std::flush;

    std::getline(std::cin, line);
    if (line == "q") {
      // delete hX;
    //   delete hY;
      break;
    }

    // delete hX;
    // delete hY;
  }
}

// ================================================================
// Wrapper: construct filename from run number
// ================================================================
void eventdisplay_FROST(int runnum, int evnum, int evnumfinal = 0, Long64_t event_start = 0, Long64_t bunch_start = 0)
{
  static const std::string PATH_ROOT     = FrostmonConfig::OUTPUT_DIR + "/rootfile_aftercalib";
  TString filename_rayraw;

  if(evnumfinal == 0){
  	filename_rayraw = TString::Format(
      		"%s/run%05d_%d_%d_lightyield.root",
      		PATH_ROOT.c_str(), runnum, evnum, evnum + 9999);
  }else{
        filename_rayraw = TString::Format(
                "%s/run%05d_%d_%d_lightyield.root",
                PATH_ROOT.c_str(), runnum, evnum, evnumfinal);
  }
  eventdisplay_FROST_(filename_rayraw, event_start, bunch_start);
}
