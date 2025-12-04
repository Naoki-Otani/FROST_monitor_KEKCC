// config.hpp
#pragma once
#include <string>

namespace FrostmonConfig {
// ----- Path configurations -----
    const std::string OUTPUT_DIR = "/group/nu/ninja/work/otani/FROST_beamdata/e71c";

    // Default chmap file name (used when rules are not available)
    const std::string CHMAP_FILE = "chmap_20251122.txt";
    // Rule file to select chmap per run (relative to <base>/chmap/)
    const std::string CHMAP_RULE_FILE = "chmap_rules.txt";
    // Default spill chmap file name (used when rules are not available)
    const std::string SPILL_CHMAP_FILE = "chmap_spillnum20251111.txt";
    // Rule file to select spill chmap per run (relative to <base>/chmap/)
    const std::string SPILL_CHMAP_RULE_FILE = "chmap_spillnum_rules.txt";

    // Default sampling index for bunch#0 (used when rules are not available)
    static const Int_t SAMPLING_FIRST_BUNCH = 10;   // integration start index for bunch#0
    // Rule file to select SAMPLING_FIRST_BUNCH per run (relative to <base>/calibration/)
    // Each line in this file defines a (run_max, sampling_first_bunch) pair.
    const std::string SAMPLING_FIRST_BUNCH_RULE_FILE = "sampling_first_bunch_rules.txt";

    const std::string REFGAIN_CSV = "ReferenceGain_fiberdif.csv";
    const std::string LIGHTYIELD_CORR_CSV = "lightyield_correctionfactor.csv";


// ----- Calibration parameters -----
    const Int_t N_RAYRAW = 11;          // plane count (1..11)
    const Int_t N_CH_PER_PLANE = 32;    // channels per plane (0..31)
    const Int_t HIST_NBIN = 200;
    const Double_t HIST_XMIN = -100.;
    const Double_t HIST_XMAX =  300.;

    // Baseline estimation
    const Int_t BL_WIN = 15;
    const Int_t BL_START_MAX = 370;
    const Int_t BL_STEP = 10;

    // Integration windows
    const Int_t INT_LEFT  = 20; // variable window: [peak-20, peak+25]
    const Int_t INT_RIGHT = 25;

    // Fixed window for 0pe: [FIX_START, FIX_START + FIX_RANGE)
    const Int_t FIX_START = 1;
    const Int_t FIX_RANGE = INT_LEFT + INT_RIGHT;

    // TSpectrum peak search (for 1pe)
    const Int_t    PEAKS_MAX    = 2;
    const Double_t SPEC_SIGMA   = 3.0;
    const Double_t SPEC_THRESH  = 0.05;

    // ADC threshold to exclude unphysical low values
    static constexpr double ADC_MIN = 505.0;

// ----- convertlightyield parameters -----
    static const Int_t INT_RANGE = 75;
    static const Int_t ADC_THRESHOLD = 520; //threshold for leading/trailing_fromadc
    static const Double_t BUNCH_INTERVAL = 43.5;    // sampling ticks

    static const Int_t NOUT         = 272;
    static const Int_t NBUNCH       = 8;

    static const Int_t BASEREF_CALIB_CABLENUM = 311;

    static const Int_t SPILL_BIT_ADC_THRESHOLD = 600;
    static const Int_t SPILL_BIT_MIN_POINTS    = 5;

    static const Double_t PILEUP_THRESHOLD = 10.0; // threshold to set pileup flag

    static constexpr int UNDERSHOOT_ADC_THRESHOLD = 500;
    static constexpr int UNDERSHOOT_MIN_POINTS    = 8;

// ----- dataqualityplot parameters -----
    // Timing hist settings
    static const Int_t    NBINS_TDC = 2048;
    static const Double_t XMIN_TDC  = 0.0;
    static const Double_t XMAX_TDC  = 8192.0;

    static const Int_t    NBINS_ADC = 400;
    static const Double_t XMIN_ADC  = 0.0;
    static const Double_t XMAX_ADC  = 400.0;

    // Lightyield settings
    static const Int_t    NBINS_LYAVG = 50;
    static const Double_t XMIN_LYAVG  = 0.0;
    static const Double_t XMAX_LYAVG  = 200.0;

    //Lightyield per channel histogram settings
    static const Int_t NBINS_LYCH = 100;
    static const Double_t XMIN_LYCH = -10.0;
    static const Double_t XMAX_LYCH = 200.0;

    // xg–yg settings
    static const Int_t    XBINS_XY  = 132;
    static const Int_t    YBINS_XY  = 140;
    static const Double_t XMIN_XY   = -660.0;
    static const Double_t XMAX_XY   =  660.0;
    static const Double_t YMIN_XY   = -700.0;
    static const Double_t YMAX_XY   =  700.0;

    static const Double_t XG_WEIGHT = 4.0;        // weight exponent
    static const Double_t LIGHTMAX_MIN = 10.0;    // threshold for xg–yg selection and threshold for neutrino event

    // 6-hour binning for LY history
    static const Double_t BINW_SEC = 6.0 * 3600.0;

// -----dataqualityplot_withBSD parameters -----
    static const Int_t SPILL_MOD = 32768;  // 2^15
    static const Int_t MAX_TIME_DIFF = 5;  // ±5 sec allowance between BSD and LY unixtime

}
