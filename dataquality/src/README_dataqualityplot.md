# Data Quality Plotter (ROOT, C++)

This program monitors a directory where ROOT files (from `convertlightyield`) are created in real time and generates several data-quality plots for the **latest run** (online monitoring), plus time-history plots accumulated across **all runs** (e.g. the 6-hour binned lightyield history and the calibration gain stability history).

In addition to the online monitoring of the latest run, the program also **backfills per-run PDFs for older runs**: if a run has “ready” ROOT files but some expected PDFs are missing, those PDFs are generated once, even if the run is not the latest.

It is a standalone C++ program that uses ROOT libraries (no ROOT macro).

---

## Features

1. **Timing histograms (latest run only)**  
   From `tree` branches:
   - `leading` (TDC)
   - `trailing` (TDC)
   - `leading_fromadc` (ADC)
   - `trailing_fromadc` (ADC)

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/tdc/runXXXXX_leading_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/tdc/runXXXXX_trailing_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/tdc/runXXXXX_leading_fromadc_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/tdc/runXXXXX_trailing_fromadc_hist.pdf`

2. **Average lightyield per channel (latest run)**  
   - For each channel (272 total) and all 8 bunches, average only values **≥ 10 p.e.**.  
   - If a channel has no value ≥ 10 p.e., its average is **0**.  
   - Histogram entries = 272.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/runXXXXX_chavg_lightyield_hist.pdf`

3. **Per-cablenum CSV (average lightyield and statistics, per run)**  
   - For each run, the program also writes a CSV listing (cablenum, average lightyield, #entries) computed from the same definition as feature 2:
     - Average uses only `lightyield >= 10 p.e.` over all events and all bunches.
     - `# of events over 10 p.e.` is the total number of (event, bunch) entries with `lightyield >= 10 p.e.` for that channel in that run.
   - The output is sorted by cablenum in ascending order.

   Output (CSV):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/runXXXXX_chavg_lightyield.csv`


4. **Per-channel lightyield histograms per RAYRAW plane (latest run)**  
   - For each RAYRAW plane (`1..11`) and local channel (`0..31`), build a lightyield histogram using all bunches of the latest run.  
   - One PDF per RAYRAW plane: an 8×4 pad canvas (32 channels) with **log-scale y axis**.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield_each_ch/runXXXXX_RAYRAW#YY.pdf`  
     (YY = `01`..`11`)

5. **Average lightyield and hit statistics vs position (latest run)**  
   Using the cable mapping defined in `CableToPosition(int cablenum, double& x, double& y)`, the program builds 1D profiles versus physical position along the X and Y axes.

   - **Average lightyield vs position**
     - X projection: X from **−660 mm to +660 mm** with **132 bins** (10 mm/bin).
     - Y projection: Y from **−700 mm to +700 mm** with **140 bins** (10 mm/bin).
     - For each channel:
       - Compute the average lightyield using only entries **≥ 10 p.e.** (same definition as in feature 2).
       - Fill the corresponding position bin (X or Y) with this average.
     - If multiple channels end up in the same bin, their averages are combined via a simple mean.
     - The y-axis minimum is fixed to **0**.

   - **Number of events with lightyield ≥ 10 p.e. vs position**
     - For each channel and run, the program counts how many times that channel records `lightyield ≥ 10 p.e.`.
     - For example, if `cablenum = 1` has **80** entries with `ly ≥ 10 p.e.`, the bin corresponding to that cable’s position will have height **80**.
     - These counts are accumulated in the X/Y position bins.
     - The y-axis minimum is fixed to **0**.

   - The two projections (X and Y) are shown on a single canvas, with:
     - **Top pad**: X projection
     - **Bottom pad**: Y projection

   Output (PDF):
   - Average lightyield vs position (X/Y):
     - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/runXXXXX_chavg_lightyield_profile_xy.pdf`
   - Number of entries with `ly ≥ 10 p.e.` vs position (X/Y):
     - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/runXXXXX_nevents_over10_profile_xy.pdf`

6. **xg–yg 2D barycenter (latest run)**  
   - Follows the previous implementation: cable mapping to x/y and weight `lightyield^XG_WEIGHT` (default 4.0).  
   - For each event, compute barycenter **per bunch** (8 times).  
   - Selection: `lightmax_x > 10` **and** `lightmax_y > 10`.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/xgyg/runXXXXX_xgyg.pdf`

7. **evnum vs unixtime graph (latest run)**

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/unixtime/runXXXXX_evnum_vs_unixtime.pdf`

8. **evnum vs spillnum graph (latest run)**

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/spillnum/runXXXXX_evnum_vs_spillnum.pdf`

9. **Average lightyield history 2D (all runs, 6-hour bins)**  
   - Accumulates per (time-bin, channel) the average of values **≥ 10 p.e.** across **all files** in the directory.  
   - Uses two on-disk caches for incremental updates:
     - `processed_files.tsv` — remember processed files and their modification times.
     - `chavg_bins.tsv` — per-bin, per-channel `(sum, count)` accumulators.
   - Fill rule (per 6-hour bin):
     - First, for each 6-hour bin, sum `count` over **all channels** for values `≥ 10 p.e.`  
       (this corresponds to summing the `cnt` column in `chavg_bins.tsv` over all channels
       at a given time bin).
     - If this **total count** in the bin is **greater than or equal to** `MIN_COUNTS_PER_BIN`
       (configured as `FrostmonConfig::MIN_COUNTS_PER_BIN` in the code), the bin is considered
       to have sufficient beam activity and is **filled**.
       For such bins, **all channels** are filled: channels with no data in that bin get
       `average = 0` (to monitor potential dead channels).
     - If the total count in a bin is **below** `MIN_COUNTS_PER_BIN`, that 6-hour bin is treated
       as a beam-off / low-statistics period and is **skipped entirely** (no entries are filled).
     - The X-axis time range of this plot is currently restricted to **events with unix time ≥ 2025-11-29 18:00 (local time)**.
     This cut is applied only to the **display range**; all bins are still accumulated internally.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/ALL_chavg_lightyield_history_2D_6h.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/ALL_chavg_lightyield_history_2D_6h_norm.pdf`  
     (Normalized by a reference average lightyield per cablenum; see “Reference lightyield CSV” below.)

10. **Normalized average lightyield history 2D (all runs, 6-hour bins)**  
   - If a reference CSV is available (see below), the program also produces a normalized version:
     `normalized = (average lightyield in the 6-hour bin) / (reference average lightyield for that cablenum)`.
   - The display range and the 6-hour bin fill/skip rule are the same as feature 9.
   - If the reference CSV or the cablenum mapping cannot be obtained, this plot is **skipped**.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/ALL_chavg_lightyield_history_2D_6h_norm.pdf`

11. **Calibration base-reference gain stability (all runs)**  
   The program monitors the gain stability of a fixed “base reference” calibration channel
   specified by `FrostmonConfig::BASEREF_CALIB_CABLENUM` (a `cablenum` value).

   - Input calibration CSV directory:
     - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/calibration/calibresult/`
   - Input calibration CSV naming convention:
     - `calib_run00003_0_9999.csv` (i.e. `calib_<runTag>_<seg0>_<seg1>.csv`)
   - CSV format (header):
     - `#cablenum,0pe,0pe_error,1pe,1pe_error,integralrange,badflag`
     - `badflag` is ignored for this plot.

   For each CSV file, the program reads the row with `cablenum == BASEREF_CALIB_CABLENUM` and computes:
   - Gain (y value): `gain = 1pe - 0pe`
   - Gain error (y error): `gain_err = sqrt( (1pe_error)^2 + (0pe_error)^2 )`

   Time assignment (x value and x error) uses the corresponding lightyield ROOT file
   `rootfile_aftercalib/<runTag>_<seg0>_<seg1>_lightyield.root`:
   - Let `u0` be the `unixtime` of the first event and `u1` be that of the last event
     (the code uses `unixtime[0]` in the first/last entry).
   - The point is placed at `x = (u0 + u1)/2` with `x_err = |u1 - u0|/2`.
   - **If the corresponding ROOT file does not exist or is invalid, the point is skipped** (no fallback to CSV timestamps).

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/gain/ALL_baseref_gain_history.pdf`
---

### Per-run plot behavior

- The **latest run** is always processed every iteration (online monitoring); its per-run PDFs are **continuously updated** as new data arrive.
- For **older runs**, if the program finds “ready” ROOT files and detects that one or more expected PDFs are **missing**, it generates the full set of per-run PDFs for that run **once** (backfill). Once all PDFs exist, that run is skipped in subsequent iterations.
- The **6-hour binned lightyield history 2D plots** are always accumulated over **all runs** and updated incrementally, independent of which run is latest:
  - `ALL_chavg_lightyield_history_2D_6h.pdf`
  - `ALL_chavg_lightyield_history_2D_6h_norm.pdf` (only if reference CSV + cablenum mapping are available)
- The **calibration base-reference gain stability plot** is also updated every iteration as an **all-runs** history plot. Points are filled only when both the calibration CSV and the corresponding lightyield ROOT file are available.
---

## Requirements

- ROOT (tested with ROOT 6).  
- A POSIX-like environment (Linux).  
- The directories referenced by the program must exist or be creatable by the user.
- **Reference lightyield CSV (required only for the normalized 6-hour history plot)**  
  The normalized plot `ALL_chavg_lightyield_history_2D_6h_norm.pdf` requires a reference CSV file:
  `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/ReferenceLightyield/Reference_chavg_lightyield.csv`  
  (The exact filename is configured by `FrostmonConfig::REF_LY_CSV` and the code reads it from
  `.../dataquality/lightyield/ReferenceLightyield/<REF_LY_CSV>`.)

---

## Build

```bash
g++ -O2 -std=c++17 dataqualityplot.cpp -o dataqualityplot $(root-config --cflags --libs)
```

> If your system uses a different standard (e.g. C++14), adjust the `-std=` flag accordingly.

---

## Run

```bash
./dataqualityplot            # infinite loop, default refresh interval (60 s)
./dataqualityplot -n 10      # run 10 iterations (refresh every 60 s)
./dataqualityplot -n -1 -r 5 # infinite loop, refresh every 5 seconds
```

### Command-line options

- `-n <int>` : number of iterations; negative means infinite (default: `-1`).
- `-r <ms>`  : refresh interval in milliseconds (default: `60000`).

---

## Data conventions

- Input ROOT files live in:  
  `/group/nu/ninja/work/otani/FROST_beamdata/e71c/rootfile_aftercalib/`
  and are named like:  
  `run00103_0_9999_lightyield.root`, `run00103_10000_19999_lightyield.root`, …

- The program determines the **latest run** by the most recently modified file and extracts its `runNNNNN` tag, aggregating all the fragments of that run, but **only** if each `*_lightyield.root` file is considered **“ready”**:
  - file size is at least **10 KB**, and
  - the size stays stable in a short internal check (~500 ms), and
  - the file has not been modified in the last **60 seconds**.

  Files that do **not** satisfy these conditions are treated as “not ready” and ignored in that iteration.

  For runs **other than the latest**, the same “ready file” conditions are applied. If a non-latest run has at least one ready ROOT file and some of its per-run PDFs are missing, the program generates the full set of per-run PDFs for that run (backfill) in one of the subsequent iterations.

- Reference lightyield CSV (for the normalized 6-hour history):
  - Path (default in this setup):  
    `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/ReferenceLightyield/Reference_chavg_lightyield.csv`
  - Format (per line): `#cablenum, average light yield [p.e.], # of events over 10 p.e.`
  - If the file is missing or empty, the normalized history plot is not produced.

- Calibration CSV files for gain monitoring live in:
  `/group/nu/ninja/work/otani/FROST_beamdata/e71c/calibration/calibresult/`
  and are named like:
  `calib_run00103_0_9999.csv`, `calib_run00103_10000_19999.csv`, …
  The gain stability plot uses only entries where the corresponding lightyield ROOT file exists.

 - Per-run CSV output:
   - The program also writes one per-run CSV in:
     `/group/nu/ninja/work/otani/FROST_beamdata/e71c/dataquality/lightyield/runXXXXX_chavg_lightyield.csv`
   - This CSV is generated together with the per-run PDFs (latest-run online monitoring and older-run backfill).

- Tree layout (as produced by `convertlightyield_rayraw_*`):
  - `evnum/I`
  - `cablenum[272]/I`
  - `rayraw[272]/I`
  - `rayrawch[272]/I`
  - `lightyield[272][8]/D`
  - `unixtime[272]/D`
  - `leading`, `trailing`, `leading_fromadc`, `trailing_fromadc` (all `vector<vector<double>>`)
  - `spillnum/I`

> Note: `unixtime[272]` is duplicated per channel; the code uses `unixtime[0]` as the representative timestamp per event.

---

## Time display (JST)

Time-axis plots (when applicable) are shown in **JST** by using:

```cpp
axis->SetTimeDisplay(1);
axis->SetTimeOffset(0, "local");   // relies on OS timezone (JST on Japanese servers)
axis->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
```

If the server is pinned to UTC, you can force UTC+9:

```cpp
axis->SetTimeOffset(9*3600, "gmt");
```

For the 6-hour binned lightyield history plot, the X-axis is additionally limited to times
later than `2025-11-29 18:00` (local time) when drawing; this threshold is hard-coded
via a call to `unixTimeLocal(2025, 11, 29, 18, 0, 0)` in the source code.

- To suppress artificial drops in the average lightyield during beam-off or very
  low-rate periods, each 6-hour bin is filled **only if** the total number of
  `lightyield ≥ 10 p.e.` entries over all channels in that bin is at least
  `MIN_COUNTS_PER_BIN`. This threshold is exposed as
  `FrostmonConfig::MIN_COUNTS_PER_BIN` in the configuration.
---

## Notes & tips

- The program skips files that are likely still being written. A `*_lightyield.root` file is used only if it is at least **10 KB**, its size is stable in a short internal check (~500 ms), and its last modification time is at least **60 seconds** in the past.  
- If you change directory paths, edit the `PATH_*` constants in the source code.
- The visible time range of the 6-hour binned lightyield history plot can be changed by modifying the threshold passed to `unixTimeLocal(...)` inside `BuildAndSaveLyAvgHistory2D_Binned()`.
- If labels overlap, consider increasing bottom margin or tweaking `SetNdivisions` and `SetLabelSize` on time axes.
- Bins in the 6-hour binned lightyield history plot that do not reach
  `MIN_COUNTS_PER_BIN` (i.e. too few `lightyield ≥ 10 p.e.` entries summed over
  all channels) are left empty. This is intended to avoid misleading low
  lightyield values during periods without beam.

---

## Author
Naoki Otani
