# Data Quality Plotter (ROOT, C++)

This program monitors a directory where ROOT files (from `convertlightyield_rayraw`) are created in real time and generates several data-quality plots for the **latest run**, plus a time-history 2D plot accumulated across **all runs** in the directory with **6-hour time bins**.

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
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/runXXXXX_leading_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/runXXXXX_trailing_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/runXXXXX_leading_fromadc_hist.pdf`
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/tdc/runXXXXX_trailing_fromadc_hist.pdf`

2. **Average lightyield per channel (latest run)**  
   - For each channel (272 total) and all 8 bunches, average only values **≥ 10 p.e.**.  
   - If a channel has no value ≥ 10 p.e., its average is **0**.  
   - Histogram entries = 272.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/runXXXXX_chavg_lightyield_hist.pdf`

3. **xg–yg 2D barycenter (latest run)**  
   - Follows the previous implementation: cable mapping to x/y and weight `lightyield^XG_WEIGHT` (default 4.0).  
   - For each event, compute barycenter **per bunch** (8 times).  
   - Selection: `lightmax_x > 10` **and** `lightmax_y > 10`.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/xgyg/runXXXXX_xgyg.pdf`

4. **evnum vs unixtime graph (latest run)**

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/unixtime/runXXXXX_evnum_vs_unixtime.pdf`

5. **evnum vs spillnum graph (latest run)**

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/spillnum/runXXXXX_evnum_vs_spillnum.pdf`

6. **Average lightyield history 2D (all runs, 6-hour bins)**  
   - Accumulates per (time-bin, channel) the average of values **≥ 10 p.e.** across **all files** in the directory.  
   - Uses two on-disk caches for incremental updates:
     - `processed_files.tsv` — remember processed files and their modification times.
     - `chavg_bins.tsv` — per-bin, per-channel `(sum, count)` accumulators.
   - Fill rule (per 6-hour bin):
     - If **any channel** has data (count > 0), **fill all channels** in the bin.  
       Channels with no data get `average = 0` (to monitor potential dead channels).  
     - If **no channel** has data (all zero), **skip** that time-bin entirely.

   Output (PDF):
   - `/group/nu/ninja/work/otani/FROST_beamdata/test/dataquality/lightyield/ALL_chavg_lightyield_history_2D_6h.pdf`

---

## Requirements

- ROOT (tested with ROOT 6).  
- A POSIX-like environment (Linux).  
- The directories referenced by the program must exist or be creatable by the user.

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
  `/group/nu/ninja/work/otani/FROST_beamdata/test/rootfile_aftercalib/`
  and are named like:  
  `run00103_0_9999_lightyield.root`, `run00103_10000_19999_lightyield.root`, …

- The program determines the **latest run** by the most recently modified file and extracts its `runNNNNN` tag, aggregating all the fragments of that run, but **only** if each `*_lightyield.root` file:
  - has size **≥ 1 MB**, and
  - passes a short internal size-stability check (no growth), and
  - has not been modified in the last **60 seconds**.

  Files that do **not** satisfy these conditions are treated as “not ready” and ignored in that iteration.


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
axis->SetTimeFormat("#splitline{%Y-%m-%d}{%H:%M}");
```

If the server is pinned to UTC, you can force UTC+9:

```cpp
axis->SetTimeOffset(9*3600, "gmt");
```

---

## Notes & tips

- The program skips files that are likely still being written. A `*_lightyield.root` file is used only if it is at least **10 KB**, its size is stable in a short internal check (~500 ms), and its last modification time is at least **60 seconds** in the past.  
- If you change directory paths, edit the `PATH_*` constants in the source code.
- If labels overlap, consider increasing bottom margin or tweaking `SetNdivisions` and `SetLabelSize` on time axes.

---

## License

MIT (or your preferred internal license).
