# dataqualityplot_withBSD

## Overview

`dataqualityplot_withBSD` is a standalone C++17 program using ROOT to monitor
beam quality and FROST detector performance.  
It automatically processes BSD beam-summary ROOT files and FROST lightyield
ROOT files, produces POT accumulation plots, daily event-rate plots (per bunch
and per spill), and generates JavaScript summary files for a web dashboard.

---

## Main Features

### 1. **Accumulated POT Plot**
- Delivered POT (from BSD): **black curve**
- Recorded POT (matched to FROST): **red curve**
- Matching uses spill number + ±5s time tolerance.
- Output PDF:
  ```
  dataquality_withBSD/pot_plot/accumulated_pot_withBSD.pdf
  ```

### 2. **Event-Rate Plots**

Common event definition:
- For each bunch (0–7), compute  
  `max(lightyield[ch][b])`  
  If this value ≥ 10, that bunch is considered to have **one event**.

Two types of daily event-rate histograms are produced:

1. **Per-bunch counting (default)**
   - Each bunch satisfying the condition counts as one event.
   - A single spill can contribute up to 8 events (one per bunch).
   - Events are accumulated per day.
   - Event rate per day:
     ```text
     (# events) / (POT_recorded / 1e15)
     ```
   - Output PDF:
     ```text
     dataquality_withBSD/eventrate_plot/eventrate_plot.pdf
     ```

2. **Per-spill counting (max 1 event per spill)**
   - If at least one bunch in a spill satisfies the condition, that spill
     contributes 1 event (even if multiple bunches fired).
   - Events are accumulated per day.
   - Event rate per day:
     ```text
     (# events) / (POT_recorded / 1e15)
     ```
   - Output PDF:
     ```text
     dataquality_withBSD/eventrate_plot/eventrate_plot_spill.pdf
     ```

Both histograms are fitted with a constant function `y = p0`.  
The fitted `χ²/ndf` and `p0 ± error` are shown in a statistics box in the
top-right corner of each plot.

### 3. **JavaScript Summary Files**
Automatically written to:
```
dataquality_withBSD/pot_info/
```
Files include:
- `period.js` – data-taking period (start–end)
- `delivered_spills.js`
- `recorded_spills.js`
- `delivered_pot.js` (scientific notation)
- `recorded_pot.js`
- `efficiency.js` (Recorded / Delivered [%])

These are used to feed an HTML dashboard.

---

## Input Data

### BSD Files
Directory:
```
/group/nu/ninja/work/otani/FROST_beamdata/e71c/bsd/
```
Tree name: **bsd**  
Used branches:
- `spillnum`
- `trg_sec[2]` (Unix time)
- `good_spill_flag` (`spill_flag` for QSD)
- `ct_np[4][0]` (POT)

### Lightyield Files
Directory:
```
/group/nu/ninja/work/otani/FROST_beamdata/e71c/rootfile_aftercalib/
```
Tree name: **tree**  
Used branches:
- `spillnum` (15-bit wrapping)
- `unixtime[0]`
- `lightyield[272][8]` (for event definition)

Only stable & fully written files are processed.

---

## Matching Logic

Every lightyield event is matched to a BSD spill via:

1. `key = spillnum & 0x7FFF`  
2. Find BSD spills with identical 15-bit spill number  
3. Choose:
   ```
   |bsd.trg_sec – round(unixtime[0])| ≤ 5 seconds
   ```
   with minimal |Δt|
4. POT from a BSD spill is counted **once per LY file**
5. Events are counted per bunch (0–7) for the **per-bunch** rate.  
   For the **per-spill** rate, a spill contributes 1 event if any acquired
   bunch satisfies the event condition.

---

## Incremental Processing

The program maintains:

- `pot_points.tsv` – POT increments (time, POT) for matched BSD spills,
  used to build the accumulated POT curve and daily recorded POT.
- `eventrate.tsv` – timestamps (Unix time) of events counted **per bunch**.
- `eventrate_spill.tsv` – timestamps (Unix time) of events counted
  **per spill** (max 1 event per spill).

These cache files are rebuilt from all available, stable lightyield files at
each iteration, so repeated executions always reflect the full dataset.

---

## Output Structure

```
dataquality_withBSD/
  pot_plot/
    pot_points.tsv
    accumulated_pot_withBSD.pdf

  pot_info/
    period.js
    delivered_spills.js
    recorded_spills.js
    delivered_pot.js
    recorded_pot.js
    efficiency.js

  eventrate_plot/
    eventrate.tsv
    eventrate_spill.tsv
    eventrate_plot.pdf
    eventrate_plot_spill.pdf
```

---

## Build Instructions

Make sure ROOT is configured (`root-config` available).  
Compile with:

```bash
g++ -O2 -std=c++17 dataqualityplot_withBSD.cpp -o dataqualityplot_withBSD $(root-config --cflags --libs)
```

---

## Running the Program

### Run continuously (default)
```bash
./dataqualityplot_withBSD
```
Refreshes every 10 minutes.

### Run once
```bash
./dataqualityplot_withBSD -n 1
```

### Continuous run with custom refresh interval (e.g., 5 sec)
```bash
./dataqualityplot_withBSD -n -1 -r 5000
```

### Help
```bash
./dataqualityplot_withBSD -h
```

---

## Author
Naoki Otani
