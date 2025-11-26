# dataqualityplot_withBSD

## Overview

`dataqualityplot_withBSD` is a standalone C++17 program using ROOT to monitor
beam quality and FROST detector performance.  
It automatically processes BSD beam-summary ROOT files and FROST lightyield
ROOT files, produces POT accumulation plots, daily neutrino event-rate plots,
and generates JavaScript summary files for a web dashboard.

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

### 2. **Neutrino Event-Rate Plot**
Neutrino-like event definition:
- For each bunch (0–7), compute  
  `max(lightyield[ch][b])`  
  If this value ≥ 10, count **1 event** for that bunch.
- Events are stored per day.
- Event rate per day =  
  ```
  (# events) / (POT_recorded / 1e14)
  ```
- Output PDF:
  ```
  dataquality_withBSD/eventrate_plot/eventrate_plot.pdf
  ```

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
- `good_spill_flag`
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
- `lightyield[272][8]` (for neutrino-like events)

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
5. Neutrino events are counted per bunch (0–7)

---

## Incremental Processing

The program maintains:

- `processed_lightyield.tsv` – which LY files have been processed
- `pot_points.tsv` – incremental POT accumulation
- `neutrino_events.tsv` – timestamps of neutrino-like events

This allows fast repeated execution.

---

## Output Structure

```
dataquality_withBSD/
  pot_plot/
    pot_points.tsv
    processed_lightyield.tsv
    accumulated_pot_withBSD.pdf

  pot_info/
    period.js
    delivered_spills.js
    recorded_spills.js
    delivered_pot.js
    recorded_pot.js
    efficiency.js

  eventrate_plot/
    neutrino_events.tsv
    eventrate_plot.pdf
```

---

## Build Instructions

Make sure ROOT is configured (`root-config` available).  
Compile with:

```bash
g++ -O2 -std=c++17 dataqualityplot_withBSD.cpp -o dataqualityplot_withBSD     $(root-config --cflags --libs)
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
