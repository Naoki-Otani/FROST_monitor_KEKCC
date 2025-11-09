# RAYRAW Calibration Batch Tool (C++ with ROOT)

This tool provides an automated **batch calibration** for RAYRAW `.root` files using ROOT, rewritten as a **compiled C++ program** (`calibration.cpp`).

It scans the `rootfile/` directory for unprocessed ROOT files, calibrates them one by one, and outputs calibration results as **CSV** and **PDF plots**.  
Optionally, it can **continuously watch** the directory and process new files as they appear.

---

## Directory Layout (default)

```
/group/nu/ninja/work/otani/FROST_beamdata/test/
  ├── rootfile/                  # input ROOT files (e.g., run00097_0_9999.root)
  ├── calibration/
  │   ├── calibresult/           # output CSV (e.g., calib_run00097_0_9999.csv)
  │   └── plot/                  # output PDFs (e.g., calib_run00097_0_9999_RAYRAW#01_0pe.pdf)
  └── chmap/
      └── chmap_20251009.txt     # channel map (default)
```

---

## Build Instructions

Make sure ROOT is set up and available (`root-config` should work in your shell).

```bash
g++ -std=c++17 -O2 calibration.cpp -o calibration $(root-config --cflags --libs) -lSpectrum
```

> If your system is using an older GCC (e.g., GCC 8.x), you may need to add:
> ```bash
> -lstdc++fs
> ```
> at the end of the command.

---

## ▶Run (One-shot Mode)

By default, the tool scans once and processes all unprocessed `.root` files.

```bash
./calibration
```

You can specify options manually:

```bash
./calibration \
  --base /group/nu/ninja/work/otani/FROST_beamdata/test \
  --chmap chmap_20251009.txt \
  --limit 0 \
  --dry-run 0
```

### Options

| Option | Description | Default |
|:--------|:-------------|:---------|
| `--base <DIR>` | Base directory (must contain `rootfile/`, `calibration/`, `chmap/`) | `/group/nu/ninja/work/otani/FROST_beamdata/test` |
| `--chmap <FILENAME>` | Channel map file inside `chmap/` | `chmap_20251009.txt` |
| `--limit <N>` | Process at most N files (0 = no limit) | `0` |
| `--dry-run <0|1>` | If 1, only list files that would be processed | `0` |
| `--watch <SEC>` | **Watch mode:** re-scan every `<SEC>` seconds (until Ctrl+C) | 60 |
| `-h`, `--help` | Show usage help | — |

---

## Watch Mode (Continuous Monitoring)

If you want the program to keep running and automatically process newly created `.root` files, use the `--watch` option:

```bash
./calibration --watch 600
```

This example re-scans every **600 seconds (10 minutes)**.

### Features

- Detects new files added to `rootfile/`
- Processes only those **without an existing CSV**
- Skips already processed files automatically
- Gracefully stops with **Ctrl-C**
- Prints timestamps for each scan

Example output:

```
[WATCH] Start watching every 600 seconds. Press Ctrl-C to stop.
[WATCH] Scan at 2025-11-09 18:30:00
[INFO] To process (2):
  - run00098_0_9999.root [run=98, segStart=0]
  - run00098_10000_19999.root [run=98, segStart=10000]
[RUN] run00098_0_9999
...
[WATCH] Waiting 600 seconds...
```

---

## Output

Each processed `.root` file produces:

- **CSV:**  
  `/calibration/calibresult/calib_<filename>.csv`

- **PDFs:**  
  `/calibration/plot/calib_<filename>_RAYRAW#XX_0pe.pdf`  
  `/calibration/plot/calib_<filename>_RAYRAW#XX_1pe.pdf`

The CSV columns:
```
#cablenum,0pe,0pe_error,1pe,1pe_error,integralrange,badflag
```

---

## Logic Summary

1. Scan the input `rootfile/` directory.
2. Sort files by run number and segment start.
3. Check if a corresponding CSV exists — skip if already processed.
4. For each unprocessed file:
   - Perform waveform integration and peak fitting (0pe and 1pe)
   - Save plots and calibration CSV
5. (Optional) Repeat automatically at intervals if `--watch` is specified.

---

## Typical Use

```bash
# One-time processing
./calibration

# Continuous monitoring (every 5 minutes)
./calibration --watch 300

# Dry run (list only)
./calibration --dry-run 1

# Help
./calibration -h
```

---

## Troubleshooting

| Problem | Possible Fix |
|----------|---------------|
| `TSpectrum` undefined reference | Add `-lSpectrum` |
| `std::filesystem` link error | Add `-lstdc++fs` |
| `root-config: command not found` | Source ROOT environment: `source /path/to/thisroot.sh` |
| Permission denied | Check write permission to calibration directories |

---

## Notes

- ROOT runs in **batch mode** (no GUI windows).  
- Safe to interrupt and restart — already processed files are skipped.  
- You can change the `watch` interval freely; the tool will remain lightweight.

---

# RAYRAW Light-Yield Conversion Tool (C++ with ROOT)

This program (`convertlightyield.cpp`) compiles into a standalone ROOT-based tool that
**converts calibrated ADC integrals into light yield** per channel and per bunch.
It automatically scans your calibration CSV directory and, **every 60 seconds by default,**
converts any runs that have a calibration CSV but **do not yet have a converted ROOT file**.

---

## Directory Layout (defaults)

```
/group/nu/ninja/work/otani/FROST_beamdata/test/
  ├── rootfile/                      # input ROOT (e.g., run00097_0_9999.root)
  ├── rootfile_aftercalib/           # output ROOT (e.g., run00097_0_9999_lightyield.root)
  ├── chmap/
  │   └── chmap_20251009.txt         # default chmap
  └── calibration/
      ├── calibresult/               # input CSV (e.g., calib_run00097_0_9999.csv)
      └── ReferenceGain/
          └── ReferenceGain_fiberdif.csv            # reference gain CSV (configurable)
```

---

## Build

Make sure the ROOT environment is set (`root-config` works). Then:

```bash
g++ -std=c++17 -O2 convertlightyield.cpp -o convertlightyield \
  $(root-config --cflags --libs) -lSpectrum
```

> On older GCC (e.g., 8.x), you might also need `-lstdc++fs` at the end.

---

## Run

By default, the tool **watches** every **60 seconds** and converts new runs automatically:

```bash
./convertlightyield
```

One-shot processing (scan once and exit):

```bash
./convertlightyield --oneshot 1
```

Customize paths and files:

```bash
./convertlightyield \
  --base /group/nu/ninja/work/otani/FROST_beamdata/test \
  --chmap chmap_20251009.txt \
  --refgain refgain.csv \
  --watch 120     # rescan interval (seconds)
```

### Options

| Option | Description | Default |
|---|---|---|
| `--base <DIR>` | Base directory containing `rootfile/`, `rootfile_aftercalib/`, `calibration/`, `chmap/` | `/group/nu/ninja/work/otani/FROST_beamdata/test` |
| `--chmap <FILE>` | Chmap file under `chmap/` | `chmap_20251009.txt` |
| `--refgain <FILE>` | Reference gain CSV under `calibration/ReferenceGain/` | `ReferenceGain_fiberdif.csv` |
| `--watch <SEC>` | Rescan interval in seconds (0 = disabled) | `60` |
| `--oneshot <0|1>` | If `1`, scan once and exit | `0` |
| `-h`, `--help` | Show help | — |

---

## What gets converted?

The tool scans `calibration/calibresult/` for files named:

```
calib_<RUNNAME>.csv
```

If all of the following hold, it runs conversion for that `<RUNNAME>`:

- `rootfile/<RUNNAME>.root` exists
- `rootfile_aftercalib/<RUNNAME>_lightyield.root` **does not** exist (yet)

---

## Output

For each run `<RUNNAME>`:

- Output ROOT:  
  `rootfile_aftercalib/<RUNNAME>_lightyield.root`

- Output TTree (named `tree`) with branches:  
  - `cablenum[272]/I`, `rayraw[272]/I`, `rayrawch[272]/I`  
  - `lightyield[272][8]/D`, `unixtime[272]/D`  
  - `leading`, `trailing` (vector<vector<double>>)  
  - `leading_fromadc`, `trailing_fromadc` (computed from waveform threshold crossing)

The calibration CSV (`calib_<RUNNAME>.csv`) and reference gain CSV are used to compute per-channel light yield.

---

## Help

```bash
./convertlightyield -h
```

---

## Troubleshooting

- **`TSpectrum` undefined reference** → add `-lSpectrum`  
- **`std::filesystem` link error** → add `-lstdc++fs` (old GCC)  
- **`root-config` not found** → source your ROOT environment (`thisroot.sh`)  
- **Missing outputs** → Check that the input ROOT file exists under `rootfile/` and the corresponding `calib_<RUNNAME>.csv` exists under `calibration/calibresult/`
