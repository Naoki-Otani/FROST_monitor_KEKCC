# RAYRAW Calibration Batch Tool (C++ with ROOT)

This tool provides an automated **batch calibration** for RAYRAW `.root` files using ROOT, rewritten as a **compiled C++ program** (`calibration.cpp`).

It scans the `rootfile/` directory for unprocessed ROOT files, calibrates them one by one, and outputs calibration results as **CSV** and **PDF plots**.  
Optionally, it can **continuously watch** the directory and process new files as they appear.

> **Safety against half-written ROOT files**  
> Calibration now processes a ROOT file **only if all of the following are true**:  
> - The file size is **at least 10 KB**  
> - The file size has not changed for **at least 60 seconds**  
> These conditions prevent processing files that are still being written by the DAQ or converter.

---

## Directory Layout (default)

```
/group/nu/ninja/work/otani/FROST_beamdata/e71c/
  ├── rootfile/                  # input ROOT files (e.g., run00097_0_9999.root)
  ├── calibration/
  │   ├── calibresult/           # output CSV (e.g., calib_run00097_0_9999.csv)
  │   └── plot/                  # output PDFs (e.g., calib_run00097_0_9999_RAYRAW#01_0pe.pdf)
  └── chmap/
      └── chmap_20251009.txt     # channel map (example)
      └── chmap_rules.txt        # optional per-run chmap rules

If `chmap_rules.txt` exists, it is used to select the chmap file per run (unless `--chmap` is given).
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
  --base /group/nu/ninja/work/otani/FROST_beamdata/e71c \
  --chmap chmap_20251009.txt \
  --limit 0 \
  --dry-run 0
```

### Options

| Option | Description | Default |
|:--------|:-------------|:---------|
| `--base <DIR>` | Base directory (must contain `rootfile/`, `calibration/`, `chmap/`) | `/group/nu/ninja/work/otani/FROST_beamdata/e71c` |
| `--chmap <FILENAME>` | Channel map file inside `chmap/`. If given, it is used for **all runs** and per-run `chmap_rules.txt` is ignored. | config default (e.g. `chmap_20251122.txt`) |
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
- Processes only files **≥10 KB** and **unchanged for ≥60 seconds**
- Processes only those **without an existing CSV**
- Skips already processed files automatically
- Gracefully stops with **Ctrl-C**
- Prints timestamps for each scan

Example output:

```
[WATCH] Start watching every 600 seconds. Press Ctrl-C to stop.
[WATCH] Scan at 2025-11-09 18:30:00
[INFO] To process (2) (only files >= 10KB and stable for >= 60s are considered):
  - run00098_0_9999.root [run=98, segStart=0]
  - run00098_10000_19999.root [run=98, segStart=10000]
  (size=123.4MB etc.)
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
2. Ignore files smaller than **10 KB** or whose size changed within the last **60 seconds**.
3. Sort remaining files by run number and segment start.
4. Check if a corresponding CSV exists — skip if already processed.
5. For each unprocessed & ready file:
   - Perform waveform integration and peak fitting (0pe and 1pe)
   - Save plots and calibration CSV
6. (Optional) Repeat automatically at intervals if `--watch` is specified.

### Per-run chmap selection

By default, `calibration` can select the chmap file **per run** using the optional
`chmap/chmap_rules.txt` file.

The format is:

```text
# run_max  chmap_filename
7     chmap_20251122.txt
9999  chmap_20251130.txt
```

For each input file whose run number is `run`, the program:

1. Parses the integer run number (e.g. `run00005_0_9999.root` → `5`).
2. Finds the first row with `run <= run_max`.
3. Uses that `chmap_filename` for calibration.

If `chmap_rules.txt` is missing or no rule matches, the **default** chmap from `config.hpp`
is used.

If you pass `--chmap <FILENAME>`, that file is used for **all runs**, and `chmap_rules.txt` is ignored.

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

> **Safety against half-written calibration CSV files**  
> `convertlightyield` now converts runs **only when the calibration CSV satisfies both**:
> - CSV size is **at least 1 KB**
> - CSV's size has not changed for **at least 60 seconds**  
>  
> This prevents conversion using incomplete or still-being-written calibration files.

---

## Directory Layout (defaults)

```
/group/nu/ninja/work/otani/FROST_beamdata/e71c/
  ├── rootfile/                      # input ROOT (e.g., run00097_0_9999.root)
  ├── rootfile_aftercalib/           # output ROOT (e.g., run00097_0_9999_lightyield.root)
  ├── chmap/
  │   ├── chmap_20251122.txt         # default RAYRAW chmap (example)
  │   ├── chmap_spillnum20251111.txt # spillnum bit → (RAYRAW, ch) map (default)
  │   ├── chmap_rules.txt            # optional per-run RAYRAW chmap rules
  │   └── chmap_spillnum_rules.txt       # optional per-run spillnum chmap rules
  └── calibration/
      ├── calibresult/               # input CSV (e.g., calib_run00097_0_9999.csv)
      ├── lightyield_correctionfactor/
      │   └── lightyield_correctionfactor.csv       # per-cable light-yield correction factors
      ├── sampling_first_bunch/
      │   └── sampling_first_bunch_rules.txt        # optional per-run sampling-start rules
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
  --base /group/nu/ninja/work/otani/FROST_beamdata/e71c \
  --chmap chmap_20251009.txt \
  --refgain refgain.csv \
  --lycorr lightyield_correctionfactor.csv \
  --spillchmap chmap_spillnum20251111.txt \
  --watch 120     # rescan interval (seconds)
```

### Options

| Option | Description | Default |
|---|---|---|
| `--base <DIR>` | Base directory containing `rootfile/`, `rootfile_aftercalib/`, `calibration/`, `chmap/` | `/group/nu/ninja/work/otani/FROST_beamdata/e71c` |
| `--chmap <FILE>` | RAYRAW chmap file under `chmap/`. If given, it is used for **all runs** and per-run `chmap_rules.txt` is ignored. | config default (e.g. `chmap_20251122.txt`) |
| `--refgain <FILE>` | Reference gain CSV under `calibration/ReferenceGain/` | `ReferenceGain_fiberdif.csv` |
| `--lycorr <FILE>` | Light-yield correction CSV under `calibration/lightyield_correctionfactor/` | `lightyield_correctionfactor.csv` |
| `--spillchmap <FILE>` | Spill-number bit mapping file under `chmap/`. If given, it is used for **all runs** and per-run `chmap_spillnum_rules.txt` is ignored. | `chmap_spillnum20251111.txt` |
| `--watch <SEC>` | Rescan interval in seconds (0 = disabled) | `60` |
| `--oneshot <0|1>` | If `1`, scan once and exit | `0` |
| `-h`, `--help` | Show help | — |

---

## What gets converted?

The tool scans `calibration/calibresult/` for files named:

```
calib_<RUNNAME>.csv
```

A run `<RUNNAME>` is **eligible for conversion** only if **all** of the following 3 conditions are satisfied:

1. **The calibration CSV is ready**  
   - CSV file size is **at least 1 KB**, **and**  
   - the size has not changed for **at least 60 seconds**  
   → If too small or recently updated, it is assumed incomplete and skipped.

2. **The corresponding input ROOT exists**  
   ```
   rootfile/<RUNNAME>.root
   ```

3. **The output (light-yield) ROOT does NOT already exist**  
   ```
   rootfile_aftercalib/<RUNNAME>_lightyield.root
   ```

If a CSV is not ready, the program reports something like:

```
[SKIP] CSV not ready (too small or recently updated): calib_run00097_0_9999.csv (size=512 bytes)
```

## Per-run configuration (chmap, spillchmap, sampling start index)

The conversion tool selects several configuration parameters **per run** based on
simple rule files.

### RAYRAW chmap rules (`chmap/chmap_rules.txt`)

If present, this file selects the RAYRAW chmap per run.  
Format:
```text
# run_max  chmap_filename
7     chmap_20251122.txt
9999  chmap_20251130.txt
```

For a run with integer run number `run` (e.g. `run00005_0_9999` → `5`),
the program finds the first line with `run <= run_max` and uses that `chmap_filename`.

If the file is missing or no rule matches, the default chmap from `config.hpp` is used.

If `--chmap` is specified on the command line, that file is used for **all runs** and
`chmap_rules.txt` is ignored.

### Spillnum chmap rules (`chmap/chmap_spillnum_rules.txt`)

Similarly, if `chmap_spillnum_rules.txt` exists, it selects the spill chmap file per run.  
Format:

```text
# run_max  spill_chmap_filename
5      chmap_spillnum20251111.txt
20     chmap_spillnum20251201.txt
999999 chmap_spillnum20260101.txt
```

If `--spillchmap` is specified, that file is used for all runs and the rules file is ignored.

### Sampling start index rules (`calibration/sampling_first_bunch/sampling_first_bunch_rules.txt`)

The integration for each bunch starts at a sample index called `sampling_first_bunch`.
By default this is a single constant value from `config.hpp` (e.g. `10`),

```cpp
const Int_t SAMPLING_FIRST_BUNCH = 10;
```

If the rule file exists, it can override this value per run.

Format:

```text
# run_max  sampling_first_bunch
5      10
20     12
999999 15
```

For each run, the program parses the integer run number from `<RUNNAME>` (e.g. `run00006_0_9999` → `6`),
selects the RAYRAW chmap, spill chmap, and `sampling_first_bunch` based on these rules (or defaults),
and runs the conversion using those settings.


---

## Output

For each run `<RUNNAME>`:

- Output ROOT:  
  `rootfile_aftercalib/<RUNNAME>_lightyield.root`

- Output TTree (named `tree`) with branches:  
  - `cablenum[272]/I`, `rayraw[272]/I`, `rayrawch[272]/I`  
  - `lightyield[272][8]/D`, `unixtime[272]/D`  
  - `spillnum/I`
  - `leading`, `trailing` (vector<vector<double>>)  
  - `leading_fromadc`, `trailing_fromadc` (computed from waveform threshold crossing)

The calibration CSV (`calib_<RUNNAME>.csv`) and reference gain CSV are used to compute **raw** per-channel light yield,
which is then multiplied by a **per-cable correction factor** from `lightyield_correctionfactor.csv`.

### Light-yield correction factor

The file `lightyield_correctionfactor.csv` has the format:

```text
#cablenum,correction_factor
1,1.04164
2,1.15779
...
```

For each `cablenum`:

- First, the raw light yield is computed from the calibrated ADC integral and reference gain information.
- Then, the stored `lightyield` is:

  \[
    \text{lightyield}_\text{out} = \text{lightyield}_\text{raw} \times \text{correction\_factor}(\text{cablenum})
  \]

### Spillnum reconstruction
The file specified by `--spillchmap` contains:

```
#bit RAYRAW ch
1 12 9
2 12 10
...
15 12 23
```

For each listed `(RAYRAW, ch)`:

- Look at the waveform of that channel  
- If **≥ 5 samples ≥ 600 ADC**, the bit is **1**  
- Otherwise, bit is **0**

Bits are interpreted as:
- bit1 → 2⁰  
- bit2 → 2¹  
- …  
- bit15 → 2¹⁴

All bits are combined into an integer and stored in `spillnum`.

Which `(RAYRAW, ch)` channels are used for each spill bit is defined by the
spill chmap selected for that run (either from `--spillchmap` or `chmap_spillnum_rules.txt`).
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
