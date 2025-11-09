# RAYRAW Calibration Batch Tool (C++ with ROOT)

This tool converts your ROOT macro into a compiled C++ program and **automates** calibration for all `.root` files in a directory, processing **unprocessed** files from **older to newer** (by run number and segment).

## Directory Layout (defaults)
```
/group/nu/ninja/work/otani/FROST_beamdata/test/
  ├── rootfile/                  # input ROOT files (e.g., run00097_0_9999.root)
  ├── calibration/
  │   ├── calibresult/           # output CSV (e.g., calib_run00097_0_9999.csv)
  │   └── plot/                  # output PDFs (e.g., calib_run00097_0_9999_RAYRAW#01_0pe.pdf, ...)
  └── chmap/
      └── chmap_20251009.txt     # channel map (default)
```

## Build

Make sure ROOT is set up (i.e., `root-config` is on your PATH). Then:
```bash
g++ -std=c++17 -O2 calibration.cpp -o calibration $(root-config --cflags --libs) -lSpectrum
```

## Run

Basic usage with defaults (matches your current tree and chmap):
```bash
./calibration
```

Custom paths & options:
```bash
./calibration \
  --base /group/nu/ninja/work/otani/FROST_beamdata/test \
  --chmap chmap_20251009.txt \
  --limit 0 \
  --dry-run 0
```

### Options
- `--base <DIR>`: Base directory containing `rootfile/`, `calibration/`, and `chmap/`.  
  Default: `/group/nu/ninja/work/otani/FROST_beamdata/test`
- `--chmap <FILENAME>`: Chmap filename **inside** `chmap/`.  
  Default: `chmap_20251009.txt`
- `--limit <N>`: Process at most N files (0 = no limit).  
  Default: `0`
- `--dry-run <0|1>`: If `1`, only list the files that would be processed.  
  Default: `0`

## What counts as "already processed"?
A file `rootfile/runXXXXY_A_B.root` is considered processed **iff** the CSV
`calibration/calibresult/calib_runXXXXY_A_B.csv` already exists. Otherwise it will be queued.

## Ordering of processing
Files are sorted by:
1. **Run number** (ascending), then
2. **Segment start** (e.g., `0` in `0_9999`, ascending),
3. If a filename doesn't match the expected pattern, it falls back to **file mtime** (older first).

## Program structure
- Implements the previous macro logic as **compiled C++** (`calibrayraw_()`).
- Adds a **batch driver** that scans `rootfile/`, filters unprocessed files, sorts them, and runs calibration.
- Produces identical outputs (CSVs and PDFs) as your macro wrapper did.

## Example
Given these files:
```
rootfile/
  run00097_0_9999.root
  run00097_10000_19999.root
  run00098_0_9999.root
```
If only `calib_run00097_0_9999.csv` exists, then running `./calibration` will process:
```
run00097_10000_19999.root
run00098_0_9999.root
```
in that order.

## Notes
- The program runs ROOT in **batch mode** (no GUI windows).
- All comments in the source are in **English** and the code is formatted for readability.
- You can safely interrupt and re-run; already produced CSVs will be skipped.
- If you prefer stricter checks, you can extend the "already processed" criteria (e.g., also checking a couple of expected PDF outputs).

## Troubleshooting
- **`root-config: command not found`**  
  Source your ROOT environment, e.g.: `source /path/to/thisroot.sh`
- **Linker errors for TSpectrum**  
  Append `-lSpectrum` at the end of the build command.
- **Permissions**  
  Ensure you have read/write permissions to the base directory and its subfolders.
