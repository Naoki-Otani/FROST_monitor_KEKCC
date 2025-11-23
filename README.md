# FROST_monitor Build Guide

This document explains where the main source files live in the `FROST_monitor` and how to build them with  `Makefile`.

## Directory Layout

The relevant parts of the directory structure are:

- `calibration/src/`
  - `calibration.cpp`
  - `convertlightyield.cpp`
- `config/`
  - `config.env`
  - `config.hpp`
- `dataquality/src/`
  - `dataqualityplot.cpp`
- `dataquality_withBSD/src/`
  - `dataqualityplot_withBSD.cpp`
  - `syn_bsd.sh`
- `divideevent/src/`
  - `extract_events.cpp`
- `monitor_latestdat/src/`
  - `update_latest_dat.sh`
- `OfflineAnalyzer/src/`
  - `convert_one_rayraw.sh`
- `OfflineAnalyzer/apps/`
  - (analysis codes rfor converting .dat files into ROOT files. 
These programs are not built by Makefile and must be compiled separately.)

## Requirements

- A C++17-capable compiler (e.g. `g++`)
- ROOT installed and `root-config` available in your `PATH`
- A suitable environment as configured via `config/config.env` and `config.hpp`

## Building the C++ Programs

From the top-level `FROST_monitor` directory, run:

```bash
make
```

This will compile the following programs and place the binaries in the same
`src` directories as their corresponding source files:

- `calibration/src/calibration`
- `calibration/src/convertlightyield`
- `dataquality/src/dataqualityplot`
- `dataquality_withBSD/src/dataqualityplot_withBSD`
- `divideevent/src/extract_events`

The build commands used are equivalent to:

```bash
g++ -std=c++17 -O2 calibration.cpp       -o calibration       $(root-config --cflags --libs) -lSpectrum
g++ -std=c++17 -O2 convertlightyield.cpp -o convertlightyield $(root-config --cflags --libs) -lSpectrum
g++ -std=c++17 -O2 dataqualityplot.cpp   -o dataqualityplot   $(root-config --cflags --libs)
g++ -std=c++17 -O2 dataqualityplot_withBSD.cpp -o dataqualityplot_withBSD $(root-config --cflags --libs)
g++ -std=c++17 -O2 extract_events.cpp    -o extract_events
```

To remove the compiled binaries, run:

```bash
make clean
```

## Configuration Files (IMPORTANT)

Before running the compiled programs on a new machine, DAQ setup, or dataset,
you **must** review and adapt the configuration files:

- `config/config.env`
- `config/config.hpp`

These files typically contain paths, detector/channel settings, thresholds,
and other environment-dependent parameters. Do **not** assume that the default
values are correct outside the original development environment.

Always check and update these files to match your local setup before using
the tools.

## OfflineAnalyzer Notes

Code under `OfflineAnalyzer/apps` is **not** built by this `Makefile`.
Those analysis programs should be compiled separately, following their own
build instructions (e.g. a dedicated `Makefile` or README within
`OfflineAnalyzer/apps`).

The shell scripts under `OfflineAnalyzer/src/` (such as
`convert_one_rayraw.sh`) and other helper scripts like
`dataquality_withBSD/src/syn_bsd.sh` and
`monitor_latestdat/src/update_latest_dat.sh` are also not compiled by this
`Makefile`; they should be executed directly as scripts once the environment
and configuration files are properly set up.
