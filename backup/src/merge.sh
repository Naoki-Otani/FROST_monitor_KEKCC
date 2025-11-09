#!/usr/bin/env bash
# merge.sh — usage: ./merge.sh 58
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <run_number>"
  echo "Example: $0 58"
  exit 1
fi

# strip leading zeros for human form (e.g., 00058 -> 58)
RUN_NO_RAW="$1"
# avoid octal: 10# forces base-10 interpretation
RUN_NO_INT=$((10#$RUN_NO_RAW))
RUN_NO_PAD=$(printf "%05d" "$RUN_NO_INT")

OUT_DIR="/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/dataaftercalib/merge"
OUT_FILE="${OUT_DIR}/run${RUN_NO_INT}merge.root"
CAEN_FILE="/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/dataaftercalib/caen/caenrun${RUN_NO_INT}_aftercalib.root"
RAYRAW_FILE="/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/dataaftercalib/rayraw/lightyield/run${RUN_NO_PAD}_lightyield.root"

# checks
command -v hadd >/dev/null 2>&1 || { echo "Error: 'hadd' not found in PATH."; exit 2; }
[[ -f "$CAEN_FILE" ]]   || { echo "Error: missing input: $CAEN_FILE"; exit 3; }
[[ -f "$RAYRAW_FILE" ]] || { echo "Error: missing input: $RAYRAW_FILE"; exit 4; }

mkdir -p "$OUT_DIR"

echo "Merging:"
echo "  CAEN  : $CAEN_FILE"
echo "  RAYRAW: $RAYRAW_FILE"
echo "    ->   $OUT_FILE"

# -f to overwrite if exists
hadd -f "$OUT_FILE" "$CAEN_FILE" "$RAYRAW_FILE"

echo "Done: $OUT_FILE"
