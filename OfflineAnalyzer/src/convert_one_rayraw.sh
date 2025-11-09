#!/usr/bin/env bash
# convert_one_rayraw.sh — worker invoked by LSF to convert one .dat → .root
# Lives on KEKCC (shared FS). No dispatcher runs on KEKCC.
set -Eeuo pipefail
IFS=$' \n\t'
export LC_ALL=C LANG=C

if (( $# != 2 )); then
  echo "Usage: $0 <input.dat> <output.root>" >&2
  exit 2
fi

DAT="$1"
ROOT="$2"

RAYRAW_BIN="/home/nu/notani/FROST_monitor/OfflineAnalyzer/apps/k18-analyzer/bin/Rayraw"
RAYRAW_CONF="/home/nu/notani/FROST_monitor/OfflineAnalyzer/apps/param/conf/ninja_rayraw.conf"
LOCK_DIR="$(dirname "$ROOT")/.locks"

mkdir -p "$LOCK_DIR" "$(dirname "$ROOT")"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }

# Skip if already done
if [[ -f "$ROOT" || -f "${ROOT}.done" ]]; then
  echo "[$(ts)] SKIP: root exists/done: $(basename "$ROOT")"
  exit 0
fi

base="$(basename "$ROOT" .root)"
lock="$LOCK_DIR/${base}.lock"

# Acquire lock
exec {fd}>"$lock"
if ! flock -n "$fd"; then
  echo "[$(ts)] BUSY: another converter holds lock for $base"
  exit 0
fi

# Double-check after lock
if [[ -f "$ROOT" || -f "${ROOT}.done" ]]; then
  echo "[$(ts)] SKIP(after-lock): root exists/done: $(basename "$ROOT")"
  exit 0
fi

echo "[$(ts)] Convert: $(basename "$DAT") -> $(basename "$ROOT")"
echo "[$(ts)] Exec: $RAYRAW_BIN $RAYRAW_CONF $DAT $ROOT"
t0=$(date +%s)

if "$RAYRAW_BIN" "$RAYRAW_CONF" "$DAT" "$ROOT"; then
  t1=$(date +%s); elapsed=$((t1 - t0))
  echo "[$(ts)] DONE: $(basename "$ROOT") (${elapsed}s)"
  : > "${ROOT}.done"
  exit 0
else
  echo "[$(ts)] ERROR: Rayraw failed for $(basename "$DAT")"
  [[ -f "$ROOT" ]] && rm -f "$ROOT"
  exit 1
fi
