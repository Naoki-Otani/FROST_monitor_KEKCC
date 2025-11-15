#!/usr/bin/env bash
set -euo pipefail

# Configuration
DATA_DIR="/group/nu/ninja/work/otani/FROST_beamdata/test/datfile"
OUT_DIR="/group/nu/ninja/work/otani/FROST_beamdata/test/latestdat_info"
INTERVAL_SEC="${INTERVAL_SEC:-300}"   # Monitoring interval in seconds (default: 5 minutes)

mkdir -p "$OUT_DIR"

# Function: Get file modification time (epoch, GNU/BSD compatible)
get_mtime_epoch() {
  local path="$1"
  if stat -c %Y "$path" >/dev/null 2>&1; then
    stat -c %Y "$path"            # GNU
  else
    stat -f %m "$path"            # BSD/macOS
  fi
}

# Function: Format timestamp to JST (YYYY/MM/DD HH:MM)
format_jst() {
  local epoch="$1"
  if TZ=Asia/Tokyo date -d "@$epoch" +'%Y/%m/%d %H:%M' >/dev/null 2>&1; then
    TZ=Asia/Tokyo date -d "@$epoch" +'%Y/%m/%d %H:%M'
  else
    TZ=Asia/Tokyo date -r "$epoch" +'%Y/%m/%d %H:%M'
  fi
}

# Function: Update JS files with the latest .dat info
update_js() {
  local latest
  if command -v find >/dev/null 2>&1 && find "$DATA_DIR" -maxdepth 1 -type f -name 'run*.dat' -printf '.' >/dev/null 2>&1; then
    latest="$(find "$DATA_DIR" -maxdepth 1 -type f -name 'run*.dat' -printf '%T@ %p\n'       | sort -nr       | head -n1       | awk '{ $1=""; sub(/^ /,""); print }')"
  else
    latest="$(ls -1t "$DATA_DIR"/run*.dat 2>/dev/null | head -n1 || true)"
  fi

  local run_name=""
  local ts_str=""

  if [[ -n "${latest:-}" && -f "$latest" ]]; then
    run_name="$(basename "$latest" .dat)"
    local mtime
    mtime="$(get_mtime_epoch "$latest")"
    ts_str="$(format_jst "$mtime")"
  else
    run_name="N/A"
    ts_str="$(TZ=Asia/Tokyo date +'%Y/%m/%d %H:%M')"
  fi

  # Write to JS files (atomic replace)
  printf 'document.write("%s");\n' "$run_name" > "$OUT_DIR/.latestdat_run.js.tmp"
  mv -f "$OUT_DIR/.latestdat_run.js.tmp" "$OUT_DIR/latestdat_run.js"

  printf 'document.write("%s");\n' "$ts_str" > "$OUT_DIR/.latestdat_date.js.tmp"
  mv -f "$OUT_DIR/.latestdat_date.js.tmp" "$OUT_DIR/latestdat_date.js"
}

# Main loop
main_loop() {
  while true; do
    update_js
    sleep "$INTERVAL_SEC"
  done
}

# Run the main loop (use INTERVAL_SEC to change interval)
main_loop
