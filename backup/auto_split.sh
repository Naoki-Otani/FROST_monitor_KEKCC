#!/usr/bin/env bash
set -Eeuo pipefail
# DO NOT exclude space from IFS; we need it for read-splitting "start end"
# If you must touch IFS, include space explicitly:
IFS=$' \n\t'

# Force C locale to avoid Unicode/locale parsing surprises
export LC_ALL=C
export LANG=C

############################################
# Config
############################################
# Incoming .dat files (runNNNNN.dat)
SRC_DIR="/group/nu/ninja/work/otani/FROST_beamdata/test/datfile"
# Output for split files
OUT_DIR="/group/nu/ninja/work/otani/FROST_beamdata/test/divided_datfile"
# Path to extractor binary
EXTRACT_CMD="/home/nu/notani/FROST_monitor/divideevent/src/extract_events"
# Wait between polls (seconds)
SLEEP_SECS=30
# Split size (fixed per requirement)
CHUNK_SIZE=10000

# Optional: directory to place probe outputs (will be auto-cleaned)
PROBE_DIR="$OUT_DIR/.probe"

# --- Logging ---
# Directory where logs are written (one file per day: YYYY-MM-DD.log)
LOG_DIR="$OUT_DIR/logs"
# If "yes", also mirror logs to console via tee. If "no", write only to file.
LOG_TO_CONSOLE="yes"

############################################
# Helpers
############################################
ts() { date "+%Y-%m-%d %H:%M:%S"; }

ensure_dirs() {
  mkdir -p "$OUT_DIR" "$PROBE_DIR" "$LOG_DIR"
}

setup_logging() {
  local logfile="$LOG_DIR/$(date +%F).log"
  if [[ "${LOG_TO_CONSOLE:-yes}" == "yes" ]]; then
    exec > >(tee -a "$logfile") 2>&1
  else
    exec >>"$logfile" 2>&1
  fi
  echo "[$(ts)] --- auto_split.sh started ---"
  echo "[$(ts)] Logging to: $logfile (mirror to console: $LOG_TO_CONSOLE)"
}

# Return basename without extension, e.g., /path/run00001.dat -> run00001
run_base_from_path() {
  local p="$1"
  local b
  b="$(basename "$p")"
  echo "${b%.*}"
}

# Find the latest .dat by mtime; echo full path or empty
latest_run_path() {
  local latest
  latest="$(ls -t "$SRC_DIR"/run*.dat 2>/dev/null | head -n1 || true)"
  [[ -n "${latest:-}" ]] && echo "$latest" || echo ""
}

# Validate unsigned integer
is_uint() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

# Parse next chunk for a run by inspecting OUT_DIR
# Echo: "<start> <end>"
next_chunk_for_run() {
  local run="$1"

  shopt -s nullglob
  local produced=("$OUT_DIR/${run}_"*_*".dat")
  shopt -u nullglob

  local max_end=-1
  if (( ${#produced[@]} > 0 )); then
    for f in "${produced[@]}"; do
      local name base start end
      name="$(basename "$f" .dat)"
      base="${name#${run}_}"               # "<start>_<end>"
      if [[ "$base" =~ ^([0-9]+)_([0-9]+)$ ]]; then
        start="${BASH_REMATCH[1]}"
        end="${BASH_REMATCH[2]}"
        if (( end > max_end )); then
          max_end="$end"
        fi
      fi
    done
  fi

  if (( max_end < 0 )); then
    printf '0 %d\n' "$((CHUNK_SIZE - 1))"
    return
  fi

  local next_start=$((max_end + 1))
  local next_end=$((next_start + CHUNK_SIZE - 1))
  printf '%d %d\n' "$next_start" "$next_end"
}

# Probe max event in a file by running extractor with an out-of-range slice.
# Prints the max event number (integer). Cleans the probe file if created.
probe_max_event() {
  local file="$1"
  local run
  run="$(run_base_from_path "$file")"

  # Use huge range that never matches any event → no writes, but prints max.
  local pstart=4294960000
  local pend=$((pstart + 1))
  local probe_out="$PROBE_DIR/${run}_${pstart}_${pend}.dat"

  local out
  if ! out="$("$EXTRACT_CMD" "$file" "$PROBE_DIR" "$pstart" "$pend" 2>&1)"; then
    echo "[$(ts)] ERROR: probe extractor failed for $file"
    echo "$out"
    return 1
  fi

  # Clean up probe file if it exists (should be empty)
  [[ -f "$probe_out" ]] && rm -f "$probe_out"

  # Parse "Max event number in file: N" robustly
  local max
  max="$(awk -F': ' '/^Max event number in file: /{print $2}' <<<"$out" | tr -cd '0-9' || true)"
  if [[ -z "${max:-}" ]]; then
    echo "[$(ts)] WARN: could not parse max event in $file; output was:"
    echo "$out"
    echo "-1"
  else
    echo "$max"
  fi
}

# Extract a chunk [start, end] for run file
extract_chunk() {
  local file="$1" start="$2" end="$3"

  # Sanitize: strip any non-digits (just in case of invisible characters)
  start="$(printf '%s' "$start" | tr -cd '0-9')"
  end="$(printf '%s' "$end" | tr -cd '0-9')"

  if ! is_uint "$start" || ! is_uint "$end"; then
    echo "[$(ts)] ERROR: non-numeric start/end after sanitize: start='$start' end='$end'"
    return 1
  fi

  printf '[%s] Extract: %s [%d-%d]\n' "$(ts)" "$(basename "$file")" "$start" "$end"
  # For debugging, print the exact command as the shell will see it:
  printf '[%s] Exec: %q %q %q %q %q\n' "$(ts)" "$EXTRACT_CMD" "$file" "$OUT_DIR" "$start" "$end"
  "$EXTRACT_CMD" "$file" "$OUT_DIR" "$start" "$end"
}

############################################
# Main loop
############################################
main() {
  ensure_dirs
  setup_logging

  local active_run="" active_path="" latest_path="" latest_run=""
  local next_start next_end max_event

  # Initialize active to the *current latest*, if any
  latest_path="$(latest_run_path)"
  if [[ -n "$latest_path" ]]; then
    active_path="$latest_path"
    active_run="$(run_base_from_path "$active_path")"
    echo "[$(ts)] Start with active run: $active_run"
  else
    echo "[$(ts)] Waiting for first run file in $SRC_DIR ..."
  fi

  while true; do
    # Detect latest by mtime
    latest_path="$(latest_run_path)"
    latest_run=""
    [[ -n "$latest_path" ]] && latest_run="$(run_base_from_path "$latest_path")"

    # If no run yet, just wait
    if [[ -z "$active_run" && -z "$latest_run" ]]; then
      sleep "$SLEEP_SECS"
      continue
    fi

    # First time we see a run
    if [[ -z "$active_run" && -n "$latest_run" ]]; then
      active_run="$latest_run"
      active_path="$latest_path"
      echo "[$(ts)] Active run set to: $active_run"
    fi

    # If latest switched to another run, drain the old run completely first
    if [[ -n "$latest_run" && "$latest_run" != "$active_run" ]]; then
      echo "[$(ts)] Detected run switch: $active_run -> $latest_run"
      # Drain remaining chunks of old run to its final max
      local old_file="$SRC_DIR/${active_run}.dat"
      if [[ -f "$old_file" ]]; then
        max_event="$(probe_max_event "$old_file")"
        if [[ "$max_event" =~ ^[0-9]+$ && $max_event -ge 0 ]]; then
          while true; do
            # IMPORTANT: ensure splitting on space
            local _line
            _line="$(next_chunk_for_run "$active_run")"
            local IFS=' '
            read -r next_start next_end <<< "$_line"
            # If already beyond the last event, break
            if (( next_start > max_event )); then
              echo "[$(ts)] Drained $active_run up to max=$max_event"
              break
            fi
            # Final chunk may be shorter than CHUNK_SIZE
            if (( next_end > max_event )); then
              next_end="$max_event"
            fi
            extract_chunk "$old_file" "$next_start" "$next_end"
          done
        else
          echo "[$(ts)] WARN: Could not determine max for $active_run; skipping drain."
        fi
      else
        echo "[$(ts)] WARN: Old run file missing: $old_file"
      fi
      # Switch to new run
      active_run="$latest_run"
      active_path="$latest_path"
      echo "[$(ts)] Switched active run to: $active_run"
    fi

    # Proceed with current active run
    local curr_file="$SRC_DIR/${active_run}.dat"
    if [[ ! -f "$curr_file" ]]; then
      echo "[$(ts)] INFO: $curr_file not found yet. Waiting..."
      sleep "$SLEEP_SECS"
      continue
    fi

    # Decide next chunk for active run (force split on space)
    local line
    line="$(next_chunk_for_run "$active_run")"
    local IFS=' '
    read -r next_start next_end <<< "$line"

    # Check if the file already contains enough events for next_end
    max_event="$(probe_max_event "$curr_file")"
    if ! [[ "$max_event" =~ ^-?[0-9]+$ ]]; then
      echo "[$(ts)] WARN: Invalid max_event result; will retry."
      sleep "$SLEEP_SECS"
      continue
    fi

    # If not enough yet, wait
    if (( max_event < next_end )); then
      echo "[$(ts)] WAIT: $(basename "$curr_file") has max=$max_event < needed end=$next_end; sleep $SLEEP_SECS"
      sleep "$SLEEP_SECS"
      continue
    fi

    # Enough data → extract the chunk
    extract_chunk "$curr_file" "$next_start" "$next_end"
    # Loop immediately to see if another chunk is also now available
  done
}

main "$@"
