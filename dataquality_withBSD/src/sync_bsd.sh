#!/usr/bin/env bash
set -euo pipefail

source /home/nu/notani/FROST_monitor/config/config.env

# ===== Configuration =====
SRC_DIR="$SRC_DIR_SYNCBSD"
DST_DIR="$DST_DIR_SYNCBSD"
INTERVAL_SECONDS=$((3 * 60))  # Check every 3 minutes

# ===== Helper functions =====
log() {
  # Log message with timestamp and level.
  # Example: 2025-11-16 12:34:56 [INFO] Message text
  local level="$1"
  shift
  local msg="$*"
  printf '%s [%s] %s
' "$(date '+%Y-%m-%d %H:%M:%S')" "$level" "$msg"
}

validate_directories() {
  # Ensure that source and destination directories are valid.
  if [[ ! -d "$SRC_DIR" ]]; then
    log "ERROR" "Source directory does not exist: $SRC_DIR"
    exit 1
  fi

  if [[ ! -e "$DST_DIR" ]]; then
    log "INFO" "Destination directory does not exist. Creating: $DST_DIR"
    mkdir -p "$DST_DIR"
  elif [[ ! -d "$DST_DIR" ]]; then
    log "ERROR" "Destination path is not a directory: $DST_DIR"
    exit 1
  fi
}

sync_once() {
  # Run one synchronization cycle using rsync.
  #
  # rsync options:
  #  -a : archive mode (recursive, preserve permissions, timestamps, etc.)
  #  -v : verbose
  #
  # Note:
  #  * We do NOT use --delete, so files removed from SRC_DIR
  #    will NOT be removed from DST_DIR.
  #
  # The trailing slash on SRC_DIR/ means "synchronize the contents"
  # of SRC_DIR into DST_DIR, not the directory itself as a subdirectory.
  rsync -av "$SRC_DIR"/ "$DST_DIR"/
}

cleanup() {
  # Handle termination signals (Ctrl+C, kill, etc.).
  log "INFO" "Received termination signal. Stopping watcher."
  exit 0
}

# ===== Main loop =====
log "INFO" "Starting BSD sync watcher."
log "INFO" "Source      : $SRC_DIR"
log "INFO" "Destination : $DST_DIR"
log "INFO" "Interval    : ${INTERVAL_SECONDS} seconds"

validate_directories

trap cleanup INT TERM

while true; do
  log "INFO" "Starting sync cycle."
  sync_once
  log "INFO" "Sync cycle completed. Sleeping until next check."
  sleep "$INTERVAL_SECONDS"
done
