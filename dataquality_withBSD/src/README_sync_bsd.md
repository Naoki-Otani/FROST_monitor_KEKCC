# BSD Sync Watcher

This repository contains a small Bash script that periodically copies new or updated files
from the BSD directory to a working directory for analysis.

The script uses `rsync` to perform incremental synchronization every 3 minutes.

---

## 1. Overview

- **Source directory (SRC_DIR)**  
  `/group/t2k/beam/exp/data/beam_summary/current/t2krun15`

- **Destination directory (DST_DIR)**  
  `/group/nu/ninja/work/otani/FROST_beamdata/e71c/bsd`

- **Behavior**
  - Every 3 minutes, the script:
    - Scans the source directory (recursively).
    - Copies **new** files to the destination.
    - Updates files in the destination if they were modified in the source.
  - Directory structure is preserved.
  - Files deleted from the source are **not** removed from the destination
    (no `--delete` option).

---

## 2. Requirements

- Bash (standard on most Linux systems)
- `rsync` installed and available in your `PATH`

You can check if `rsync` is available by running:

```bash
which rsync
rsync --version
```

If `rsync` is not installed, please ask your system administrator or install it
using your system's package manager (e.g. `apt`, `yum`, etc., depending on the environment).

---

## 3. Script File

The main script is:

- **File name:** `sync_bsd.sh`

Key configuration variables inside the script:

```bash
SRC_DIR="/group/t2k/beam/exp/data/beam_summary/current/t2krun15"
DST_DIR="/group/nu/ninja/work/otani/FROST_beamdata/e71c/bsd"
INTERVAL_SECONDS=$((3 * 60))  # Check every 3 minutes
```

You can adjust these paths or the interval if needed.

---

## 4. How to Use

### 4.1. Save the script

Save the script as `sync_bsd.sh` in a suitable directory, for example:

```bash
/home/your_user/scripts/sync_bsd.sh
```

(or any directory you prefer).

If you downloaded this file directly, just place it where you want to run it.

### 4.2. Make the script executable

Run:

```bash
chmod +x sync_bsd.sh
```

### 4.3. Run in the foreground

To run the watcher in the foreground (you see logs in the terminal):

```bash
./sync_bsd.sh
```

You can stop it with:

```bash
Ctrl + C
```

### 4.4. Run in the background

If you want it to keep running after you log out, you can use `nohup`:

```bash
nohup ./sync_bsd.sh > sync_bsd.log 2>&1 &
```

- Logs will be written to `sync_bsd.log`.
- The process will continue running in the background.

You can check running processes with, for example:

```bash
ps aux | grep sync_bsd.sh
```

And stop it with:

```bash
kill <PID>
```

(where `<PID>` is the process ID shown by `ps`).

---

## 5. Customization

### 5.1. Change the check interval

To change the check interval, edit this line in the script:

```bash
INTERVAL_SECONDS=$((3 * 60))
```

Examples:

- Every 1 minute: `INTERVAL_SECONDS=$((1 * 60))`
- Every 5 minutes: `INTERVAL_SECONDS=$((5 * 60))`

### 5.2. Make destination a "mirror" of source

By default, files that are removed from `SRC_DIR` are **not** removed from `DST_DIR`.

If you want `DST_DIR` to always exactly match `SRC_DIR` (including deletions),
you can add `--delete` to the `rsync` command in `sync_once()`:

```bash
rsync -av --delete "$SRC_DIR"/ "$DST_DIR"/
```

> ⚠️ **Warning:** With `--delete`, files that no longer exist in `SRC_DIR`
> will be deleted from `DST_DIR`. Use with care.

---

## 6. Logs

The script prints log messages with timestamps and log levels:

- `[INFO]` for normal messages
- `[ERROR]` for error messages

Example log lines:

```text
2025-11-16 12:34:56 [INFO] Starting sync cycle.
2025-11-16 12:34:57 [INFO] Sync cycle completed. Sleeping until next check.
```

When using `nohup`, these will be redirected to `sync_bsd.log` if you use:

```bash
nohup ./sync_bsd.sh > sync_bsd.log 2>&1 &
```

---

## 7. Stopping the Watcher

If the script is running in the foreground, simply press:

```bash
Ctrl + C
```

If it is running in the background, find the process with:

```bash
ps aux | grep sync_bsd.sh
```

and stop it with:

```bash
kill <PID>
```

---

## 8. Troubleshooting

- **Q: The destination directory is not created automatically.**  
  A: The script should create it if it does not exist. Check that you have
  write permissions to the parent directory of `DST_DIR`.

- **Q: I get an error about `rsync` not found.**  
  A: Install `rsync` or contact your system administrator.

- **Q: I want to use different directories.**  
  A: Edit `SRC_DIR` and `DST_DIR` at the top of the script and restart the watcher.

---

## Author
Naoki Otani
