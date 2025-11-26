# Latest .dat File Monitor Script

This script monitors the directory containing `.dat` files in real time and updates two JavaScript files:
- `latestdat_run.js` — contains the latest run name (e.g., `document.write("run00106");`)
- `latestdat_date.js` — contains the last modification date/time of that file (e.g., `document.write("2025/11/12 23:55");`)

Both files are continuously updated (default: every 5 minutes).

---

## Directory Structure

```
/group/nu/ninja/work/otani/FROST_beamdata/e71c/datfile         # Source .dat files
/group/nu/ninja/work/otani/FROST_beamdata/e71c/latestdat_info  # Output JS files
```

---

## Installation

1. **Save the script:**
   ```bash
   cat > ~/update_latest_dat.sh <<'EOS'
   # (Paste the full script here)
   EOS
   chmod +x ~/update_latest_dat.sh
   ```

2. **Run continuously (every 5 minutes):**
   ```bash
   nohup ~/update_latest_dat.sh >/tmp/update_latest_dat.log 2>&1 &
   ```

3. **Change interval (e.g., every 60 seconds):**
   ```bash
   INTERVAL_SEC=60 nohup ~/update_latest_dat.sh >/tmp/update_latest_dat.log 2>&1 &
   ```

---

## Alternative: Cron Job (every 5 minutes)

If you prefer to use `cron` instead of a continuous loop, modify the script to run `update_js` once (not in a loop), then add this to your crontab:

```bash
crontab -e
# Add this line:
*/5 * * * * /bin/bash -lc '~/update_latest_dat.sh >> /tmp/update_latest_dat.log 2>&1'
```

---

## Notes

- The script automatically detects the latest `.dat` file based on modification time.
- It handles both GNU and BSD/macOS versions of `stat` and `date` commands.
- Output files are written atomically to avoid race conditions.
- Timezone is fixed to **Asia/Tokyo**.
