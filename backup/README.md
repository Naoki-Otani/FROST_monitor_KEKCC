# extract_events

This tool extracts a specified range of events from a RAYRAW `.dat` binary file  
and saves them into a new `.dat` file.  
It also **detects the maximum event number** in the input file and **warns** if the specified `end_event` is larger.

---

## 1. Build

Compile with `g++`:

### Linux / macOS
```bash
g++ -O2 -std=c++17 extract_events.cpp -o extract_events
```

### Windows (MinGW)
```bash
g++ -O2 -std=c++17 extract_events.cpp -o extract_events.exe
```

---

## 2. Usage

```bash
./extract_events <input.dat> <output_directory> <start_event> <end_event>
```

### Example

```bash
./extract_events /hoge/input/run00001.dat /hoge/output 10000 20000
```

- Extracts **events 10000 to 20000** from `run00001.dat`
- Output file is automatically named:

```
/hoge/output/run00001_10000_20000.dat
```

After extraction, the program prints:

```
Output file: /hoge/output/run00001_10000_20000.dat
Extraction complete.
Max event number in file: <max_event_num>
Warning: Specified end_event (<end_event>) exceeds the maximum event number in file (<max_event_num>).
```

*(Warning is shown only if `end_event` is larger than the actual maximum event number in the file.)*

---

## 3. RAYRAW Data Format

- A valid event in the `.dat` file has the following structure:

1. Starts with **`0x45564E54`** (`EVNT`) → **Magic Number**
2. If **8 words (32-bit × 8)** later is **`0x45564E54`** (`EVNT`),  
   then this block is a **valid event header**
3. **Event number** is located **2 words after the magic number** (3rd word)
4. Event data is **variable-length** and continues until the next valid event header

Example:

```
[0]  0x45564E54  ← MAGIC
[1]  0x0000191e
[2]  0x00000005  ← Event Number
[3]  0x00000002
[4]  0x00000001
[5]  0x00000000
[6]  0x00000002
[7]  0x6653f5bf
[8]  0x45564E54  ← MAGIC (8 words after)
...
[Next Event] 0x45564E54 ...
```

- Events are **in ascending order**: `0, 1, 2, 3 ...`

---

## 4. Output Specification

- Only events in the range `[start_event, end_event]` are written
- Output file is automatically named:

```
run<run#>_<start_event>_<end_event>.dat
```

Example:

```bash
./extract_events /hoge/input/run00001.dat /hoge/output 10000 20000
# => /hoge/output/run00001_10000_20000.dat
```

If `end_event` exceeds the maximum event number in the file,  
the program **outputs all available events up to the end of the file** and prints a warning.

---

## 5. Important Notes

1. **Endianness**
   - RAYRAW `.dat` files are generated on Intel/AMD PCs → **Little Endian**
   - The program reads the file as little-endian 32-bit words

2. **Maximum Event Detection**
   - The program **detects and prints the maximum event number in the input file**

3. **Large Files**
   - Streams the file word by word to avoid memory issues

4. **Early Exit for Speed**
   - Processing stops once the event number exceeds `end_event`

---

# auto_split.sh

`auto_split.sh` is a Bash script that **automates 10,000-event chunk extraction** from RAYRAW `.dat` files  
using the `extract_events` program.

It monitors your DAQ directory (e.g., `/group/nu/ninja/work/otani/FROST_beamdata/test/datfile`) where files like `run00001.dat`, `run00002.dat`, ... keep being updated.  
While a run is active, it **waits until the next 10,000-event boundary exists**, then extracts:

```
0–9999 → 10000–19999 → 20000–29999 → ...
```

**When the latest file switches to the next run** (e.g., `run00002.dat` appears),  
the script **finishes the remaining chunks of the previous run up to its last event** (e.g., `20000–27000`)  
**before** starting `run00002.dat` extraction.

---

## 1. Directory Configuration

Edit variables at the top of the script as needed:

```bash
SRC_DIR="/group/nu/ninja/work/otani/FROST_beamdata/test/datfile"
OUT_DIR="/group/nu/ninja/work/otani/FROST_beamdatatest/divided_datfile"
EXTRACT_CMD="/home/nu/notani/FROST_monitor/divideevent/src/extract_events"
SLEEP_SECS=30
CHUNK_SIZE=10000

# Optional probe directory (auto-cleaned)
PROBE_DIR="$OUT_DIR/.probe"
```

---

## 2. **Logging**

The script **dumps all output to log files** under `LOG_DIR` (default: `$OUT_DIR/logs`).  
One log file per day is created: `YYYY-MM-DD.log`.

```bash
LOG_DIR="$OUT_DIR/logs"
LOG_TO_CONSOLE="yes"   # "yes" → log & console (tee), "no" → log only
```

---

## 3. Robustness (NEW)

To avoid issues like `std::invalid_argument: stoul` (caused by unexpected non-digit characters),
the script now includes:

- `export LC_ALL=C` to force C locale
- Numeric **validation & sanitization** for `start` / `end` before calling the extractor
- More robust parsing of `"Max event number in file: N"` from extractor output
- Extra logging of the **exact command line** used to invoke the extractor

If the log shows something like:
```
ERROR: non-numeric start/end after sanitize: start='...' end='...'
```
please share that line to diagnose upstream formatting.

---

## 4. Usage

```bash
# 1) Build extractor first
g++ -O2 -std=c++17 extract_events.cpp -o extract_events

# 2) Put auto_split.sh next to (or point EXTRACT_CMD to) the extractor and make it executable
chmod +x auto_split.sh

# 3) Run normally (logs to $OUT_DIR/logs/DATE.log, also shown on console)
./auto_split.sh

# (Optional) Run in background; still logs to files
nohup ./auto_split.sh >/dev/null 2>&1 &
```

---

## 5. Notes

- Output files are named:  
  `$OUT_DIR/runNNNNN_<start>_<end>.dat`
- The script is **idempotent**: it derives progress from existing output files; already-created chunks are not re-created.
- The “probe” to read `Max event number` **does not** produce partial chunks; any temp file is deleted.
- If a run file temporarily disappears (e.g., move/rotate), the script waits until it reappears or until a newer run is present.
- Logs rotate **daily** by filename (no size-based rotation).

---

## 6. License

This script is free to use and modify for your experiments.

# extract_events

(omitted identical sections for brevity)

# auto_split.sh — Logging & Robustness Notes

- The script **relies on splitting `"start end"` by a space**.  
  Therefore, **do not remove space from `IFS`**. The script sets
  ```bash
  IFS=$' \n\t'
  ```
  to include a **space** explicitly, avoiding the bug where `"0 9999"`
  is read as a single token (which led to `start='09999' end=''`).
- Arguments to `extract_events` are **sanitized to digits** and validated.
- Locale is forced to **C** to avoid non-ASCII surprises.
