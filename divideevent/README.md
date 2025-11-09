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
