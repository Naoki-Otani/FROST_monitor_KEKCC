# Local-only Dispatcher (no dispatcher on KEKCC)

**Architecture:** Your **local PC** runs the dispatcher loop; KEKCC only runs **worker jobs**.
No dispatcher job ever runs on KEKCC.

## Files
- `local_dispatcher_submit.sh` — run on your **local PC**; periodically submits conversion jobs.
- `convert_one_rayraw.sh` — worker wrapper on KEKCC that converts one `.dat` → `.root`, using locks and `.done` markers.

## How it works
1. Local dispatcher SSHes to `kekcc` (login.cc.kek.jp via your ProxyCommand).
2. It executes in a **login shell** and explicitly sources `. /opt/lsf/conf/profile.lsf` so `bsub/bjobs` exist.
3. It lists remote `.dat` (newest first), checks **stability** (file size unchanged), and de-duplicates:
   - skip if `<base>.root` or `<base>.root.done` exists
   - skip if a job named `rayraw_<base>` is present in `bjobs`
4. It submits each job:
   ```bash
   bsub -q l -J rayraw_<base> -o <LSF_OUT>/<base>.out convert_one_rayraw.sh <dat> <root>
   ```

## Configure
Edit the top of `local_dispatcher_submit.sh`:
```bash
REMOTE_HOST="kekcc"    # ssh alias in ~/.ssh/config
DIV_DIR="/group/.../divided_datfile"
ROOT_DIR="/group/.../rootfile"
WORKER_SCRIPT="/group/.../src/convert_one_rayraw.sh"
LSF_OUT_DIR="$ROOT_DIR/lsf_out"
QUEUE="l"
JOB_NAME_PREFIX="rayraw"
INTERVAL=30
WAIT_STABLE_SECS=60
MAX_JOBS=0            # 0=unlimited
EXTRA_BSUB_OPTS=""
```

## Install on KEKCC
Place `convert_one_rayraw.sh` on the shared filesystem and make it executable:
```bash
chmod +x /group/.../src/convert_one_rayraw.sh
```

## Run on Local PC
```bash
chmod +x local_dispatcher_submit.sh
./local_dispatcher_submit.sh
```
It prints periodic status. Stop with Ctrl‑C.

## Notes
- LSF environment is loaded on every SSH call via `/opt/lsf/conf/profile.lsf`.
- Only worker jobs live on KEKCC; no dispatcher chains or resident processes.
- Concurrency throttling applies to jobs with name prefix `rayraw_`.
