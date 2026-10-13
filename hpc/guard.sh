# hpc/guard.sh -- sourced as the first line of every server run script (calc01 is shared). It
#  1. re-runs the calling script inside its own systemd scope, placed in the user slice pursue.slice:
#       MEMTOTAL  shared cap for ALL our runs together (the slice; default 100G) -- set on every launch, so the most
#                 recent launch's value holds for every run in the slice;
#       MEMMAX    optional cap for this run alone (default: none beyond MEMTOTAL); MEMMAX=off / MEMTOTAL=off: no scope.
#     One fresh scope per launch, gone when the run ends; nothing to close or rejoin. If the slice outgrows its cap the
#     kernel kills our workers only (finished cells are checkpointed; relaunch resumes);
#  2. makes our processes the OOM killer's first choice should the machine as a whole run out (oom_score_adj 1000);
#  3. pins BLAS to one thread (multithreaded BLAS changes msq / mid draws: srv22) -- after anything that counts cores,
#     since GNU nproc honours OMP_NUM_THREADS (p14 ran one cell at a time).
# Usage, before any cd:   source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
# Launch:   MEMTOTAL=130G CORES=16 nohup nice -n 10 bash hpc/run_x.sh ...
# Watch:    systemctl --user show pursue.slice -p MemoryCurrent -p MemoryMax
if [ -z "${PURSUE_GUARDED:-}" ] && [ "${MEMMAX:-}" != "off" ] && [ "${MEMTOTAL:-}" != "off" ]; then
  export PURSUE_GUARDED=1
  if command -v systemd-run >/dev/null 2>&1 && systemd-run --user --scope -q true >/dev/null 2>&1; then
    tot="${MEMTOTAL:-100G}"; ud="${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user"
    if [ ! -f "$ud/pursue.slice" ]; then                 # a unit file keeps the slice loadable between runs
      mkdir -p "$ud" && printf '[Unit]\nDescription=PURSUE runs on calc01 (shared memory cap)\n\n[Slice]\nMemoryMax=%s\n' "$tot" > "$ud/pursue.slice"
      systemctl --user daemon-reload 2>/dev/null || true
    fi
    systemctl --user set-property --runtime pursue.slice MemoryMax="$tot" 2>/dev/null \
      || echo ">> guard: could not set pursue.slice MemoryMax (continuing with the per-run cap only)"
    echo ">> guard: scope in pursue.slice -- all PURSUE runs together capped at $tot${MEMMAX:+, this run at $MEMMAX}"
    exec systemd-run --user --scope -q --slice=pursue.slice ${MEMMAX:+-p MemoryMax="$MEMMAX"} bash "$0" "$@"
  else
    echo ">> guard: systemd user scopes unavailable -- running without a memory cap"
  fi
fi
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
