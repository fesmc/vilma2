#!/usr/bin/env bash
#
# run_block_A.sh — stage and submit block A, the community benchmarks. Each row
# is one of the model's pass/fail test drivers, run through bin/exe.sh so that
# it also writes its profile and time-series tables into the run directory.
#
#   ./experiments/run_block_A.sh            # stage + submit everything
#   SUBMIT=0 ./experiments/run_block_A.sh   # stage only
#   ONLY="A1 A3" ./experiments/run_block_A.sh
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_stage.sh
. "$HERE/lib_stage.sh"

RUNS="${RUNS:-$VALIDATION/runs/A}"
STAGE="$VALIDATION/experiments/.stage"        # generated shims (untracked)
SUBMIT=${SUBMIT:-1}
ONLY=${ONLY:-}                                # space-separated ids to run
PAR="$MODEL/vilma.nml"                        # runme wants a namelist; the shim discards it

mkdir -p "$RUNS" "$STAGE"

# shim <name> <exe> [args...] — bin/exe.sh for one driver
shim() {
  local name=$1 exe=$2; shift 2
  local out="$STAGE/$name.sh"
  sed -e "s|@MODEL@|$MODEL|g" -e "s|@EXE@|$exe|g" -e "s|@ARGS@|$*|g" \
      "$HERE/bin/exe.sh" > "$out"
  chmod +x "$out"
  printf '%s' "$out"
}

# launch <id> <wall> <threads> <exe> [args...]
launch() {
  local id=$1 wall=$2 omp=$3 exe=$4; shift 4
  if [ -n "$ONLY" ] && [[ " $ONLY " != *" ${id%%.*} "* ]] && [[ " $ONLY " != *" $id "* ]]; then
    return 0
  fi
  echo ">>> $id  ($exe $*, wall=$wall, threads=$omp)"
  v_stage "$RUNS/$id" "$(shim "$id" "$exe" "$@")" "$PAR" "$wall" "$omp" close
  if [ "$SUBMIT" = "1" ]; then
    ( cd "$RUNS/$id" && sbatch --job-name="$id" job.submit )
  fi
}

# The drivers run from $MODEL/bin; `make check` does not build all of them.
DRIVERS=(test_benchmark_love test_benchmark_disc test_benchmark_martinec
         test_benchmark_sle test_rotation test_rotation_sle test_benchmark_lvz
         test_response_3d test_rotinv test_restart)
( cd "$MODEL" && make openmp=1 "${DRIVERS[@]}" >/dev/null )
for d in "${DRIVERS[@]}"; do
  [ -x "$MODEL/bin/$d.x" ] || { echo "run_block_A: bin/$d.x was not built" >&2; exit 1; }
done

echo "model   : $MODEL ($(cd "$MODEL" && git rev-parse --short HEAD))"
echo "machine : $V_MACHINE    runs: $RUNS    submit: $SUBMIT"
echo

# Wall times have generous headroom over the measured run times (under a minute
# each, F1 about five). Thread counts are not tuned.

# A1: Love numbers, degrees 1..256, elastic and fluid h/l/k against Spada (2011)
launch A1.love          00:20:00 16 test_benchmark_love.x

# A2: Spada (2011) disc load, u and N profiles and the centre transient
launch A2.disc          00:20:00 32 test_benchmark_disc.x

# A3: sea-level equation, Martinec (2018). Case A (pure loading) and the four
# migrating-coast cases
launch A3.caseA         00:20:00 32 test_benchmark_martinec.x
for CASE in C2 D3 E2 F1; do
  launch "A3.$CASE"     00:30:00 64 test_benchmark_sle.x "$CASE"
done

# A4: polar motion (Spada 2011 test 3/2), and the rotation coupled into the SLE
launch A4.rotation      00:20:00 16 test_rotation.x
launch A4b.rotation_sle 00:30:00 32 test_rotation_sle.x

# A5: low-viscosity zone (Weerdesteijn 2023, sect. 5.2), lmax 512
launch A5.lvz           00:30:00 64 test_benchmark_lvz.x

# A6: 3-D consistency. Uniform-eta 3-D path against the 1-D path; off-pole
# rotational invariance at lmax 128
launch A6.response_3d   00:20:00 32 test_response_3d.x
launch A6.rotinv128     00:30:00 64 test_rotinv.x 128

# A7: restart (bit-for-bit, and coarse to fine)
launch A7.restart       00:20:00 16 test_restart.x

echo
echo "done. squeue -u $USER"
