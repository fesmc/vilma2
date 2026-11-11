#!/usr/bin/env bash
#
# run.sh — the validation suite, one stage per argument (no argument: this list).
#
#   ./run.sh build     compile the model, the benchmark drivers, the VILMA1 binary
#   ./run.sh check     submit `make check` and `make check-slow` as one job
#   ./run.sh A         block A: community benchmarks                  (13 jobs, minutes)
#   ./run.sh B         block B: lateral-viscosity benchmark, bench_visc3d (6 chains of 5 jobs; C-i ~20 h each)
#   ./run.sh ref       the LGM reference: 26 ka load held 30 kyr, lmax 170
#   ./run.sh prep      export the LGM reference, write the shared viscosity files
#   ./run.sh D         deglaciation, VILMA2, 3-D and 1Db, lmax 170
#   ./run.sh V1        the same through the VILMA1 backend
#   ./run.sh T         D's 3-D run with the toroidal field on (V7)
#   ./run.sh E         sensitivity: lmax 64/128/256 and visc3d_tol
#   ./run.sh reduce    runs/ -> data/published/
#   ./run.sh analyse   figures (docs/benchmarks/figures, figures/) and tables/
#   ./run.sh submit <stages...>   run the given stages in one SLURM job (prep, reduce, analyse)
#
# The run stages submit to SLURM and return; submit them from a login node
# (staging records the model commit with git). `ref` must have finished before
# prep, and prep before D, V1 or E:
#
#     ./run.sh build A B ref      # wait for the queue
#     ./run.sh submit prep        # wait (a few minutes)
#     ./run.sh D V1               # wait again
#     ./run.sh submit reduce analyse
#
# `reduce` skips runs that are missing or incomplete, and `analyse` reads only
# data/published/. The machine (queue, account, paths) is machines/<name>.sh,
# chosen by host name or VALIDATION_MACHINE.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
# shellcheck source=experiments/lib_stage.sh
. "$HERE/experiments/lib_stage.sh"
# shellcheck source=experiments/earth_lid.sh
. "$HERE/experiments/earth_lid.sh"
# shellcheck source=analysis/julia_env.sh
. "$HERE/analysis/julia_env.sh"

REF_NC="$HERE/data/reference/vilma2_lgm_eq.nc"
REF_RUN="$HERE/runs/blockD_ref/D-V2-3D.l170.cold.ref"
V1_EXE="$HERE/experiments/bin/vilma_v1.x"
LMAX=170

say()  { printf '\n== %s\n' "$*"; }
die()  { printf '\n!! %s\n' "$*" >&2; exit 1; }
jl()   { resolve_julia || exit 1; "$JULIA" --project="$HERE" "$@"; }   # validation/Project.toml

# ---------------------------------------------------------------- build
# The VILMA1 binary needs its own objects (-DVILMA_V1), so it is built first
# and the tree is then rebuilt for the default binary.
stage_build() {
  say "building the VILMA1-backend binary"
  [ -f "$V_VILMA1_ROOT/lib/vega_pism.a" ] || die "no VILMA1 library under $V_VILMA1_ROOT (machines/$V_MACHINE.sh)"
  ( cd "$MODEL"
    make clean
    make vilma vilma_v1=1 VILMA_V1_ROOT="$V_VILMA1_ROOT" openmp=1 -j8 )
  mkdir -p experiments/bin && cp "$MODEL/bin/vilma.x" "$V1_EXE"
  say "building the model and the benchmark drivers"
  ( cd "$MODEL"
    make clean
    make vilma openmp=1 -j8
    make openmp=1 -j8 test_benchmark_love test_benchmark_disc test_benchmark_martinec \
         test_benchmark_sle test_rotation test_rotation_sle test_benchmark_lvz \
         test_response_3d test_rotinv test_restart bench_visc3d )
}

# ---------------------------------------------------------------- check
stage_check() {
  say "submitting make check and make check-slow"
  mkdir -p runs/check
  sbatch --job-name=vilma-check --partition="$V_PARTITION" --account="$V_ACCOUNT" \
         --qos="$(v_queue 02:00:00)" --nodes=1 --ntasks=1 --cpus-per-task=32 --time=02:00:00 \
         --output="$HERE/runs/check/check.out" \
         --wrap="cd $MODEL && export OMP_NUM_THREADS=32 && make openmp=1 -k check; make openmp=1 check-slow"
}

# ---------------------------------------------------------------- blocks A, B
stage_A() { say "block A — community benchmarks"; ./experiments/run_block_A.sh; }
stage_B() { say "block B — lateral-viscosity benchmark"; ./experiments/run_block_B.sh; }

# ---------------------------------------------------------------- LGM reference
# One run at i_eq = 1: hold the 26 ka load against RTopo-2 for a fixed 30 kyr,
# from rest (equil_rate_tol = 0, so the spin-up never stops early), then one
# 100 yr step (the driver refuses a window of fewer than two slices). The state
# at the end of the hold, exported by make_lgm_equilibrium.jl, is the i_eq = 2
# reference both backends read.
stage_ref() {
  say "the LGM reference (30 kyr hold, lmax $LMAX)"
  RUNS_BLOCKD="$HERE/runs/blockD_ref" SCHEME=fe LMAX_LIST_D="$LMAX" SUFFIX=".ref" DEG1=cm \
    EQUIL=30000.0 EQTOL=0.0 PRE1D=.false. T1_D=-25900.0 \
    ./experiments/run_deglac.sh blockD-res
}

ensure_reference() {
  [ -f "$REF_NC" ] && [ "$REF_NC" -nt "$REF_RUN/job.submit" ] && return 0
  say "exporting the LGM reference"
  [ -f "$REF_RUN/job.submit" ] || die "no reference run under $REF_RUN -- run './run.sh ref' first"
  v_busy "$REF_RUN" && die "the reference run is still queued or running: $REF_RUN"
  # spinup/vilma_restart.nc is written when the hold ends. It must postdate the
  # job script: staging does not clear an earlier run's output.
  [ -f "$REF_RUN/spinup/vilma_restart.nc" ] && \
    [ "$REF_RUN/spinup/vilma_restart.nc" -nt "$REF_RUN/job.submit" ] || \
    die "the reference run has no spinup/vilma_restart.nc of its own -- see $REF_RUN/out.err"
  jl experiments/make_lgm_equilibrium.jl "$REF_RUN" "$REF_NC"
  [ -f "$REF_NC" ] || die "make_lgm_equilibrium.jl did not produce $REF_NC"
}

# ---------------------------------------------------------------- matched inputs
# The viscosity files both backends share, all derived from the clamped Bagge
# field: the laterally uniform field (1Db), and the two VILMA1 files on the
# element centres of vilma1_input_vilma2grid/densi.inp.
ensure_matched_inputs() {
  local VF=experiments/vilma1_input_vilma2grid f src want lidargs=()
  [ -f "$V_VISC_BAGGE" ] || die "missing $V_VISC_BAGGE (V_VISC_BAGGE, machines/$V_MACHINE.sh; see README.md)"
  mkdir -p data/reference
  [ -f data/reference/bagge2021_vilma2_clamp.nc ] || \
    VISC_SRC="$V_VISC_BAGGE" jl experiments/make_vilma1_visc.jl
  [ -f data/reference/bagge2021_uniform_vilma2.nc ] && [ -f data/reference/bagge2021_uniform_vilma1.nc ] || \
    jl experiments/make_uniform_visc.jl
  # The on-grid files record the lid rule they were written with and are rebuilt
  # when it is not the current one.
  if [ "$(echo "$LID_DEPTH > 0" | bc -l)" = 1 ]; then
    want=$(printf 'elements wholly above %.1f km with log10 eta > %.2f set to %.2f' \
           "$(echo "$LID_DEPTH/1000" | bc -l)" "$LID_LOG10MAX" "$LID_CAP")
    lidargs=("$LID_DEPTH" "$LID_LOG10MAX" "$LID_CAP")
  else
    want=none
  fi
  for pair in unif_vilma2grid.nc:data/reference/bagge2021_uniform_vilma1.nc \
              bagge_vilma2grid.nc:data/reference/bagge2021_vilma2_clamp.nc; do
    f=${pair%%:*}; src=${pair##*:}
    if [ ! -f "data/reference/$f" ] || \
       [[ "$(ncdump -h "data/reference/$f")" != *"lid_rule = \"$want\""* ]]; then
      jl experiments/make_vilma1_visc_ongrid.jl "$VF/densi.inp" \
        "$src" "data/reference/$f" "${lidargs[@]}"
    fi
    ln -sfn "../../data/reference/$f" "$VF/$f"
  done
  # Static VILMA1 inputs (tint.inp, load_hist.inp, vilma_grid.nc, ...): links
  # into the machine's VILMA1 input directory, for whatever is not here already.
  [ -d "$V_VILMA1_INPUT" ] || die "no VILMA1 input directory $V_VILMA1_INPUT (machines/$V_MACHINE.sh)"
  for src in "$V_VILMA1_INPUT"/*; do
    f=$(basename "$src")
    [ -e "$VF/$f" ] || [ -L "$VF/$f" ] || ln -s "$src" "$VF/$f"
  done
  ln -sfn ../../data/reference/visko_bagge.inp "$VF/visko_bagge.inp"
}

stage_prep() { ensure_reference; ensure_matched_inputs; }

# D, V1 and E only check that prep has been done: its Julia steps read and
# write files of a few hundred MB and belong in a job (./run.sh submit prep).
require_prep() {
  [ -f "$REF_NC" ] && [ "$REF_NC" -nt "$REF_RUN/job.submit" ] || \
    die "no current LGM reference $REF_NC -- run './run.sh submit prep' after 'ref' has finished"
  local f
  for f in bagge2021_uniform_vilma2.nc unif_vilma2grid.nc bagge_vilma2grid.nc; do
    [ -f "data/reference/$f" ] || die "data/reference/$f missing -- run './run.sh submit prep'"
  done
}

# ---------------------------------------------------------------- D, V1, E
# i_eq = 2 throughout: the reference is the glacial state, so there is no spin-up.
stage_D() {
  require_prep
  say "deglaciation, VILMA2 (3D and 1Db, lmax $LMAX)"
  SCHEME=fe LMAX_D=$LMAX STRUCTS="3D 1Db" DEG1=cm \
    IEQ=2 ZEQ="$REF_NC" HEQ="$REF_NC" EQUIL=0.0 PRE1D=.false. \
    ./experiments/run_deglac.sh blockD
}

stage_V1() {
  require_prep
  [ -x "$V1_EXE" ] || die "$V1_EXE missing -- run './run.sh build'"
  # runme records the model's current commit whatever binary runs; refuse a
  # VILMA1 binary older than the last commit that changed the model's sources.
  [ "$(stat -c %Y "$V1_EXE")" -ge "$(cd "$MODEL" && git log -1 --format=%ct HEAD -- src config)" ] || \
    die "$V1_EXE predates the model's last source change -- run './run.sh build'"
  ensure_matched_inputs                   # links only, once prep has run
  say "deglaciation, VILMA1 backend (3D and 1Db, matched to the D-V2 rows)"
  MODEL_TAG=D-V1 NML_TAG=dV1 EXE="$V1_EXE" LMAX_D=$LMAX STRUCTS="3D 1Db" \
    IEQ=2 ZEQ="$REF_NC" HEQ="$REF_NC" \
    ./experiments/run_deglac.sh blockDV1
}

stage_T() {
  require_prep
  say "deglaciation, VILMA2 3-D with the toroidal field (V7)"
  SCHEME=fe LMAX_D=$LMAX DEG1=cm \
    IEQ=2 ZEQ="$REF_NC" HEQ="$REF_NC" EQUIL=0.0 PRE1D=.false. \
    ./experiments/run_deglac.sh blockD-tor
}

stage_E() {
  require_prep
  say "sensitivity: truncation and visc3d_tol"
  RUNS_BLOCKD="$HERE/runs/blockE" SCHEME=fe LMAX_LIST_D="64 128 256" DEG1=cm \
    IEQ=2 ZEQ="$REF_NC" HEQ="$REF_NC" EQUIL=0.0 PRE1D=.false. \
    ./experiments/run_deglac.sh blockD-res
  SCHEME=fe LMAX_E3=$LMAX DEG1=cm \
    IEQ=2 ZEQ="$REF_NC" HEQ="$REF_NC" EQUIL=0.0 PRE1D=.false. \
    ./experiments/run_deglac.sh blockE-tol
}

# ---------------------------------------------------------------- reduce, analyse
# make_reduced.jl is the only script that opens runs/. Re-run it whenever a run
# is re-made: the figures read data/published/ only. One block: 
#     julia --project=. analysis/make_reduced.jl A        (from validation/)
stage_reduce()  { say "reduce — runs/ -> data/published/"; jl analysis/make_reduced.jl; }

# script:published input it needs; a script whose input is absent is skipped.
ANALYSIS=(run_blockA.jl:blockA.nc run_blockB.jl:bench_visc3d
          run_deglac_maps.jl:blockD run_blockD_diags.jl:blockD
          check_v2v1.jl:v2v1.nc v2v1_series_maps.jl:v2v1.nc
          run_toroidal.jl:blockD/D-V2-3D.l170.tor.nc run_blockE.jl:blockE.nc)
stage_analyse() {
  local e s failed=()
  [ -d data/published ] || die "data/published/ is missing -- run './run.sh reduce'"
  for e in "${ANALYSIS[@]}"; do
    s=${e%%:*}
    if [ ! -e "data/published/${e#*:}" ]; then
      say "analysis/$s: skipped (no data/published/${e#*:})"; continue
    fi
    say "analysis/$s"
    jl "analysis/$s" || failed+=("$s")
  done
  [ ${#failed[@]} -eq 0 ] || die "these did not complete: ${failed[*]}"
}

# ---------------------------------------------------------------- submit
stage_submit() {
  [ $# -gt 0 ] || die "submit: name the stages to run in the job"
  mkdir -p runs/jobs
  sbatch --job-name="vilma-val" --partition="$V_PARTITION" --account="$V_ACCOUNT" \
         --qos="$(v_queue 02:00:00)" --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G --time=02:00:00 \
         --output="$HERE/runs/jobs/%j.out" \
         --wrap="export JULIA_NUM_THREADS=8; $HERE/run.sh $*"
}

[ $# -gt 0 ] || { awk 'NR>1 && !/^#/{exit} NR>1' "$0" | sed 's/^# \?//'; exit 1; }
while [ $# -gt 0 ]; do
  arg=$1; shift
  case "$arg" in
    build)   stage_build ;;
    check)   stage_check ;;
    A)       stage_A ;;
    B)       stage_B ;;
    ref)     stage_ref ;;
    prep)    stage_prep ;;
    D)       stage_D ;;
    V1)      stage_V1 ;;
    T)       stage_T ;;
    E)       stage_E ;;
    reduce)  stage_reduce ;;
    analyse) stage_analyse ;;
    submit)  stage_submit "$@"; break ;;
    *)       die "unknown stage '$arg'" ;;
  esac
done
