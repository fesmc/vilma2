#!/usr/bin/env bash
#
# Generic runme shim for the VILMA2 benchmark drivers.
#
# runme is namelist-oriented: .runme/info.json sets par_path_as_argument=true, so
# it stages a run directory (with `data` and `input` symlinked back into the model
# tree), copies the "executable" in, and invokes `./<exe> <staged.nml>` with
# cwd = the run directory. Two things make that a poor fit for the benchmark
# drivers, and this shim fixes both:
#
#   1. Some drivers take a real CLI argument — test_benchmark_sle.x <CASE>,
#      test_rotinv.x <LMAX> — in the argv(1) slot runme has already claimed for
#      the namelist path. The shim ignores its own arguments and passes @ARGS@.
#   2. test_restart.x writes scratch restart files to ./libvilma/include/ (the
#      build's object directory), which a staged run directory does not have.
#      The shim creates it.
#
# It also turns the drivers' data dumps on. VILMA_BENCH_DUMP names the directory the
# profile/time-series tables are written to; the drivers skip the dump entirely
# when it is unset, which is what keeps `make check` unchanged. Pointing it at the
# run directory keeps each run's data next to its log.
#
# cwd stays the run directory, so the drivers' relative data/benchmarks/... paths
# resolve through the staged `data` symlink.
#
# @MODEL@, @EXE@ and @ARGS@ are substituted by validation/experiments/run_block_A.sh.
set -euo pipefail
mkdir -p libvilma/include
export VILMA_BENCH_DUMP="$PWD"
exec "@MODEL@/bin/@EXE@" @ARGS@
